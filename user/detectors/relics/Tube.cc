#include "Tube.hh"
#include "PandaXSensitiveDetector.hh"

#include <G4Material.hh>
#include <G4LogicalVolume.hh>
#include <G4Tubs.hh>

DetectorRegister<Tube, std::string, BambooParameters> Tube::reg("Tube");

Tube::Tube (const std::string &n, const BambooParameters &pars)
  : Monoblock(n, pars) {
    G4cout << "create detector Tube..." << G4endl;
}

bool Tube::constructMainLV(const BambooParameters &) {
    using namespace CLHEP;
    auto soliname = parameters.getParameter("soliname");
    auto logivol = parameters.getParameter("logivol");
    auto height = parameters.evaluateParameter("height");
    auto outer_radius = parameters.evaluateParameter("outer_radius");
    auto inner_radius = parameters.evaluateParameter("inner_radius");
    auto material = parameters.getParameter("material");
    
    // 设置默认值
    if (height == 0) {
        height = 2 * m;
    }
    if (outer_radius == 0) {
        outer_radius = 1 * m;
    }
    if (inner_radius == 0) {
        inner_radius = 0 * m;  // 实心圆柱
    }
    
    // 创建管状固体
    G4VSolid* solid = new G4Tubs(soliname, 
                                inner_radius,  // 内半径
                                outer_radius,  // 外半径  
                                height / 2,    // 半高
                                0,             // 起始角度
                                2 * M_PI);     // 结束角度
    
    auto medium = G4Material::GetMaterial(material);
    mainLV = new G4LogicalVolume(solid, medium, logivol, 0, 0, 0);
    return true;
}