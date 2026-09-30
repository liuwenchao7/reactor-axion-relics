#include "Cylinder.hh"
#include "PandaXSensitiveDetector.hh"

#include <G4Material.hh>
#include <G4LogicalVolume.hh>
#include <G4Tubs.hh>

DetectorRegister<Cylinder, std::string, BambooParameters> Cylinder::reg("Cylinder");

Cylinder::Cylinder (const std::string &n, const BambooParameters &pars)
  : Monoblock(n, pars) {
    G4cout << "create detector Cylinder..." << G4endl;
}

bool Cylinder::constructMainLV(const BambooParameters &) {
    using namespace CLHEP;
    auto soliname = parameters.getParameter("soliname");
    auto logivol = parameters.getParameter("logivol");
    auto height = parameters.evaluateParameter("height");
    auto radius = parameters.evaluateParameter("radius");
    auto material = parameters.getParameter("material");
    auto optical_sister = parameters.getParameter("optical_sister");
    if (height == 0) {
        height = 2 * m;
    }
    if (radius == 0) {
        radius = 1 * m;
    }
    auto solid =
        new G4Tubs(soliname, 0, radius, height / 2, 0, 2 * M_PI);
    auto medium = G4Material::GetMaterial(material);
    mainLV =
        new G4LogicalVolume(solid, medium, logivol, 0, 0, 0);
    return true;
}

