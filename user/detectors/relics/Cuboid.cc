#include "Cuboid.hh"
#include "PandaXSensitiveDetector.hh"

#include <G4Material.hh>
#include <G4LogicalVolume.hh>
#include <G4Box.hh>

DetectorRegister<Cuboid, std::string, BambooParameters> Cuboid::reg("Cuboid");

Cuboid::Cuboid (const std::string &n, const BambooParameters &pars)
  : Monoblock(n, pars) {
    G4cout << "create detector Cuboid..." << G4endl;
}

bool Cuboid::constructMainLV(const BambooParameters &) {
    using namespace CLHEP;
    auto soliname = parameters.getParameter("soliname");
    auto logivol = parameters.getParameter("logivol");
    auto width_x = parameters.evaluateParameter("width_x");
    auto width_y = parameters.evaluateParameter("width_y");
    auto width_z = parameters.evaluateParameter("width_z");
    auto material = parameters.getParameter("material");
    auto optical_sister = parameters.getParameter("optical_sister");
    if (width_x == 0) {
        width_x = 4 * m;
    }
    if (width_y == 0) {
        width_y = 4 * m;
    }
    if (width_z == 0) {
        width_z = 4 * m;
    }
    auto solid =
        new G4Box(soliname, width_x / 2, width_y / 2, width_z / 2);
    auto medium = G4Material::GetMaterial(material);
    mainLV =
        new G4LogicalVolume(solid, medium, logivol, 0, 0, 0);
    return true;
}

