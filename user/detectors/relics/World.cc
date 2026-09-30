#include "World.hh"

#include <G4Material.hh>
#include <G4ThreeVector.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4VisAttributes.hh>
#include <G4Box.hh>

DetectorRegister<World, std::string, BambooParameters> World::reg("World");

World::World (const std::string &n, const BambooParameters &pars)
  : BambooDetector(n, pars) {
    G4cout << "create detector World..." << G4endl;
}

bool World::construct (const BambooParameters &, BambooDetector *) {
    // add construction code here
    using namespace CLHEP;
    G4Material *air = G4Material::GetMaterial("G4_AIR");
    auto half_x = parameters.evaluateParameter("half_x");
    auto half_y = parameters.evaluateParameter("half_y");
    auto half_z = parameters.evaluateParameter("half_z");
    if (half_x == 0) {
        half_x = 10 * m;
    }
    if (half_y == 0) {
        half_y = 10 * m;
    }
    if (half_z == 0) {
        half_z = 10 * m;
    }
    auto worldBox = new G4Box("WorldBox", half_x, half_y, half_z);
    mainLV = new G4LogicalVolume(worldBox, air, "WorldLog", 0, 0, 0);
    mainPV =
        new G4PVPlacement(0, G4ThreeVector(), mainLV, "World", 0, false, 0);
    containerLV = mainLV;
    containerPV = mainPV;
    auto vis = new G4VisAttributes();
    vis->SetVisibility(false);
    mainLV->SetVisAttributes(vis);
    return true;
}

