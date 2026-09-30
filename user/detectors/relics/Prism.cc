#include "Prism.hh"
#include "PandaXSensitiveDetector.hh"

#include <G4Material.hh>
#include <G4LogicalVolume.hh>
#include <G4Polyhedra.hh>
#include <G4GenericPolycone.hh>

DetectorRegister<Prism, std::string, BambooParameters> Prism::reg("Prism");

Prism::Prism (const std::string &n, const BambooParameters &pars)
  : Monoblock(n, pars) {
    G4cout << "create detector Prism..." << G4endl;
}

bool Prism::constructMainLV(const BambooParameters &) {
    using namespace CLHEP;
    auto soliname = parameters.getParameter("soliname");
    auto logivol = parameters.getParameter("logivol");
    auto height = parameters.evaluateParameter("height");
    // radius of inscribed circle
    auto outer_radius = parameters.evaluateParameter("outer_radius");
    auto inner_radius = parameters.evaluateParameter("inner_radius");
    int sides = parameters.getParameter<int>("sides");
    auto material = parameters.getParameter("material");
    if (height == 0) {
        height = 2 * m;
    }
    if (outer_radius == 0) {
        outer_radius = 1 * m;
    }
    // G4double outer_r = outer_radius / std::cos(M_PI / sides);
    // G4double inner_r = inner_radius / std::cos(M_PI / sides);
    // G4double r[4] = {inner_r, outer_r, outer_r, inner_r};
    G4double z[4] = {-height / 2, -height / 2, height / 2, height / 2};
    G4VSolid* solid;
    if (sides == 0) {
        G4double r[4] = {inner_radius, outer_radius, outer_radius, inner_radius};
        solid = new G4GenericPolycone(soliname, 0, 2 * M_PI, 4, r, z);
    } else {
        G4double outer_r = outer_radius / std::cos(M_PI / sides);
        G4double inner_r = inner_radius / std::cos(M_PI / sides);
        G4double r[4] = {inner_r, outer_r, outer_r, inner_r};
        solid = new G4Polyhedra(soliname, 0, 2 * M_PI, sides, 4, r, z);
    }
    auto medium = G4Material::GetMaterial(material);
    mainLV =
        new G4LogicalVolume(solid, medium, logivol, 0, 0, 0);
    return true;
}

