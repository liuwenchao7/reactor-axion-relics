#include "Polyhedra.hh"
#include "PandaXSensitiveDetector.hh"

#include <G4Material.hh>
#include <G4LogicalVolume.hh>
#include <G4Polyhedra.hh>
#include <G4GenericPolycone.hh>

DetectorRegister<Polyhedra, std::string, BambooParameters> Polyhedra::reg("Polyhedra");

Polyhedra::Polyhedra (const std::string &n, const BambooParameters &pars)
  : Monoblock(n, pars) {
    G4cout << "create detector Polyhedra..." << G4endl;
}

bool Polyhedra::constructMainLV(const BambooParameters &) {
    using namespace CLHEP;
    auto soliname = parameters.getParameter("soliname");
    auto logivol = parameters.getParameter("logivol");
    std::string source = parameters.getParameter("source");
    // default value of sides is 0
    int sides = parameters.getParameter<int>("sides");
    auto material = parameters.getParameter("material");
    auto optical_sister = parameters.getParameter("optical_sister");
    std::ifstream source_stream(source);
    G4int index, numRZ;
    double R, Z;
    if (source_stream >> numRZ) {
        /* Unimplemented at the moment */
    } else {
        throw std::runtime_error("Array: error read data");
    }
    // G4double r[numRZ];
    // G4double z[numRZ];
    G4double *r = new G4double[numRZ];
    G4double *z = new G4double[numRZ];
    for (int il = 0; il < numRZ; ++il) {
        source_stream >> index >> R >> Z;
        if (source_stream.fail())
            throw std::runtime_error("Array: error read data");
        r[il] = R * cm;
        z[il] = Z * cm;
    }
    G4VSolid* solid;
    if (sides == 0) {
        solid = new G4GenericPolycone(soliname, 0, 2 * M_PI, numRZ, r, z);
    } else {
        solid = new G4Polyhedra(soliname, 0, 2 * M_PI, sides, numRZ, r, z);
    }
    auto medium = G4Material::GetMaterial(material);
    mainLV =
        new G4LogicalVolume(solid, medium, logivol, 0, 0, 0);
    delete[] r;
    delete[] z;
    return true;
}

