#include "GetMass.hh"

#include <sstream>

#include <G4String.hh>
#include <G4Material.hh>
#include <G4UnitsTable.hh>
#include <G4VPhysicalVolume.hh>
#include <G4TransportationManager.hh>

std::set<G4String> GetMass::calculated = {};

G4VPhysicalVolume* GetMass::getWorld() {
    auto world = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume();
    return world;
}

void GetMass::printMass (G4VPhysicalVolume* mainPV) {
    auto mainLV = mainPV->GetLogicalVolume();
    auto name = mainPV->GetName();
    if (calculated.find(name) != calculated.end()) {
        return;
    }
    auto material = mainLV->GetMaterial()->GetName();
    // https://github.com/Geant4/geant4/blob/geant4-11.0-release/source/geometry/management/src/G4VSolid.cc#L190
    auto mass = mainLV->GetMass(false, false, nullptr);
    auto density = mainLV->GetMaterial()->GetDensity();
    auto volume = mainLV->GetSolid()->GetCubicVolume();
    auto area = mainLV->GetSolid()->GetSurfaceArea();
    // G4cout << name << " : " << material << " : " << G4BestUnit(mass, "Mass") << G4endl;
    std::stringstream ss;
    ss << std::setw(25) << name << ":"\
    "" << std::setw(25) << material << ":"\
    "" << std::setw(25) << G4BestUnit(mass,    "Mass") << ":"\
    "" << std::setw(25) << G4BestUnit(density, "Volumic Mass") << ":"\
    "" << std::setw(25) << G4BestUnit(volume,  "Volume") << ":"\
    "" << std::setw(25) << G4BestUnit(area,    "Surface");
	G4String str = ss.str();
    // ss >> str;
    G4cout << str << G4endl;
    auto nDaughters = mainLV->GetNoDaughters();
    G4VPhysicalVolume* sisterPV;
    for (int id = 0; id < nDaughters; ++id) {
        sisterPV = mainLV->GetDaughter(id);
        printMass(sisterPV);
    }
    calculated.insert(name);
}

