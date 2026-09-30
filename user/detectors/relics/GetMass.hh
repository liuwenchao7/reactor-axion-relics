#pragma once

#include <set>

#include <G4String.hh>
#include <G4VPhysicalVolume.hh>

class GetMass {

  public:

    GetMass() = default;

    void printMass(G4VPhysicalVolume* mainPV);

    G4VPhysicalVolume* getWorld();

  private:

    static std::set<G4String> calculated;
};

