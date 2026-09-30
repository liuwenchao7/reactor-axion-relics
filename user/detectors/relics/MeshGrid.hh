#pragma once

#include "BambooControl.hh"
#include "BambooDetector.hh"
#include "BambooFactory.hh"
#include "Monoblock.hh"

#include <G4TwoVector.hh>

class MeshGrid : public Monoblock {

  public:
    MeshGrid(const std::string &n, const BambooParameters &pars);

    bool constructMainLV(const BambooParameters &global_pars);

    static DetectorRegister<MeshGrid, std::string, BambooParameters> reg;

  private:
    // define additional parameters here
    std::vector<G4TwoVector> GetPolygon(
      const G4int &sides,
      const G4double &radius
    );

  protected:
};

