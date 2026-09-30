#pragma once

#include "BambooControl.hh"
#include "BambooDetector.hh"
#include "BambooFactory.hh"

class World : public BambooDetector {

  public:
    World(const std::string &n, const BambooParameters &pars);

    virtual bool construct(const BambooParameters &global_pars, BambooDetector *parent);

    static DetectorRegister<World, std::string, BambooParameters> reg;

    auto getContainerPV() { return containerPV; }

  private:
    // define additional parameters here

  protected:
    G4VPhysicalVolume *containerPV = nullptr;
};

