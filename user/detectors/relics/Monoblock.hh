#pragma once

#include "BambooControl.hh"
#include "BambooDetector.hh"
#include "BambooFactory.hh"

class Monoblock : public BambooDetector {

  public:

    Monoblock (const std::string &n, const BambooParameters &pars) : BambooDetector(n, pars) {}

    virtual bool construct(const BambooParameters &global_pars, BambooDetector *parent);

    virtual bool constructMainLV(const BambooParameters &global_pars) = 0;

    auto getContainerPV() { return containerPV; }

  private:
    // define additional parameters here

  protected:
    G4VPhysicalVolume *containerPV = nullptr;
};

