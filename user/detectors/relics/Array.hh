#pragma once

#include "BambooControl.hh"
#include "BambooDetector.hh"
#include "BambooFactory.hh"

class Array : public BambooDetector {

  public:
    Array(const std::string &n, const BambooParameters &pars);

    virtual bool construct(const BambooParameters &global_pars, BambooDetector *parent);

    static DetectorRegister<Array, std::string, BambooParameters> reg;

  private:
    // define additional parameters here
    G4LogicalVolume *CombineLV;
    G4LogicalVolume *SensitiveLV;
};

