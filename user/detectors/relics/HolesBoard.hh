#pragma once

#include "BambooControl.hh"
#include "BambooDetector.hh"
#include "BambooFactory.hh"
#include "Monoblock.hh"

class HolesBoard : public Monoblock {

  public:
    HolesBoard(const std::string &n, const BambooParameters &pars);

    bool constructMainLV(const BambooParameters &global_pars);

    static DetectorRegister<HolesBoard, std::string, BambooParameters> reg;

  private:
    // define additional parameters here

  protected:
};

