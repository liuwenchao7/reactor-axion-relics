#pragma once

#include "BambooControl.hh"
#include "BambooDetector.hh"
#include "BambooFactory.hh"
#include "Monoblock.hh"

class Polyhedra : public Monoblock {

  public:
    Polyhedra(const std::string &n, const BambooParameters &pars);

    bool constructMainLV(const BambooParameters &global_pars);

    static DetectorRegister<Polyhedra, std::string, BambooParameters> reg;

  private:
    // define additional parameters here

  protected:
};

