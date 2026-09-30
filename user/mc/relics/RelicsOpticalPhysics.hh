#pragma once

#include "BambooPhysics.hh"
#include "BambooFactory.hh"

class RelicsOpticalPhysics : public BambooPhysics {
  public:
    RelicsOpticalPhysics(const BambooParameters &pars);

    ~RelicsOpticalPhysics() = default;

    static PhysicsRegister<RelicsOpticalPhysics> reg;
};

