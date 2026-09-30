#pragma once

#include <G4UserRunAction.hh>

class G4Run;

class RelicsOpticalRunAction : public G4UserRunAction {
  public:
    RelicsOpticalRunAction(
      const std::string &name = "mc_out.root",
      bool get_geo_mass = false
    ) : outputName{name}, getMass{get_geo_mass} {};

    virtual ~RelicsOpticalRunAction();

    virtual void BeginOfRunAction(const G4Run *aRun);

    virtual void EndOfRunAction(const G4Run *aRun);

  private:
    std::string outputName;
    bool getMass = false;
    int round = 0;
};
