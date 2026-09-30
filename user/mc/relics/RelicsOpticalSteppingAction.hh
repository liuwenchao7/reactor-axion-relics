#pragma once

#include <G4UserSteppingAction.hh>

#include <string>

class G4Step;

class RelicsOpticalSteppingAction : public G4UserSteppingAction {
  public:
    RelicsOpticalSteppingAction(const std::string &s) : sdName{s} {}

    ~RelicsOpticalSteppingAction() = default;

    virtual void UserSteppingAction(const G4Step *aStep);

  private:
    std::string sdName;
};
