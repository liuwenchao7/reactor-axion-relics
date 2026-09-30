#pragma once

#include <globals.hh>
#include <G4UserStackingAction.hh>

class G4Track;

class RelicsOpticalStackingAction : public G4UserStackingAction {
  public:
    RelicsOpticalStackingAction() = default;

    ~RelicsOpticalStackingAction() {}

    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track *aTrack);

    virtual void NewStage();

    virtual void PrepareNewEvent();

  private:
    G4int fScintillationCounter = 0;
    G4int fCerenkovCounter = 0;
};
