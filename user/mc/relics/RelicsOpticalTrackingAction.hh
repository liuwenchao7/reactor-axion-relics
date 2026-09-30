#pragma once

#include <G4UserTrackingAction.hh>

class G4Track;

class RelicsOpticalTrackingAction : public G4UserTrackingAction {
  public:
    RelicsOpticalTrackingAction() = default;

    ~RelicsOpticalTrackingAction() = default;

    virtual void PreUserTrackingAction(const G4Track *aTrack);
    virtual void PostUserTrackingAction(const G4Track *aTrack);

  private:
};
