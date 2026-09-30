#pragma once

#include <G4UserEventAction.hh>


class G4Event;

class RelicsOpticalEventAction : public G4UserEventAction {
  public:
    RelicsOpticalEventAction() = default;

    virtual ~RelicsOpticalEventAction();

    virtual void BeginOfEventAction(const G4Event *aEvent);

    virtual void EndOfEventAction(const G4Event *aEvent);

    void setUserSeed(long seed) { user_seed = seed; }

  private:
    long user_seed = -1;
};
