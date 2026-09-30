#pragma once

#include <G4ThreeVector.hh>
#include <G4VUserEventInformation.hh>

class RelicsOpticalEventInformation : public G4VUserEventInformation {
  public:
    RelicsOpticalEventInformation() {};

    inline virtual void Print() const {};

    const G4int &GetScintillation() const { return fScintillationCounter; }

    void SetScintillation(const G4int &i) { fScintillationCounter = i; }

    const G4int &GetCerenkov() const { return fCerenkovCounter; }

    void SetCerenkov(const G4int &i) { fCerenkovCounter = i; }

  private:
    G4int fScintillationCounter;
    G4int fCerenkovCounter;
};
