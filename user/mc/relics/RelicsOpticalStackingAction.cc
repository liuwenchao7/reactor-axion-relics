#include "RelicsOpticalStackingAction.hh"
#include "RelicsOpticalDataManager.hh"
#include "RelicsOpticalEventInformation.hh"

#include <G4Track.hh>
#include <G4EventManager.hh>
#include <G4OpticalPhoton.hh>

G4ClassificationOfNewTrack
RelicsOpticalStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
  { // particle is optical photon
    if(aTrack->GetParentID() > 0)
    { // particle is secondary
      if(aTrack->GetCreatorProcess()->GetProcessName() == "Scintillation")
        fScintillationCounter++;
      if(aTrack->GetCreatorProcess()->GetProcessName() == "Cerenkov")
        fCerenkovCounter++;
    }
  }
  return fUrgent;
}

void RelicsOpticalStackingAction::NewStage() {
    auto aEvent = G4EventManager::GetEventManager()->GetNonconstCurrentEvent();
    auto info = new RelicsOpticalEventInformation();
    info->SetScintillation(fScintillationCounter);
    info->SetCerenkov(fCerenkovCounter);
    aEvent->SetUserInformation(info);
}

void RelicsOpticalStackingAction::PrepareNewEvent() {
    fScintillationCounter = 0;
    fCerenkovCounter = 0;
}
