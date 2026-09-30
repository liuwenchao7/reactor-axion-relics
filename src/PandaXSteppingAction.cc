#include "PandaXSteppingAction.hh"
#include "PandaXDataManager.hh"

#include <G4Ions.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4VProcess.hh>
#include <G4Version.hh>
#include <G4ThreeVector.hh>

PandaXSteppingAction::~PandaXSteppingAction() {}

void PandaXSteppingAction::UserSteppingAction(const G4Step *aStep) {
    G4Track *aTrack = aStep->GetTrack();
    auto *def = aTrack->GetParticleDefinition();

    // ===== NEW: record gamma Compton scatterings inside the LXe TPC =====
    // The gamma's Compton step has ZERO energy deposition, so it never reaches
    // the sensitive detector (ProcessHits). Its PostStepPoint momentum is the
    // outgoing (scattered) photon direction -> capture it here.
    if (def->GetParticleName() == "gamma") {
        const G4StepPoint *post = aStep->GetPostStepPoint();
        const G4StepPoint *pre = aStep->GetPreStepPoint();
        const G4VProcess *proc = post->GetProcessDefinedStep();
        if (proc != nullptr && proc->GetProcessName() == "compt") {
            const G4VPhysicalVolume *pv = post->GetPhysicalVolume();
            G4String volName = pv ? pv->GetName() : "unknown";
            // record ALL gamma Comptons (any volume) for diagnosis
            if (true) {
                G4ThreeVector pos = post->GetPosition();
                G4ThreeVector mom = post->GetMomentum();
                G4ThreeVector pmom = pre->GetMomentum();
                PandaXDataManager::getInstance().addGammaCompRecord(
                    pos.x(), pos.y(), pos.z(),
                    mom.x(), mom.y(), mom.z(),
                    post->GetGlobalTime(),
                    aTrack->GetTrackID(),
                    pmom.x(), pmom.y(), pmom.z(),
                    volName);
            }
        }
        return;
    }
    // ===== end NEW =====

    if (def->GetParticleType() != "nucleus" || aTrack->GetTrackID() == 1) {
        return;
    }
    if (aTrack->GetTrackStatus() != fStopButAlive)
        return;
    if (def->GetIonLifeTime() > chainSplittingLifeTime) {
        aTrack->SetTrackStatus(fPostponeToNextEvent);
        aTrack->SetGlobalTime(0.0 * second);
    }
}

