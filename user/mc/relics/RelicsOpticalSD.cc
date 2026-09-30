#include "RelicsOpticalSD.hh"
#include "RelicsOpticalTrackInformation.hh"
#include "PandaXTrackInformation.hh"

#include <G4OpticalPhoton.hh>
#include <G4VProcess.hh>

RelicsOpticalSD::RelicsOpticalSD(const G4String &name, G4bool recordDeposition)
    : G4VSensitiveDetector(name), recordEnergyDeposition(recordDeposition) {
    collectionName.push_back("OpticalHits");
    collectionName.push_back("EnergyDepositionHits");
}

void RelicsOpticalSD::Initialize(G4HCofThisEvent *hc) {
    opHitsCollection = new RelicsOpticalHitsCollection(SensitiveDetectorName,
                                                       collectionName[0]);
    hc->AddHitsCollection(GetCollectionID(0), opHitsCollection);
    eDHitsCollection = new PandaXEnergyDepositionHitsCollection(
        SensitiveDetectorName, collectionName[1]);
    hc->AddHitsCollection(GetCollectionID(1), eDHitsCollection);
}

void RelicsOpticalSD::EndOfEvent(G4HCofThisEvent *) {}

G4bool RelicsOpticalSD::ProcessHits(G4Step *aStep, G4TouchableHistory *) {
    if (!recordEnergyDeposition) {
        return false;
    }
    auto track = aStep->GetTrack();
    if (track->GetParticleDefinition() ==
        G4OpticalPhoton::OpticalPhotonDefinition())
        return false;
    G4double edep = aStep->GetTotalEnergyDeposit();
    if (edep == 0)
        return false;
    PandaXEnergyDepositionHit *hit = new PandaXEnergyDepositionHit();
    int hitId = track->GetTrackID();
    hit->setTrackId(hitId);
    int parentId = track->GetParentID();
    hit->setParentId(parentId);
    if (parentId != 0) {
        auto trackInfo =
            static_cast<PandaXTrackInformation *>(track->GetUserInformation());
        hit->setParent(trackInfo->GetParent());
    }
    hit->setEnergy(edep);
    hit->setX(aStep->GetPostStepPoint()->GetPosition().x());
    hit->setY(aStep->GetPostStepPoint()->GetPosition().y());
    hit->setZ(aStep->GetPostStepPoint()->GetPosition().z());
    hit->setT(aStep->GetPostStepPoint()->GetGlobalTime());
    hit->setType(aStep->GetTrack()->GetParticleDefinition()->GetParticleName());
    if (aStep->GetTrack()->GetCreatorProcess()) {
        hit->setCreatorProcess(
            aStep->GetTrack()->GetCreatorProcess()->GetProcessName());
    }
    hit->setDepositionProcess(
        aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());
    eDHitsCollection->insert(hit);
    return true;
}

G4bool RelicsOpticalSD::ProcessHits_constStep(const G4Step *aStep,
                                              G4TouchableHistory *) {
    auto track = aStep->GetTrack();
    if (track->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
        return false;
    auto pmt_number =
        aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(2);
        // we need to know that here `2` is just hard code
        // we can use `G4VTouchable::GetVolume` to search the PV name is geometry tree

    auto info = static_cast<RelicsOpticalTrackInformation *>(
        track->GetUserInformation());
    auto hit = new RelicsOpticalHit();
    hit->setPMTNumber(pmt_number);
    hit->setEnergy(track->GetKineticEnergy());
    hit->setPosition(aStep->GetPostStepPoint()->GetPosition());
    if (track->GetCreatorProcess())
        hit->setCreatorProcess(track->GetCreatorProcess()->GetProcessName());
    else
        hit->setCreatorProcess("unknown");
    hit->setVelocity(track->GetVelocity());
    hit->setGlobalTime(track->GetGlobalTime());
    hit->setLocalTime(track->GetLocalTime());
    hit->setProperTime(track->GetProperTime());
    if (info) {
        hit->setSourcePos(info->GetSourcePos());
        hit->setParent(info->GetParent());
        hit->setParentId(track->GetParentID());
    }
    opHitsCollection->insert(hit);
    return true;
}
