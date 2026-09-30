#include "RelicsOpticalDataManager.hh"
#include "PandaXEnergyDepositionHit.hh"
#include "RelicsOpticalHit.hh"
#include "RelicsOpticalEventInformation.hh"

#include <TTree.h>

#include <G4Event.hh>
#include <G4ParticleDefinition.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4THitsMap.hh>

RelicsOpticalDataManager::~RelicsOpticalDataManager() {
    G4cout << "RelicsOpticalDataManager destructed." << G4endl;
}

void RelicsOpticalDataManager::book(const std::string &name) {
    rootFile = new TFile(name.c_str(), "RECREATE");
    if (rootFile->IsZombie()) {
        std::cerr << "RelicsOpticalDataManager::book: ";
        std::cerr << "problem creating ROOT file." << std::endl << std::endl;
        rootFile = nullptr;
        return;
    }
    mcTree = new TTree("mcTree", "Tree with Optical MC info");
    mcTree->Branch("runId", &runId, "runId/I");
    mcTree->Branch("eventId", &eventId, "eventId/I");
    mcTree->Branch("seed", &seed, "seed/L");
    mcTree->Branch("nHits", &nHits, "nHits/I");
    mcTree->Branch("nScintillation", &nScintillation, "nScintillation/I");
    mcTree->Branch("nCerenkov", &nCerenkov, "nCerenkov/I");
    mcTree->Branch("pmtNumber", &pmtNumber);
    mcTree->Branch("creatorProcess", &creatorProcess);
    mcTree->Branch("x", &x);
    mcTree->Branch("y", &y);
    mcTree->Branch("z", &z);
    mcTree->Branch("tGlobal", &tGlobal);
    mcTree->Branch("tLocal", &tLocal);
    mcTree->Branch("tProper", &tProper);
    mcTree->Branch("energy", &energy);
    mcTree->Branch("velocity", &velocity);
    mcTree->Branch("ox", &ox);
    mcTree->Branch("oy", &oy);
    mcTree->Branch("oz", &oz);
    mcTree->Branch("parent", &parent);
    mcTree->Branch("photonParentId", &photonParentId);

    if (recordEnergyDeposition) {
        mcTree->Branch("nEnergyHits", &nEnergyHits, "nEnergyHits/I");
        mcTree->Branch("trackId", &trackId);
        mcTree->Branch("parentId", &parentId);
        mcTree->Branch("type", &type);
        mcTree->Branch("parentType", &parentType);
        mcTree->Branch("depositionProcess", &depositionProcess);
        mcTree->Branch("volume", &volume);
        mcTree->Branch("totalEnergy", &totalEnergy);
        mcTree->Branch("xd", &xd);
        mcTree->Branch("yd", &yd);
        mcTree->Branch("zd", &zd);
        mcTree->Branch("td", &td);
        mcTree->Branch("hitEnergy", &hitEnergy);
    }
    if (recordPrimaryParticle) {
        mcTree->Branch("nPrimaries", &nPrimaries);
        mcTree->Branch("primaryType", &primaryType);
        mcTree->Branch("primaryId", &primaryId);
        mcTree->Branch("primaryEnergy", &primaryEnergy);
        mcTree->Branch("primaryPx", &primaryPx);
        mcTree->Branch("primaryPy", &primaryPy);
        mcTree->Branch("primaryPz", &primaryPz);
        mcTree->Branch("primaryX", &primaryX);
        mcTree->Branch("primaryY", &primaryY);
        mcTree->Branch("primaryZ", &primaryZ);
    }
}

void RelicsOpticalDataManager::save() {
    mcTree->Write("", TObject::kWriteDelete);
    G4cout << "RelicsOpticalDataManager: Write data to tree..." << G4endl;
}

void RelicsOpticalDataManager::close() {
    if (rootFile != nullptr) {
        rootFile->Close();
        delete rootFile;
        rootFile = nullptr;
    }
}

void RelicsOpticalDataManager::fillEvent(const G4Event *aEvent) {
    auto hCthis = aEvent->GetHCofThisEvent();
    int nHitCollections = hCthis->GetNumberOfCollections();
    resetData();
    if (recordEnergyDeposition) {
        for (int i = 0; i < nHitCollections; ++i) {
            G4VHitsCollection *hitsCollection = hCthis->GetHC(i);
            if (hitsCollection->GetName().contains("EnergyDepositionHits")) {
                PandaXEnergyDepositionHitsCollection *hC =
                    (PandaXEnergyDepositionHitsCollection *)hitsCollection;
                for (size_t j = 0; j < hitsCollection->GetSize(); ++j) {
                    PandaXEnergyDepositionHit *hit =
                        (PandaXEnergyDepositionHit *)hC->GetHit(j);
                    trackId.push_back(hit->getTrackId());
                    parentId.push_back(hit->getParentId());
                    type.push_back(hit->getType());
                    parentType.push_back(hit->getParent());
                    depositionProcess.push_back(hit->getDepositionProcess());
                    volume.push_back(hitsCollection->GetSDname());
                    xd.push_back(hit->getX() / mm);
                    yd.push_back(hit->getY() / mm);
                    zd.push_back(hit->getZ() / mm);
                    td.push_back(hit->getT() / s);
                    hitEnergy.push_back(hit->getEnergy() / keV);
                    totalEnergy += (hit->getEnergy()) / keV;
                    nEnergyHits++;
                }
            }
        }
    }
    for (int i = 0; i < nHitCollections; ++i) {
        auto hC = hCthis->GetHC(i);
        if (hC->GetName().contains("OpticalHits")) {
            auto hCo = static_cast<RelicsOpticalHitsCollection *>(hC);
            for (auto j = 0u; j < hC->GetSize(); ++j) {
                auto hit = static_cast<RelicsOpticalHit *>(hCo->GetHit(j));
                pmtNumber.push_back(hit->getPMTNumber());
                x.push_back(hit->getPosition().x() / mm);
                y.push_back(hit->getPosition().y() / mm);
                z.push_back(hit->getPosition().z() / mm);
                energy.push_back(hit->getEnergy() / keV);
                velocity.push_back(hit->getVelocity() / c_light);
                creatorProcess.push_back(hit->getCreatorProcess());
                tGlobal.push_back(hit->getGlobalTime() / ns);
                tLocal.push_back(hit->getLocalTime() / ns);
                tProper.push_back(hit->getProperTime() / ns);
                ox.push_back(hit->getSourcePos().x() / mm);
                oy.push_back(hit->getSourcePos().y() / mm);
                oz.push_back(hit->getSourcePos().z() / mm);
                parent.push_back(hit->getParent());
                photonParentId.push_back(hit->getParentId());
                nHits++;
            }
        }
    }
    auto eventInfo =
        static_cast<RelicsOpticalEventInformation *>(aEvent->GetUserInformation());
    nCerenkov = eventInfo->GetCerenkov();
    nScintillation = eventInfo->GetScintillation();
    if (recordPrimaryParticle) {
        // loop over all primary particles...
        int nVertex = aEvent->GetNumberOfPrimaryVertex();
        for (int iV = 0; iV < nVertex; ++iV) {
            G4PrimaryVertex *vertex = aEvent->GetPrimaryVertex(iV);
            int nParticles = vertex->GetNumberOfParticle();
            double vx = vertex->GetX0();
            double vy = vertex->GetY0();
            double vz = vertex->GetZ0();
            for (int ip = 0; ip < nParticles; ++ip) {
                G4PrimaryParticle *particle = vertex->GetPrimary(ip);
                auto pdef = particle->GetParticleDefinition();
                if (pdef != nullptr) {
                    primaryType.push_back(pdef->GetParticleName());
                } else {
                    primaryType.push_back("unknown");
                }
                primaryId.push_back(particle->GetTrackID());
                primaryEnergy.push_back((particle->GetTotalEnergy() - particle->GetMass()) / keV);
                primaryPx.push_back(particle->GetPx() / keV);
                primaryPy.push_back(particle->GetPy() / keV);
                primaryPz.push_back(particle->GetPz() / keV);
                primaryX.push_back(vx / mm);
                primaryY.push_back(vy / mm);
                primaryZ.push_back(vz / mm);
                nPrimaries++;
            }
        }
    }
    if ((nHits > 0) ||
        (recordNullEvent && recordPrimaryParticle))
        mcTree->Fill();
}

void RelicsOpticalDataManager::resetData() {
    nHits = 0;
    nScintillation = 0;
    nCerenkov = 0;
    x.clear();
    y.clear();
    z.clear();
    energy.clear();
    creatorProcess.clear();
    tGlobal.clear();
    tLocal.clear();
    tProper.clear();
    pmtNumber.clear();
    ox.clear();
    oy.clear();
    oz.clear();
    parent.clear();
    photonParentId.clear();
    velocity.clear();

    if (recordEnergyDeposition) {
        nEnergyHits = 0;
        trackId.clear();
        parentId.clear();
        type.clear();
        parentType.clear();
        volume.clear();
        depositionProcess.clear();
        totalEnergy = 0;
        xd.clear();
        yd.clear();
        zd.clear();
        td.clear();
        hitEnergy.clear();
    }
    if (recordPrimaryParticle) {
        nPrimaries = 0;
        primaryType.clear();
        primaryId.clear();
        primaryEnergy.clear();
        primaryPx.clear();
        primaryPy.clear();
        primaryPz.clear();
        primaryX.clear();
        primaryY.clear();
        primaryZ.clear();
    }
}
