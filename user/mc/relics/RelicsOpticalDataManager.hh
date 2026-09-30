#pragma once

#include <string>
#include <vector>

#include <TFile.h>
#include <G4Event.hh>

class TTree;

// class G4Event;

class RelicsOpticalDataManager {

    friend class RelicsOpticalRunAction;
    friend class RelicsOpticalEventAction;
    // friend class RelicsAnalysis;

  private:
    RelicsOpticalDataManager() = default;
    RelicsOpticalDataManager(const RelicsOpticalDataManager &) = delete;
    RelicsOpticalDataManager &
    operator=(const RelicsOpticalDataManager &) = delete;

  public:
    static RelicsOpticalDataManager &getInstance() {
        static RelicsOpticalDataManager instance;
        return instance;
    }

    ~RelicsOpticalDataManager();

    void book(const std::string &name = "relicsout.root");

    void save();

    void close();

    void setRecordEnergyDeposition(bool t) { recordEnergyDeposition = t; }
    void setRecordPrimaryParticle(bool t) { recordPrimaryParticle = t; }
    void setRecordNullEvent(bool t) { recordNullEvent = t; }

    void fillEvent(const G4Event *aEvent);

  private:
    void resetData();

    TFile *rootFile = nullptr;
    TTree *mcTree = nullptr;

    int runId;
    int eventId;
    int64_t seed;

    int nHits;
    int nScintillation;
    int nCerenkov;
    std::vector<int> pmtNumber;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> energy;
    std::vector<double> velocity;
    std::vector<std::string> creatorProcess;
    std::vector<double> tGlobal;
    std::vector<double> tLocal;
    std::vector<double> tProper;
    std::vector<double> ox;
    std::vector<double> oy;
    std::vector<double> oz;
    std::vector<std::string> parent;
    std::vector<int> photonParentId;

    // variables for primary particle
    int nPrimaries;
    std::vector<std::string> primaryType;
    std::vector<int> primaryId;
    std::vector<double> primaryEnergy;
    std::vector<double> primaryPx;
    std::vector<double> primaryPy;
    std::vector<double> primaryPz;
    std::vector<double> primaryX;
    std::vector<double> primaryY;
    std::vector<double> primaryZ;

    // variables related to energy deposition
    int nEnergyHits;
    double totalEnergy;
    std::vector<int> trackId;
    std::vector<int> parentId;
    std::vector<std::string> type;
    std::vector<std::string> parentType;
    std::vector<std::string> volume;
    std::vector<std::string> depositionProcess;
    std::vector<double> xd;
    std::vector<double> yd;
    std::vector<double> zd;
    std::vector<double> td;
    std::vector<double> hitEnergy;

    bool recordEnergyDeposition = true;
    bool recordPrimaryParticle = true;
    bool recordNullEvent = true;
};
