#include "RelicsOpticalAnalysis.hh"

#include "RelicsOpticalDataManager.hh"
#include "RelicsOpticalEventAction.hh"
#include "RelicsOpticalRunAction.hh"
#include "RelicsOpticalSteppingAction.hh"
#include "RelicsOpticalTrackingAction.hh"
#include "RelicsOpticalStackingAction.hh"

AnalysisRegister<RelicsOpticalAnalysis>
    RelicsOpticalAnalysis::reg("RelicsOpticalAnalysis");

RelicsOpticalAnalysis::RelicsOpticalAnalysis(const BambooParameters &pars)
    : BambooAnalysis{pars} {
    bool enable_energy_deposition = true;
    if (analysisParameters.getParameters().find("EnableEnergyDeposition") !=
        analysisParameters.getParameters().end())
        enable_energy_deposition =
            analysisParameters.getParameter<bool>("EnableEnergyDeposition");
    bool enable_primary_particle =
        analysisParameters.getParameter<bool>("EnablePrimaryParticle");
    bool save_null_events =
        analysisParameters.getParameter<bool>("SaveNullEvents");
    bool get_geo_mass =
        analysisParameters.getParameter<bool>("GetGeoMass");

    auto &dm = RelicsOpticalDataManager::getInstance();
    dm.setRecordEnergyDeposition(enable_energy_deposition);
    dm.setRecordPrimaryParticle(enable_primary_particle);
    dm.setRecordNullEvent(save_null_events);
    auto sdName = analysisParameters.getParameter("OpticalSdName");
    long user_seed = analysisParameters.getParameter<long>("UserSeed");

    runAction = new RelicsOpticalRunAction(
        analysisParameters.getParameter("DataFileName"), get_geo_mass);
    eventAction = new RelicsOpticalEventAction();
    if (user_seed > 0) {
        dynamic_cast<RelicsOpticalEventAction *>(eventAction)->setUserSeed(user_seed);
    }

    stackingAction = new RelicsOpticalStackingAction();

    steppingAction = new RelicsOpticalSteppingAction(sdName);

    trackingAction = new RelicsOpticalTrackingAction();
}
