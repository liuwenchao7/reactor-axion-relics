#include "RelicsOpticalRunAction.hh"
#include "RelicsOpticalDataManager.hh"
#include "GetMass.hh"

#include <G4Run.hh>
#include <G4UnitsTable.hh>

#include <CLHEP/Random/Randomize.h>

#include <chrono>

RelicsOpticalRunAction::~RelicsOpticalRunAction() {
    G4cout << "RelicsOpticalRunAction destructed." << G4endl;
}

void RelicsOpticalRunAction::BeginOfRunAction(const G4Run *aRun) {
    G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
    auto &dm = RelicsOpticalDataManager::getInstance();
    if (outputName == "")
        outputName = "bamboo_out.root";
    if (round == 0) {
        dm.book(outputName);
    } else {
        dm.book(outputName + "." + std::to_string(round));
    }
    dm.runId = aRun->GetRunID();
}

void RelicsOpticalRunAction::EndOfRunAction(const G4Run *aRun) {
    G4cout << "\n----------End of Run " << aRun->GetRunID() << "----------\n";
    RelicsOpticalDataManager::getInstance().save();
    RelicsOpticalDataManager::getInstance().close();
    ++round;

    if (getMass) {
        // G4UnitDefinition::PrintUnitsTable();
        G4cout << "\n----------Printing Mass of Geometry ----------\n";
        auto handle = new GetMass();
        handle->printMass(handle->getWorld());
        G4cout << "\n-----------Printed Mass of Geometry ----------\n";
    }
}
