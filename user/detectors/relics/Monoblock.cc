#include "Monoblock.hh"
#include "Material.hh"
#include "PandaXSensitiveDetector.hh"

#include <G4ThreeVector.hh>
#include <G4PVPlacement.hh>
#include <G4SDManager.hh>

bool Monoblock::construct (const BambooParameters &pars,
                          BambooDetector *parent) {
    // add construction code here
    if (parent == nullptr || parent->getContainerLV() == nullptr) {
        return false;
    }
    using namespace CLHEP;
    auto shift_x = parameters.evaluateParameter("shift_x");
    auto shift_y = parameters.evaluateParameter("shift_y");
    auto shift_z = parameters.evaluateParameter("shift_z");
    auto physvol = parameters.getParameter("physvol");
    auto sendname = parameters.getParameter("sendname");
    auto optical_sister = parameters.getParameter("optical_sister");
    bool entrack = parameters.getParameter<bool>("entrack");
    bool surface = parameters.getParameter<bool>("surface");
    bool check_overlap = parameters.getParameter<bool>("check_overlap");
    bool is_construct = constructMainLV(pars);
    if (!is_construct) {
        return false;
    }
    mainPV =
        new G4PVPlacement(0, G4ThreeVector(shift_x, shift_y, shift_z), mainLV,
                          physvol, parent->getContainerLV(), false, 0, check_overlap);
    containerLV = mainLV;
    containerPV = mainPV;

    if (surface) {
        Material* pMaterial = new Material();
        pMaterial->AddLogicalSurface(mainPV, parent->getMainPV(), parameters);
    }

    if (sendname != "") {
        auto sd = new PandaXSensitiveDetector(sendname, true, entrack);
        G4SDManager::GetSDMpointer()->AddNewDetector(sd);
        mainLV->SetSensitiveDetector(sd);
    }
    return true;
}