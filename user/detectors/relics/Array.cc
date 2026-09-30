#include "Array.hh"
#include "PMTsR8520.hh"
#include "Cuboid.hh"
#include "Prism.hh"
#include "PandaXSensitiveDetector.hh"
#include "PandaXOpticalSD.hh"

#include <G4Material.hh>
#include <G4ThreeVector.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4VisAttributes.hh>
#include <G4Transform3D.hh>
#include <G4SDManager.hh>
#include <G4VSensitiveDetector.hh>

DetectorRegister<Array, std::string, BambooParameters> Array::reg("Array");

Array::Array (const std::string &n, const BambooParameters &pars)
  : BambooDetector(n, pars) {
    G4cout << "create detector Array..." << G4endl;
}

bool Array::construct (const BambooParameters &global_pars, BambooDetector *parent) {
    // add construction code here
    if (parent == nullptr || parent->getContainerLV() == nullptr) {
        return false;
    }

    using namespace CLHEP;
    std::string objects = parameters.getParameter("objects");
    if (objects == "PMTsR8520") {
      PMTsR8520 *pPMTsR8520 = new PMTsR8520(name, parameters);
      PMTLVs pmtlvs = pPMTsR8520->construct(global_pars, parent);
      CombineLV = pmtlvs.CombineLV;
      SensitiveLV = pmtlvs.PhotoCathodeLV;
    } else if (objects == "Cuboid") {
      Cuboid *pCuboid = new Cuboid(name, parameters);
      pCuboid->constructMainLV(global_pars);
      CombineLV = pCuboid->getMainLV();
      SensitiveLV = pCuboid->getMainLV();
    } else if (objects == "Prism") {
      Prism *pPrism = new Prism(name, parameters);
      pPrism->constructMainLV(global_pars);
      CombineLV = pPrism->getMainLV();
      SensitiveLV = pPrism->getMainLV();
    }
    else {
      return false;
    }

    auto physvol = parameters.getParameter("physvol");
    std::string source = parameters.getParameter("source");
    auto sendname = parameters.getParameter("sendname");
    bool opticalsd = parameters.getParameter<bool>("opticalsd");
    bool check_overlap = parameters.getParameter<bool>("check_overlap");
    std::ifstream source_stream(source);
    int index, nObjects;
    double X, Y, Z, RotY, RotZ;
    if (source_stream >> nObjects) {
      for (int il = 0; il < nObjects; ++il) {
        source_stream >> index >> X >> Y >> Z >> RotY >> RotZ;
        if (source_stream.fail())
          throw std::runtime_error("Array: error read data");
        G4RotationMatrix *mRot = new G4RotationMatrix();
        mRot->rotateZ(-RotZ * deg);
        mRot->rotateY(-RotY * deg);
        new G4PVPlacement(mRot,
                          G4ThreeVector(X * cm, Y * cm, Z * cm),
                          CombineLV, physvol,
                          parent->getContainerLV(), false, index, check_overlap);
      }
    } else {
      throw std::runtime_error("Array: error read data");
    }
    containerLV = CombineLV;
    if (sendname != "") {
      G4VSensitiveDetector* sd = NULL;
      if (opticalsd) {
        sd = new PandaXOpticalSD(sendname);
      } else {
        bool entrack = parameters.getParameter<bool>("entrack");
        sd = new PandaXSensitiveDetector(sendname, true, entrack);
      }
      G4SDManager::GetSDMpointer()->AddNewDetector(sd);
      SensitiveLV->SetSensitiveDetector(sd);
    }
    return true;
}

