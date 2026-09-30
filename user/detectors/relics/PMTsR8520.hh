#pragma once

#include "BambooControl.hh"
#include "BambooDetector.hh"

#include <G4LogicalVolume.hh>

struct PMTLVs
{
   G4LogicalVolume *CombineLV;
   G4LogicalVolume *WindowLV;
   G4LogicalVolume *PhotoCathodeLV;
   G4LogicalVolume *CasingLV;
};

class PMTsR8520 {

  public:
    PMTsR8520(const std::string &n, const BambooParameters &pars);

    PMTLVs construct(const BambooParameters &global_pars, BambooDetector *parent);

  private:
    // define additional parameters here
    BambooParameters parameters;

    G4double WindowWidth;
    G4double WindowThickness;
    G4double PhotoCathodeWidth;
    G4double PhotoCathodeThickness;
    G4double CasingWidth;
    G4double CasingHeight;
    G4double CasingThickness;
    G4double BaseThickness;
    G4double PMTToPMTBase;

    G4LogicalVolume *CombineLogicalVolume;
    G4LogicalVolume *WindowLogicalVolume;
    G4LogicalVolume *CasingLogicalVolume;
    G4LogicalVolume *InteriorLogicalVolume;
    G4LogicalVolume *PhotoCathodeLogicalVolume;
    G4LogicalVolume *BaseLogicalVolume;

    G4VPhysicalVolume *WindowPhysicalVolume;
    G4VPhysicalVolume *CasingPhysicalVolume;
    G4VPhysicalVolume *InteriorPhysicalVolume;
    G4VPhysicalVolume *PhotoCathodePhysicalVolume;
    G4VPhysicalVolume *BasePhysicalVolume;
};

