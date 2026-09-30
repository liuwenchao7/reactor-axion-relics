#pragma once

#include "BambooMaterial.hh"
#include "BambooFactory.hh"

#include <G4Version.hh>
#include <G4OpticalSurface.hh>

struct BorderPVs
{
    G4VPhysicalVolume *vol1;
    G4VPhysicalVolume *vol2;
};

class Material : public BambooMaterial {

  public:
    Material() = default;

    virtual void defineMaterials(const BambooParameters &pars);

    static MaterialRegister<Material> reg;

    void AddLogicalSurface(
      G4VPhysicalVolume* mainPV,
      G4VPhysicalVolume* containerPV,
      const BambooParameters &parameters,
      const G4String &prefix = ""
    );

    void AddLogicalBorderSurface(
      G4VPhysicalVolume* mainPV,
      G4VPhysicalVolume* containerPV,
      G4OpticalSurface* opSurface,
      const BambooParameters &parameters,
      const G4String &prefix
    );

    void AddLogicalSkinSurface(
      G4LogicalVolume* mainLV,
      G4OpticalSurface* opSurface,
      const BambooParameters &parameters,
      const G4String &prefix
    );

    G4OpticalSurface* EvalOpParameters(
      const BambooParameters &parameters,
      G4VPhysicalVolume* mainPV,
      const G4String &prefix
    );

    void SetProperties(
      const BambooParameters &parameters,
      const std::string &material,
      G4MaterialPropertiesTable *table
    );

    void SetVectorProperty(
      G4MaterialPropertiesTable *table,
      const std::string &property,
      G4String newValue
    );

    void SetConstProperty(
      G4MaterialPropertiesTable *table,
      const std::string &property,
      G4String newValue
    );

  private:
    std::map<std::string, G4OpticalSurfaceModel> Models = {
      {"", unified},  // default value
      {"glisur", glisur},
      {"unified", unified}
    };

    std::map<std::string, G4OpticalSurfaceFinish> Finishes = {
      {"", ground},  // default value
      {"polished", polished},
      {"ground", ground}
    };

    std::map<std::string, G4SurfaceType> Types = {
      {"", dielectric_dielectric},  // default value
      {"dielectric_metal", dielectric_metal},
      {"dielectric_dielectric", dielectric_dielectric}
    };

    // energy in eV!
    const std::set<std::string> vector_properties = {
        "RINDEX",
        "ABSLENGTH",                // mm
        "RAYLEIGH",                 // mm
#if G4VERSION_NUMBER>=1100
        "SCINTILLATIONCOMPONENT1",
        "SCINTILLATIONCOMPONENT2",
#else
        "FASTCOMPONENT",
        "SLOWCOMPONENT",
#endif
        "EFFICIENCY",
        "REFLECTIVITY",
        "SPECULARLOBECONSTANT",
        "SPECULARSPIKECONSTANT",
        "BACKSCATTERCONSTANT",
        "TRANSMITTANCE"
    };

    const std::set<std::string> const_properties = {
        "SCINTILLATIONYIELD",          // 1/MeV
        "RESOLUTIONSCALE",
#if G4VERSION_NUMBER>=1100
        "SCINTILLATIONTIMECONSTANT1",  // ns
        "SCINTILLATIONTIMECONSTANT2",  // ns
        "SCINTILLATIONYIELD1",
        "SCINTILLATIONYIELD2",
#else
        "FASTTIMECONSTANT",            // ns
        "SLOWTIMECONSTANT",            // ns
        "YIELDRATIO",
#endif
    };
};

