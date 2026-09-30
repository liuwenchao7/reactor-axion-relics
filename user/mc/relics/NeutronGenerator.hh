#pragma once

#include "BambooFactory.hh"
#include "BambooGenerator.hh"

#include <G4ParticleGun.hh>

#include <memory>

class G4ParticleGun;
class G4Event;

class NeutronGenerator : public BambooGenerator {
  public:
    NeutronGenerator(const BambooParameters &pars);

    ~NeutronGenerator() = default;

    virtual void GeneratePrimaries(G4Event *anEvent);

    static GeneratorRegister<NeutronGenerator> reg;

  private:
    G4double shielding_x;
    G4double shielding_z;
    void NeutronLoad(const std::string& energySpectrumFile, const std::string& angularDistributionFile);
    std::vector<double> energyBins;
    std::vector<double> energySpectrum;
    std::vector<double> energyValue;
    std::vector<double> angleBins;
    std::vector<double> angularDistribution;

    G4ParticleDefinition* getneutron();
    void LoadEnergySpectrum(const std::string& energySpectrumFile);
    void LoadAngularDistribution(const std::string& angularDistributionFile);

    std::unique_ptr<G4ParticleGun> gun;
};
