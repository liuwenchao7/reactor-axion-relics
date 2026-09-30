#include "MuonGenerator.hh"
#include "BambooUtils.hh"

#include <G4Event.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <Randomize.hh>

using namespace CLHEP;

GeneratorRegister<MuonGenerator>
    MuonGenerator::reg("MuonGenerator");

MuonGenerator::MuonGenerator(const BambooParameters &pars)
    : BambooGenerator{pars}, gun{new G4ParticleGun} {
        // http://arxiv.org/abs/1606.06907
        // Energy and angular distributions of atmospheric muons at the Earth
        // constants in Shukla's function
        n = 3.01;
        E0 = 4.29 * GeV;
        epsilon = 854. * GeV;
        Rd = 174.;
        // http://arxiv.org/abs/1005.5332
        // Measurement of the charge ratio of atmospheric muons with the CMS detector
        charge_ratio = 1.2766;
        charge_thres = 1 / (1 + charge_ratio);
        theta_max = 0.616327;

        // Energy range in simulation
        E_low = 0.5 * GeV;
        E_up = 1000. * GeV;

        const auto &pmap = generatorParameters.getParameters();
        if (pmap.find("E_low") != pmap.end()) {
            E_low = BambooUtils::evaluate(
                generatorParameters.getParameter("E_low"));
        }
        if (pmap.find("E_up") != pmap.end()) {
            E_up = BambooUtils::evaluate(
                generatorParameters.getParameter("E_up"));
        }
        // Shielding size, muon will shower cling to the top
        if (pmap.find("shielding_x") != pmap.end()) {
            shielding_x = BambooUtils::evaluate(
                generatorParameters.getParameter("shielding_x"));
        } else {
            throw std::runtime_error("MuonGenerator: should provide `shielding_x`");
        }
        if (pmap.find("shielding_z") != pmap.end()) {
            shielding_z = BambooUtils::evaluate(
                generatorParameters.getParameter("shielding_z"));
        } else {
            throw std::runtime_error("MuonGenerator: should provide `shielding_z`");
        }
    }

G4double MuonGenerator::thetaSpectrum(G4double theta)
{
    G4double Dtheta = sqrt(Rd * Rd * cos(theta) * cos(theta) + 2 * Rd + 1) - Rd * cos(theta);
    G4double Dtheta_n = pow(Dtheta, -(n - 1)) * sin(theta);
    return Dtheta_n;
}

G4double MuonGenerator::energySpectrum(G4double energy)
{
    G4double IE = pow(E0 + energy, -n) / (1 + energy / epsilon);
    return IE;
}

G4double MuonGenerator::getRandomTheta() {
    G4double theta, rvs=1., upper=0.;
    G4double thetaSpectrum_max = thetaSpectrum(theta_max);
    while (rvs > upper)
    {
        theta = G4UniformRand() * M_PI_2;
        rvs = G4UniformRand();
        upper = thetaSpectrum(theta) / thetaSpectrum_max;
    }
    return theta;
}

G4double MuonGenerator::getRandomEnergy() {
    G4double energy, rvs=1., upper=0.;
    G4double energySpectrum_max = energySpectrum(E_low);
    while (rvs > upper)
    {
        energy = G4UniformRand() * (E_up - E_low) + E_low;
        rvs = G4UniformRand();
        upper = energySpectrum(energy) / energySpectrum_max;
    }
    return energy;
}

G4double MuonGenerator::getRandomPhi() {
    G4double phi = G4UniformRand() * M_PI * 2;
    return phi;
}

G4ParticleDefinition* MuonGenerator::getMuonCharge() {
    G4ParticleDefinition* particle;
    auto particleTable = G4ParticleTable::GetParticleTable();
    if (G4UniformRand() < charge_thres) {
        particle = particleTable->FindParticle("mu-");
    } else {
        particle = particleTable->FindParticle("mu+");
    }
    return particle;
}

void MuonGenerator::GeneratePrimaries(G4Event *anEvent) {
    // G4cout << "Cosmic ray muon generator used!" << G4endl;

    G4double posX, posY, posZ;
    G4double dirX, dirY, dirZ;
    G4double theta = getRandomTheta();
    G4double phi = getRandomPhi();
    G4double energy = getRandomEnergy();

    dirX = sin(theta) * cos(phi);
    dirY = sin(theta) * sin(phi);
    dirZ = -cos(theta);
    posX = (G4UniformRand() - 0.5) * shielding_x;
    posY = (G4UniformRand() - 0.5) * shielding_x;
    posZ = shielding_z / 2;

    auto particle = getMuonCharge();
    G4double mass = particle->GetPDGMass();

    gun->SetParticleDefinition(particle);
    gun->SetParticleEnergy(energy - mass);
    gun->SetParticlePosition(G4ThreeVector(posX, posY, posZ));
    gun->SetParticleMomentumDirection(G4ThreeVector(dirX, dirY, dirZ));

    gun->GeneratePrimaryVertex(anEvent);
}
