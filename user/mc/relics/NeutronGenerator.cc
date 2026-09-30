#include "NeutronGenerator.hh"
#include "BambooUtils.hh"

#include <G4Event.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <Randomize.hh>
#include <random>

using namespace CLHEP;

GeneratorRegister<NeutronGenerator>
    NeutronGenerator::reg("NeutronGenerator");

NeutronGenerator::NeutronGenerator(const BambooParameters &pars)
    : BambooGenerator{pars}, gun{new G4ParticleGun}{
        const auto &pmap = generatorParameters.getParameters();

        // Shielding size, CRN will shower cling to the top
        if (pmap.find("shielding_x") != pmap.end()) {
            shielding_x = BambooUtils::evaluate(
                generatorParameters.getParameter("shielding_x"));
        } else {
            throw std::runtime_error("NeutronGenerator: should provide `shielding_x`");
        }
        if (pmap.find("shielding_z") != pmap.end()) {
            shielding_z = BambooUtils::evaluate(
                generatorParameters.getParameter("shielding_z"));
        } else {
            throw std::runtime_error("NeutronGenerator: should provide `shielding_z`");
        }
}

void NeutronGenerator::NeutronLoad(const std::string &energySpectrumFile, const std::string &angularDistributionFile)
{
    LoadEnergySpectrum(energySpectrumFile);
    LoadAngularDistribution(angularDistributionFile);
}


void NeutronGenerator::LoadEnergySpectrum(const std::string& energySpectrumFile) {
    energyBins.clear();
    energySpectrum.clear();
    std::ifstream file(energySpectrumFile);
    if (!file.is_open()) {
        std::cerr << "Error: Failed to open energy spectrum file: " << energySpectrumFile << std::endl;
        return;
    }

    double bin, value;
    while (file >> bin >> value) {
        energyBins.push_back(bin);
        energySpectrum.push_back(value);
    }
    int n = energySpectrum.size();
    
    for (int i = 0; i < n - 1; ++i) {
        double diff = energyBins[i + 1] - energyBins[i];
        double result = diff * energySpectrum[i];
        energyValue.push_back(result);
    }
    energyValue.push_back(0);
    file.close();
    // // 输出读取成功 测试能谱
    // std::cout << "读取成功" << std::endl;

    // // 打印 angleBins 内容
    // std::cout << "energyValue 中的内容：" << std::endl;
    // for (const auto& test : energyValue) {
    //     std::cout << test << " ";
    // }
    // std::cout << std::endl;
}


void NeutronGenerator::LoadAngularDistribution(const std::string& angularDistributionFile) {
    angleBins.clear();
    angularDistribution.clear();
    std::ifstream file(angularDistributionFile);
    if (!file.is_open()) {
        std::cerr << "Error: Failed to open angular distribution file: " << angularDistributionFile << std::endl;
        return;
    }

    double Theta, value;
    while (file >> Theta >> value) {
        angleBins.push_back(Theta);
        angularDistribution.push_back(value);
    }

    file.close();

    // // 输出读取成功
    // std::cout << "读取成功" << std::endl;

    // // 打印 angleBins 内容
    // std::cout << "angleBins 中的内容：" << std::endl;
    // for (const auto& angle : angleBins) {
    //     std::cout << angle << " ";
    // }
    // std::cout << std::endl;
}

G4ParticleDefinition* NeutronGenerator::getneutron() {
    G4ParticleDefinition* particle;
    auto particleTable = G4ParticleTable::GetParticleTable();
    particle = particleTable->FindParticle("neutron");
    // if (particle) {
    //     std::cout << "Particle Name: " << particle->GetParticleName() << std::endl;
    //     std::cout << "Particle Mass: " << particle->GetPDGMass() << std::endl;
    // } else {
    //     std::cout << "Particle not found!" << std::endl;
    // }
    return particle;
}

void NeutronGenerator::GeneratePrimaries(G4Event *anEvent) {
    LoadEnergySpectrum("./data/neutronE.txt");
    LoadAngularDistribution("./data/neutronA.txt");
    double energy, directionX, directionY, directionZ;
    std::random_device rd;
    std::mt19937 gen(rd());

    std::piecewise_constant_distribution<double> disEnergy(
        energyBins.begin(), energyBins.end(), energyValue.begin()
    );
    energy = disEnergy(gen);
    // // 打开文件流 检验能量分布
    // std::ofstream outputFile("Energy.txt", std::ios::app);

    // // 生成并写入数据到文件
    // outputFile << energy << '\n';  // 每个数据后添加换行符

    // // 关闭文件流
    // outputFile.close();

    // // Generate a random energy based on the spectrum
    // double randomEnergy = G4UniformRand() * (energyBins.back() - energyBins.front()) + energyBins.front();
    // energy = randomEnergy;

    // // Generate a random direction based on the angular distribution
    // double randomAngle = G4UniformRand() * (angleBins.back() - angleBins.front()) + angleBins.front();
    // double cosTheta = randomAngle;

    std::piecewise_constant_distribution<double> disAngle(
        angleBins.begin(), angleBins.end(), angularDistribution.begin()
    );
    double cosTheta = disAngle(gen);

    // // 打开文件流 检验角度分布
    // std::ofstream outputFile1("Angle.txt", std::ios::app);

    // // 生成并写入数据到文件
    // outputFile1 << cosTheta << '\n';  // 每个数据后添加换行符

    // // 关闭文件流
    // outputFile1.close();

    // Convert cosTheta to direction components
    double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
    double phi = G4UniformRand() * M_PI * 2;
    directionX = sinTheta * std::cos(phi);
    directionY = sinTheta * std::sin(phi);
    directionZ = -cosTheta;

    G4double posX, posY, posZ;
    posX = (G4UniformRand() - 0.5) * shielding_x;
    posY = (G4UniformRand() - 0.5) * shielding_x;
    posZ = shielding_z / 2;

    auto particle = getneutron();
    G4double mass = particle->GetPDGMass();

    gun->SetParticleDefinition(particle);
    gun->SetParticleEnergy(energy);
    gun->SetParticlePosition(G4ThreeVector(posX, posY, posZ));
    gun->SetParticleMomentumDirection(G4ThreeVector(directionX, directionY, directionZ));

    gun->GeneratePrimaryVertex(anEvent);
}

