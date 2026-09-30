#pragma once
#include "G4UnitsTable.hh"
#include "BambooFactory.hh"
#include "BambooGenerator.hh"
#include <G4ParticleGun.hh>
#include <fstream>
#include <memory>
#include "G4SystemOfUnits.hh"
#include <G4ParticleTable.hh>
#include <string>

#include <chrono>
#include <vector>
using namespace CLHEP;
using namespace std;

class G4ParticleGun;
class G4Event;


class FileGenerator : public BambooGenerator {
  public:
    FileGenerator(const BambooParameters &pars);

    // ~FileGenerator();

    virtual void GeneratePrimaries(G4Event *anEvent);

    static GeneratorRegister<FileGenerator> reg;
    int start0 = 0;
    int i =1;
    int event_number;

    vector<vector<double>> dataarray; // energy, px, py, pz, x, y, z 
    vector<string> typearray; //particle type
    int counts;


  
  
 
  private:
    
    

    G4int getRandomLineNum(int);

    std::unique_ptr<G4ParticleGun> fParticleGun;


    
    std::string source;
    std::ifstream source_stream;
};
