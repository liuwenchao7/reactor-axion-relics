#include "FileGenerator.hh"
// #include "FileGeneratorMessenger.hh"
#include <string>
#include "BambooUtils.hh"
#include "G4UnitsTable.hh"
#include <G4Event.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include "G4SystemOfUnits.hh"
#include <Randomize.hh>
#include <iostream>
#include <ctime>
#include <chrono> 
#include <vector>

#include <cstdio>

using namespace CLHEP;
using namespace std;



GeneratorRegister<FileGenerator> FileGenerator::reg("FileGenerator");
vector<vector<double>> dataarray; // energy, px, py, pz, x, y, z 
vector<string> typearray; //particle type

int event_number;
int counts = 0;


FileGenerator::FileGenerator(const BambooParameters &pars)
    : BambooGenerator{pars},
      source{pars.getParameter("source")}, source_stream{source},fParticleGun{new G4ParticleGun} {

    if (!source_stream.is_open()){
        throw std::runtime_error("FileGenerator: Cannot open source file " +
                                 source);
                                 

                                 }
    source_stream.seekg(0,ios::beg);


    string particletype;
    double px, py, pz, x, y, z,energy;
    int runId,eventId;  
    
    
    if (source_stream >> event_number) {
             
       for (int i = 0; i< event_number; ++ i)
        {
        vector<double> temp;
        source_stream >> runId >> eventId >> energy >> particletype >> px >> py >> pz >> x >> y >> z ;
        temp.push_back(energy);
        temp.push_back(px);
        temp.push_back(py);
        temp.push_back(pz);
        
        temp.push_back(x);
        temp.push_back(y);
        temp.push_back(z);
        dataarray.push_back(temp);
        temp.clear();
        typearray.push_back(particletype);

        }

        
    } 
    else {
       throw std::runtime_error("FileGenerator: error reading data");
    }


      }


// FileGenerator::~FileGenerator()
// {
//     delete fParticleGun;
// }


G4int FileGenerator::getRandomLineNum(int LineNum) {
    G4double temp;
    G4int ReadLineNum;
    temp = G4UniformRand() * (LineNum);
    ReadLineNum = int(temp);
    // std::cout << "reading Line: " << ReadLineNum
    //               << std::endl;
    return ReadLineNum;
}





void FileGenerator::GeneratePrimaries(G4Event *anEvent) {

    using CLHEP::keV;
    using CLHEP::mm;

     
    
    double px, py, pz, x, y, z,energy;



    int N;
    N = getRandomLineNum(event_number);
    energy = dataarray[N][0];
    px = dataarray[N][1];
    py = dataarray[N][2];
    pz = dataarray[N][3];
    x = dataarray[N][4];
    y = dataarray[N][5];
    z = dataarray[N][6];
    


    G4ThreeVector pos(x*mm,y*mm,z*mm);
    G4ThreeVector dir(px/energy,py/energy,pz/energy);
            
            

    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();

    G4ParticleDefinition *particle = particleTable->FindParticle(typearray[N]);


    //Fire ParticleGun
    fParticleGun->SetParticleDefinition(particle);
    // G4double mass = particle->GetPDGMass();
    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(dir);
    fParticleGun->SetParticleEnergy(energy *keV);
    // std::cout << "mass:" << mass << endl;
    fParticleGun->SetParticleDefinition(particle);

    fParticleGun->GeneratePrimaryVertex(anEvent);
    counts ++ ;
    //print information ,check if there's any error 
    if ((counts-1)%100000 ==0 ){
        
        
        
        std::cout << "input type: " << typearray[N]<<"ENERGY:"<<energy<<"x"<<x<<"y"<<y<<"z"<<z<<"input Line:"<<N<<endl;  

        chrono::system_clock::time_point now = chrono::system_clock::now();
        auto t = chrono::system_clock::to_time_t(now);
	    auto localTm = gmtime(&t);
	    char buff[32];
	    strftime(buff, 32, "%H:%M:%S", localTm);
	    std::cout << "time is:" << buff << endl;
    }


            

        
        
        
        
     
}
 