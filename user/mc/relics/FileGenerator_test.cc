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

    // =========== 添加调试信息 ===========
    std::cout << "\n=== FileGenerator构造函数开始 ===" << std::endl;
    
    if (!source_stream.is_open()){
        throw std::runtime_error("FileGenerator: Cannot open source file " + source);
    }
    
    source_stream.seekg(0,ios::beg);

    string particletype;
    double px, py, pz, x, y, z,energy;
    int runId,eventId;  
    
    if (source_stream >> event_number) {
        std::cout << "从文件读取到事件总数: " << event_number << std::endl;
        
        // 测试随机数生成
        std::cout << "\n=== 随机数生成测试（生成5个随机数） ===" << std::endl;
        for(int i = 0; i < 5; i++) {
            G4double rand_val = G4UniformRand();
            std::cout << "随机数 " << i << ": " << rand_val << std::endl;
        }
        std::cout << "=== 随机数生成测试结束 ===\n" << std::endl;
        
        for (int i = 0; i< event_number; ++ i) {
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
        std::cout << "成功加载 " << dataarray.size() << " 个粒子数据" << std::endl;
    } 
    else {
       throw std::runtime_error("FileGenerator: error reading data");
    }
    
    std::cout << "=== FileGenerator构造函数结束 ===\n" << std::endl;
}

G4int FileGenerator::getRandomLineNum(int LineNum) {
    G4double temp;
    G4int ReadLineNum;
    temp = G4UniformRand() * (LineNum);
    ReadLineNum = int(temp);
    
    // 添加调试输出，每100次调用打印一次
    static int function_call_count = 0;
    function_call_count++;
    
    if (function_call_count % 100 == 1) {  // 每100次打印一次，从第1次开始
        std::cout << "[DEBUG getRandomLineNum] 调用次数: " << function_call_count 
                  << ", 随机数: " << temp 
                  << ", 总行数: " << LineNum 
                  << ", 选择行: " << ReadLineNum << std::endl;
    }
    
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
    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(dir);
    fParticleGun->SetParticleEnergy(energy *keV);
    fParticleGun->GeneratePrimaryVertex(anEvent);
    
    counts++;
    
    // 修改打印条件：每100次调用打印一次（原来每100000次）
    if (counts % 100 == 0) {  
        // 获取当前时间
        chrono::system_clock::time_point now = chrono::system_clock::now();
        auto t = chrono::system_clock::to_time_t(now);
        auto localTm = gmtime(&t);
        char buff[32];
        strftime(buff, 32, "%H:%M:%S", localTm);
        
        // 打印详细信息
        std::cout << "\n[DEBUG GeneratePrimaries]" << std::endl;
        std::cout << "  事件计数: " << counts << std::endl;
        std::cout << "  粒子类型: " << typearray[N] << std::endl;
        std::cout << "  能量: " << energy << " keV" << std::endl;
        std::cout << "  位置: (" << x << ", " << y << ", " << z << ") mm" << std::endl;
        std::cout << "  动量方向: (" << px/energy << ", " << py/energy << ", " << pz/energy << ")" << std::endl;
        std::cout << "  选择的数据行: " << N << " / " << event_number << std::endl;
        std::cout << "  当前时间: " << buff << std::endl;
        std::cout << "[DEBUG结束]\n" << std::endl;
    }
    
    // 保持原有的每100000次打印的逻辑（如果需要可以保留）
    if ((counts-1)%100000 == 0) {
        std::cout << "input type: " << typearray[N]<<"ENERGY:"<<energy<<"x"<<x<<"y"<<y<<"z"<<z<<"input Line:"<<N<<endl;  
        chrono::system_clock::time_point now = chrono::system_clock::now();
        auto t = chrono::system_clock::to_time_t(now);
        auto localTm = gmtime(&t);
        char buff[32];
        strftime(buff, 32, "%H:%M:%S", localTm);
        std::cout << "time is:" << buff << endl;
    }
}