#pragma once

#include <G4Allocator.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include <string>

using std::string;

class RelicsOpticalHit : public G4VHit {
  public:
    RelicsOpticalHit() = default;
    ~RelicsOpticalHit() = default;

    const RelicsOpticalHit &operator=(const RelicsOpticalHit &right);
    G4int operator==(const RelicsOpticalHit &right) const;

    inline void *operator new(size_t);

    inline void operator delete(void *aHit);

    void Draw() {}
    void Print() {}

    void setPMTNumber(int n) { pmt_number = n; }
    int getPMTNumber() const { return pmt_number; }

    void setEnergy(double v) { energy = v; }
    double getEnergy() const { return energy; }

    void setPosition(const G4ThreeVector &v) { pos = v; }
    const G4ThreeVector &getPosition() const { return pos; }

    void setCreatorProcess(const string &s) { creatorProcess = s; }
    const string &getCreatorProcess() const { return creatorProcess; }

    void setGlobalTime(double v) { t_global = v; }
    double getGlobalTime() const { return t_global; }

    void setLocalTime(double v) { t_local = v; }
    double getLocalTime() const { return t_local; }

    void setProperTime(double v) { t_proper = v; }
    double getProperTime() const { return t_proper; }

    void setVelocity(double v) { velocity = v; }
    double getVelocity() const { return velocity; }

    void setSourcePos(const G4ThreeVector &v) { source_pos = v; }
    const G4ThreeVector &getSourcePos() const { return source_pos; }

    void setParent(const std::string &s) { parent = s; }
    const std::string &getParent() { return parent; }

    void setParentId(int n) { parentId = n; }
    int getParentId() { return parentId; }

  private:
    int pmt_number = -1;
    G4ThreeVector pos = {0, 0, 0};
    double energy = 0;
    string creatorProcess = "";
    double t_global = 0;
    double t_local = 0;
    double t_proper = 0;
    double velocity = 0;

    G4ThreeVector source_pos = {0, 0, 0};
    string parent = "";
    int parentId = -1;
};

using RelicsOpticalHitsCollection = G4THitsCollection<RelicsOpticalHit>;

extern G4ThreadLocal G4Allocator<RelicsOpticalHit> *RelicsOpticalHitAllocator;

inline void *RelicsOpticalHit::operator new(size_t) {
    if (!RelicsOpticalHitAllocator)
        RelicsOpticalHitAllocator = new G4Allocator<RelicsOpticalHit>;
    return (void *)RelicsOpticalHitAllocator->MallocSingle();
}

inline void RelicsOpticalHit::operator delete(void *aHit) {
    RelicsOpticalHitAllocator->FreeSingle((RelicsOpticalHit *)aHit);
}
