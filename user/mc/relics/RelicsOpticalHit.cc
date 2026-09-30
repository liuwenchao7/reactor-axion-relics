#include "RelicsOpticalHit.hh"

G4ThreadLocal G4Allocator<RelicsOpticalHit> *RelicsOpticalHitAllocator = nullptr;

const RelicsOpticalHit &RelicsOpticalHit::
operator=(const RelicsOpticalHit &right) {
    pmt_number = right.pmt_number;
    pos = right.pos;
    energy = right.energy;
    creatorProcess = right.creatorProcess;
    t_global = right.t_global;
    t_local = right.t_local;
    t_proper = right.t_proper;
    velocity = right.velocity;
    return *this;
}

G4int RelicsOpticalHit::operator==(const RelicsOpticalHit &right) const {
    return ((this == &right) ? 1 : 0);
}
