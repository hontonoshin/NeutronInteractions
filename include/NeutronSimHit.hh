#pragma once
// NeutronSimHit.hh
// Lightweight hit: stores Edep and secondary particle info per step.

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4Allocator.hh"

#include <vector>
#include <string>

struct SecondaryInfo {
    G4int    pdg;
    G4double kineticEnergy;  // MeV
    std::string name;
};

class NeutronSimHit : public G4VHit
{
public:
    NeutronSimHit() = default;
    ~NeutronSimHit() override = default;

    // Allocator (performance)
    void* operator new(size_t);
    void  operator delete(void* hit);

    // Data members (public for simplicity)
    G4double              edep          = 0.0;   // MeV
    G4double              stepLength    = 0.0;   // mm
    G4ThreeVector         position;
    G4int                 copyNo        = -1;
    G4int                 trackID       = -1;
    G4int                 parentID      = -1;
    G4int                 pdg           = 0;
    G4double              preKineticEnergy = 0.0;
    std::string           processName;
    std::vector<SecondaryInfo> secondaries;  // created at this step
};

using NeutronSimHitsCollection = G4THitsCollection<NeutronSimHit>;

extern G4ThreadLocal G4Allocator<NeutronSimHit>* gNeutronSimHitAllocator;

inline void* NeutronSimHit::operator new(size_t)
{
    if (!gNeutronSimHitAllocator)
        gNeutronSimHitAllocator = new G4Allocator<NeutronSimHit>;
    return gNeutronSimHitAllocator->MallocSingle();
}
inline void NeutronSimHit::operator delete(void* hit)
{
    gNeutronSimHitAllocator->FreeSingle(static_cast<NeutronSimHit*>(hit));
}
