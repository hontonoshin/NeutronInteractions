#pragma once
// CrystalSD.hh — sensitive detector for GAGG bulk

#include "G4VSensitiveDetector.hh"
#include "NeutronSimHit.hh"

class CrystalSD : public G4VSensitiveDetector
{
public:
    CrystalSD(const G4String& name, const G4String& hcName);
    ~CrystalSD() override = default;

    void   Initialize(G4HCofThisEvent* hce) override;
    G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override;
    void   EndOfEvent(G4HCofThisEvent* hce) override;

private:
    NeutronSimHitsCollection* fHitsCollection = nullptr;
    G4int fHCID = -1;
};
