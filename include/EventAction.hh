#pragma once
// EventAction.hh

#include "G4UserEventAction.hh"
#include "G4Event.hh"

class EventAction : public G4UserEventAction
{
public:
    EventAction();
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

    // Primary energy now retrieved via PrimaryGeneratorAction::GetCurrentEnergy()
    // (thread-local) — no setter needed here.

private:
    G4int fCrystalHCID = -1;
    G4int fBiHCID      = -1;
};
