#pragma once
// RunAction.hh — books histograms and ntuple, opens/closes ROOT file

#include "G4UserRunAction.hh"
#include "G4Run.hh"

class RunAction : public G4UserRunAction
{
public:
    RunAction();
    ~RunAction() override = default;

    void BeginOfRunAction(const G4Run* run) override;
    void EndOfRunAction(const G4Run* run) override;
};
