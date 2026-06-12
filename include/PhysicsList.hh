#pragma once
// PhysicsList.hh
// INCL++/ABLA07 hadronic physics + NeutronHP for En < 20 MeV.
// No optical physics.

#include "G4VModularPhysicsList.hh"

class PhysicsList : public G4VModularPhysicsList
{
public:
    PhysicsList();
    ~PhysicsList() override = default;

    void ConstructParticle() override;
    void ConstructProcess() override;
    void SetCuts() override;
};
