#pragma once
/*===========================================================
  PrimaryGeneratorAction.hh

  Atmospheric neutron spectrum (Gordon et al. 2004) with
  importance sampling: the high-energy component (>20 MeV)
  is oversampled by factor ~50 relative to the true spectrum.
  Each event carries an event weight stored thread-locally.
  
  EventAction retrieves both energy and weight via:
    PrimaryGeneratorAction::GetCurrentEnergy()
    PrimaryGeneratorAction::GetCurrentWeight()
  
  The weight corrects for the biased sampling when filling
  histograms — multiply all histogram fills by GetCurrentWeight().
===========================================================*/

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "G4Types.hh"

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* event) override;

    static G4double GetCurrentEnergy() { return fgCurrentEnergy; }
    static G4double GetCurrentWeight() { return fgCurrentWeight; }

private:
    G4ParticleGun* fParticleGun = nullptr;

    static G4ThreadLocal G4double fgCurrentEnergy;
    static G4ThreadLocal G4double fgCurrentWeight;

    G4double fSourceZ    = 26.0*CLHEP::cm;  // upstream of calorimeter front face

    G4double      SampleEnergy(G4double& weight) const;
    G4ThreeVector SampleHemisphere() const;
};
