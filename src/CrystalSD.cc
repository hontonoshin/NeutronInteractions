/*===========================================================
  CrystalSD.cc
  Records every step with Edep > 0 in the GAGG crystal.
  Captures:
    - deposited energy per step
    - incident particle PDG + pre-step KE
    - process name that caused the step
    - all secondaries created at this step (PDG + KE)
  This lets us reconstruct: what entered the crystal, what
  process fired, and what fragments/gammas were produced.
===========================================================*/

#include "CrystalSD.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"

CrystalSD::CrystalSD(const G4String& name, const G4String& hcName)
    : G4VSensitiveDetector(name)
{
    collectionName.insert(hcName);
}

void CrystalSD::Initialize(G4HCofThisEvent* hce)
{
    fHitsCollection = new NeutronSimHitsCollection(
        SensitiveDetectorName, collectionName[0]);
    if (fHCID < 0)
        fHCID = G4SDManager::GetSDMpointer()
                    ->GetCollectionID(collectionName[0]);
    hce->AddHitsCollection(fHCID, fHitsCollection);
}

G4bool CrystalSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    const G4double edep = step->GetTotalEnergyDeposit();
    // Record step only if it deposits energy OR if it's a neutron
    // (neutrons can scatter without depositing — we still want to
    //  track what processes fired on them inside the crystal).
    const G4Track* track = step->GetTrack();
    const G4bool isNeutron = (track->GetDefinition()->GetPDGEncoding() == 2112);
    if (edep <= 0.0 && !isNeutron) return false;

    auto* hit = new NeutronSimHit();
    hit->edep             = edep;
    hit->stepLength       = step->GetStepLength();
    hit->position         = step->GetPostStepPoint()->GetPosition();
    hit->copyNo           = step->GetPreStepPoint()
                                ->GetTouchable()->GetCopyNumber(1); // Bi shell copy
    hit->trackID          = track->GetTrackID();
    hit->parentID         = track->GetParentID();
    hit->pdg              = track->GetDefinition()->GetPDGEncoding();
    hit->preKineticEnergy = step->GetPreStepPoint()->GetKineticEnergy();

    const G4VProcess* proc = step->GetPostStepPoint()->GetProcessDefinedStep();
    hit->processName = proc ? proc->GetProcessName() : "Transportation";

    // Collect secondaries created at this step
    const std::vector<const G4Track*>* secs = step->GetSecondaryInCurrentStep();
    if (secs) {
        for (const auto* sec : *secs) {
            SecondaryInfo si;
            si.pdg           = sec->GetDefinition()->GetPDGEncoding();
            si.kineticEnergy = sec->GetKineticEnergy() / MeV;
            si.name          = sec->GetDefinition()->GetParticleName();
            hit->secondaries.push_back(si);
        }
    }

    fHitsCollection->insert(hit);
    return true;
}

void CrystalSD::EndOfEvent(G4HCofThisEvent*) {}
