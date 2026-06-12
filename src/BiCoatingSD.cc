/*===========================================================
  BiCoatingSD.cc
  Records interactions in the 5-µm Bi-209 shell.
  
  Key physics captured here:
    - Neutron elastic/inelastic scattering off Bi-209
    - INCL++/ABLA07 spallation/fragmentation products
      (heavy recoils, protons, alphas, light ions)
    - Energy deposited by fragments before they enter the crystal
  
  Every step (even zero-edep neutron steps) is recorded so we
  can see *what process* fires on the Bi layer.
===========================================================*/

#include "BiCoatingSD.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"

BiCoatingSD::BiCoatingSD(const G4String& name, const G4String& hcName)
    : G4VSensitiveDetector(name)
{
    collectionName.insert(hcName);
}

void BiCoatingSD::Initialize(G4HCofThisEvent* hce)
{
    fHitsCollection = new NeutronSimHitsCollection(
        SensitiveDetectorName, collectionName[0]);
    if (fHCID < 0)
        fHCID = G4SDManager::GetSDMpointer()
                    ->GetCollectionID(collectionName[0]);
    hce->AddHitsCollection(fHCID, fHitsCollection);
}

G4bool BiCoatingSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    const G4double edep   = step->GetTotalEnergyDeposit();
    const G4Track* track  = step->GetTrack();
    const G4int    pdg    = track->GetDefinition()->GetPDGEncoding();

    // For Bi layer: record ALL steps of primary neutron (even if no edep)
    // and all steps with energy deposition (secondaries, fragments, etc.)
    const bool isPrimaryNeutron = (pdg == 2112 && track->GetParentID() == 0);
    if (edep <= 0.0 && !isPrimaryNeutron) return false;

    auto* hit = new NeutronSimHit();
    hit->edep             = edep;
    hit->stepLength       = step->GetStepLength();
    hit->position         = step->GetPreStepPoint()->GetPosition();
    hit->copyNo           = step->GetPreStepPoint()
                                ->GetTouchable()->GetCopyNumber(); // Bi copy#
    hit->trackID          = track->GetTrackID();
    hit->parentID         = track->GetParentID();
    hit->pdg              = pdg;
    hit->preKineticEnergy = step->GetPreStepPoint()->GetKineticEnergy();

    const G4VProcess* proc = step->GetPostStepPoint()->GetProcessDefinedStep();
    hit->processName = proc ? proc->GetProcessName() : "Transportation";

    // Secondaries: the interesting fragmentation products appear here
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

void BiCoatingSD::EndOfEvent(G4HCofThisEvent*) {}
