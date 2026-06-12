/*===========================================================
  EventAction.cc
  
  Per-event aggregation and ROOT output.
  
  What we fill:
  
  Histograms (TH1D via G4Analysis):
    H0  : Total Edep in GAGG per event [MeV]
    H1  : Total Edep in Bi per event [MeV]
    H2  : Edep ratio Bi / (Bi + GAGG)
    H3  : Primary neutron energy [GeV]
    H4  : Number of secondaries created in Bi per event
    H5  : Number of secondaries created in GAGG per event
    H6  : KE spectrum of ALL secondaries from Bi [MeV]
    H7  : KE spectrum of ALL secondaries from GAGG [MeV]
    H8  : Secondary PDG code histogram (Bi) — fragmentation fingerprint
    H9  : Secondary PDG code histogram (GAGG)
  
  Ntuple (one row per event):
    col 0 : event ID
    col 1 : primary neutron energy [GeV]
    col 2 : total Edep in GAGG [MeV]
    col 3 : total Edep in Bi [MeV]
    col 4 : number of hits in GAGG SD
    col 5 : number of hits in Bi SD
    col 6 : number of unique secondary types in Bi (fragmentation richness)
    col 7 : max secondary KE in Bi [MeV] (largest fragment energy)
    col 8 : flag: did a hadronic inelastic process fire in Bi? (0/1)
    col 9 : flag: did a hadronic inelastic process fire in GAGG? (0/1)
    col 10: Edep in GAGG from heavy ions (A>4) [MeV]
    col 11: Edep in GAGG from protons/neutrons [MeV]
    col 12: Edep in GAGG from e±/gamma [MeV]
===========================================================*/

#include "EventAction.hh"
#include "NeutronSimHit.hh"
#include "PrimaryGeneratorAction.hh"   // for GetCurrentEnergy()

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include "G4Event.hh"

#include <set>
#include <cmath>

// ── Helper: returns true if process name indicates hadronic inelastic ────────
// Covers all Geant4 naming variants:
//   "neutronInelastic", "NeutronInelastic", "protonInelastic",
//   "pi-Inelastic", "hBertiniInelastic", "hFritiofInelastic", etc.
static bool IsHadronicInelastic(const std::string& proc)
{
    if (proc.find("Inelastic") != std::string::npos) return true;
    if (proc.find("inelastic") != std::string::npos) return true;
    if (proc.find("INCLXX")    != std::string::npos) return true;
    if (proc.find("Fritiof")   != std::string::npos) return true;
    if (proc.find("Bertini")   != std::string::npos) return true;
    return false;
}

// ── Helper: classify PDG into broad particle category ────────────────────────
// Returns: 0=other/unknown, 1=neutron, 2=proton, 3=light ion (A≤4, Z≥1),
//          4=heavy ion (A>4), 5=e±/gamma, 6=pion/kaon, 7=alpha
static int ClassifyPDG(G4int pdg)
{
    if (pdg == 2112) return 1;           // neutron
    if (pdg == 2212) return 2;           // proton
    if (pdg == 11 || pdg == -11) return 5;
    if (pdg == 22)  return 5;            // gamma
    if (pdg == 211 || pdg == -211 || pdg == 111) return 6; // pion
    if (pdg == 321 || pdg == -321) return 6;  // kaon
    if (pdg == 1000020040) return 7;     // alpha
    // Ion PDG: 10LZZZAAAI — check if ion
    if (std::abs(pdg) > 1000000000) {
        G4int A = (std::abs(pdg) / 10) % 1000;
        if (A > 4) return 4;  // heavy ion
        return 3;              // light ion (d, t, He-3)
    }
    return 0;
}

EventAction::EventAction() : G4UserEventAction() {}

void EventAction::BeginOfEventAction(const G4Event*) {}

void EventAction::EndOfEventAction(const G4Event* event)
{
    // Retrieve primary energy from thread-local set by PrimaryGeneratorAction
    const G4double fPrimaryEkin = PrimaryGeneratorAction::GetCurrentEnergy();

    auto* sdMgr = G4SDManager::GetSDMpointer();
    auto* hce   = event->GetHCofThisEvent();
    if (!hce) return;

    // Cache HC IDs on first call
    if (fCrystalHCID < 0)
        fCrystalHCID = sdMgr->GetCollectionID("CrystalSD/CrystalHitsCollection");
    if (fBiHCID < 0)
        fBiHCID = sdMgr->GetCollectionID("BiCoatingSD/BiHitsCollection");

    auto* crystalHC = static_cast<NeutronSimHitsCollection*>(
                          hce->GetHC(fCrystalHCID));
    auto* biHC      = static_cast<NeutronSimHitsCollection*>(
                          hce->GetHC(fBiHCID));

    // ── Accumulate GAGG hits ─────────────────────────────────────────────────
    G4int    evtID             = event->GetEventID();  // needed by ntuple 1
    G4double edepGAGG         = 0.0;
    G4double edepGAGG_heavyIon = 0.0;
    G4double edepGAGG_nP      = 0.0;
    G4double edepGAGG_em      = 0.0;
    G4int    nSecGAGG         = 0;
    G4bool   hadronicInGAGG   = false;

    if (crystalHC) {
        G4int nhits = crystalHC->entries();
        auto* am = G4AnalysisManager::Instance();

        for (G4int i = 0; i < nhits; ++i) {
            auto* hit = (*crystalHC)[i];
            edepGAGG += hit->edep;

            // Classify Edep by particle type
            int cat = ClassifyPDG(hit->pdg);
            if (cat == 4) edepGAGG_heavyIon += hit->edep;
            else if (cat == 1 || cat == 2) edepGAGG_nP += hit->edep;
            else if (cat == 5) edepGAGG_em += hit->edep;

            if (IsHadronicInelastic(hit->processName)) hadronicInGAGG = true;

            // Secondaries — fill histos and ntuple 1
            for (const auto& sec : hit->secondaries) {
                ++nSecGAGG;
                am->FillH1(7, sec.kineticEnergy);
                am->FillH1(9, (G4double)sec.pdg);
                // Ntuple 1: per-secondary record (volume=1 → GAGG)
                am->FillNtupleIColumn(1, 0, evtID);
                am->FillNtupleIColumn(1, 1, 1);
                am->FillNtupleDColumn(1, 2, (G4double)sec.pdg);
                am->FillNtupleDColumn(1, 3, sec.kineticEnergy);
                am->FillNtupleIColumn(1, 4, hit->pdg);
                am->FillNtupleDColumn(1, 5, fPrimaryEkin / GeV);
                am->AddNtupleRow(1);
            }
        }
    }

    // ── Accumulate Bi hits ───────────────────────────────────────────────────
    G4double edepBi       = 0.0;
    G4double maxSecKE_Bi  = 0.0;
    G4int    nSecBi       = 0;
    G4bool   hadronicInBi = false;
    std::set<G4int> secTypesInBi;

    if (biHC) {
        G4int nhits = biHC->entries();
        auto* am = G4AnalysisManager::Instance();

        for (G4int i = 0; i < nhits; ++i) {
            auto* hit = (*biHC)[i];
            edepBi += hit->edep;

            if (IsHadronicInelastic(hit->processName)) hadronicInBi = true;

            for (const auto& sec : hit->secondaries) {
                ++nSecBi;
                secTypesInBi.insert(sec.pdg);
                if (sec.kineticEnergy > maxSecKE_Bi)
                    maxSecKE_Bi = sec.kineticEnergy;
                am->FillH1(6, sec.kineticEnergy);
                am->FillH1(8, (G4double)sec.pdg);
                // Ntuple 1: per-secondary record (volume=0 → Bi)
                am->FillNtupleIColumn(1, 0, evtID);
                am->FillNtupleIColumn(1, 1, 0);
                am->FillNtupleDColumn(1, 2, (G4double)sec.pdg);
                am->FillNtupleDColumn(1, 3, sec.kineticEnergy);
                am->FillNtupleIColumn(1, 4, hit->pdg);
                am->FillNtupleDColumn(1, 5, fPrimaryEkin / GeV);
                am->AddNtupleRow(1);
            }
        }
    }

    // ── Fill histograms ───────────────────────────────────────────────────────
    auto* am = G4AnalysisManager::Instance();
    am->FillH1(0, edepGAGG / MeV);
    am->FillH1(1, edepBi   / MeV);
    if ((edepGAGG + edepBi) > 0.)
        am->FillH1(2, edepBi / (edepGAGG + edepBi));
    am->FillH1(3, fPrimaryEkin / GeV);
    am->FillH1(4, (G4double)nSecBi);
    am->FillH1(5, (G4double)nSecGAGG);

    // ── Fill ntuple 0 (per-event summary) ────────────────────────────────────
    am->FillNtupleIColumn(0, 0,  evtID);
    am->FillNtupleDColumn(0, 1,  fPrimaryEkin / GeV);
    am->FillNtupleDColumn(0, 2,  edepGAGG / MeV);
    am->FillNtupleDColumn(0, 3,  edepBi   / MeV);
    am->FillNtupleIColumn(0, 4,  crystalHC ? crystalHC->entries() : 0);
    am->FillNtupleIColumn(0, 5,  biHC      ? biHC->entries()      : 0);
    am->FillNtupleIColumn(0, 6,  (G4int)secTypesInBi.size());
    am->FillNtupleDColumn(0, 7,  maxSecKE_Bi);
    am->FillNtupleIColumn(0, 8,  hadronicInBi   ? 1 : 0);
    am->FillNtupleIColumn(0, 9,  hadronicInGAGG ? 1 : 0);
    am->FillNtupleDColumn(0, 10, edepGAGG_heavyIon / MeV);
    am->FillNtupleDColumn(0, 11, edepGAGG_nP       / MeV);
    am->FillNtupleDColumn(0, 12, edepGAGG_em        / MeV);
    am->AddNtupleRow(0);
}
