/*===========================================================
  RunAction.cc
  Books all analysis objects (histograms + ntuple) and
  opens/closes the output ROOT file.
===========================================================*/

#include "RunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Run.hh"

RunAction::RunAction() : G4UserRunAction()
{
    auto* am = G4AnalysisManager::Instance();
    am->SetVerboseLevel(1);
    am->SetNtupleMerging(true);   // merge worker threads in MT mode

    // ── Histograms ────────────────────────────────────────────────────────────
    // H0: Total Edep in GAGG per event
    am->CreateH1("EdepGAGG",  "Edep in GAGG per event",
                 200, 0., 5000., "MeV");

    // H1: Total Edep in Bi per event
    am->CreateH1("EdepBi",    "Edep in Bi-coating per event",
                 200, 0., 500., "MeV");

    // H2: Fraction of Edep in Bi
    am->CreateH1("EdepFracBi", "Edep fraction in Bi / total",
                 100, 0., 1.);

    // H3: Primary neutron energy spectrum
    am->CreateH1("PrimaryE",  "Primary neutron energy",
                 100, 0., 110., "GeV");

    // H4: Number of secondaries per event in Bi
    am->CreateH1("NSecBi",    "# secondaries per event (Bi)",
                 100, 0., 200.);

    // H5: Number of secondaries per event in GAGG
    am->CreateH1("NSecGAGG",  "# secondaries per event (GAGG)",
                 100, 0., 200.);

    // H6: KE spectrum of secondaries produced in Bi [MeV]
    am->CreateH1("SecKEBi",   "Secondary KE in Bi",
                 200, 0., 10000., "MeV");

    // H7: KE spectrum of secondaries produced in GAGG [MeV]
    am->CreateH1("SecKEGAGG", "Secondary KE in GAGG",
                 200, 0., 10000., "MeV");

    // H8: PDG code of secondaries in Bi (fragmentation fingerprint)
    // Range covers standard particles and light ions
    am->CreateH1("SecPDGBi",  "Secondary PDG (Bi)",
                 400, -200., 1000200400.);

    // H9: PDG code of secondaries in GAGG
    am->CreateH1("SecPDGGAGG","Secondary PDG (GAGG)",
                 400, -200., 1000200400.);

    // ── Ntuple ────────────────────────────────────────────────────────────────
    // ── Ntuple 0: per-event summary ───────────────────────────────────────────
    am->CreateNtuple("NeutronStudy", "Per-event summary");
    am->CreateNtupleIColumn("eventID");          // col 0
    am->CreateNtupleDColumn("primaryE_GeV");     // col 1
    am->CreateNtupleDColumn("edepGAGG_MeV");     // col 2
    am->CreateNtupleDColumn("edepBi_MeV");       // col 3
    am->CreateNtupleIColumn("nHitsGAGG");        // col 4
    am->CreateNtupleIColumn("nHitsBi");          // col 5
    am->CreateNtupleIColumn("nSecTypesBi");      // col 6
    am->CreateNtupleDColumn("maxSecKE_Bi_MeV");  // col 7
    am->CreateNtupleIColumn("hadronicInBi");     // col 8
    am->CreateNtupleIColumn("hadronicInGAGG");   // col 9
    am->CreateNtupleDColumn("edepHeavyIon_MeV"); // col 10
    am->CreateNtupleDColumn("edepNP_MeV");       // col 11
    am->CreateNtupleDColumn("edepEM_MeV");       // col 12
    am->FinishNtuple(0);

    // ── Ntuple 1: per-secondary record (one row per secondary particle) ───────
    // This gives us proper PDG spectra for Bi fragments without histogram
    // range/binning problems. Ion PDG codes (~10^9) are stored as doubles.
    am->CreateNtuple("Secondaries", "Per-secondary particle record");
    am->CreateNtupleIColumn("eventID");       // col 0
    am->CreateNtupleIColumn("volume");        // col 1: 0=Bi, 1=GAGG
    am->CreateNtupleDColumn("pdg");           // col 2: as double (ion PDG ~10^9)
    am->CreateNtupleDColumn("ke_MeV");        // col 3: kinetic energy
    am->CreateNtupleIColumn("parentPDG");     // col 4: who produced it
    am->CreateNtupleDColumn("primaryE_GeV");  // col 5: primary neutron energy
    am->FinishNtuple(1);
}

void RunAction::BeginOfRunAction(const G4Run* run)
{
    auto* am = G4AnalysisManager::Instance();
    G4String fname = "neutron_study_run"
                   + std::to_string(run->GetRunID())
                   + ".root";
    am->OpenFile(fname);
    G4cout << "[RunAction] Output file: " << fname << G4endl;
}

void RunAction::EndOfRunAction(const G4Run* run)
{
    auto* am = G4AnalysisManager::Instance();
    am->Write();
    am->CloseFile();
    G4cout << "[RunAction] Run " << run->GetRunID()
           << "  (" << run->GetNumberOfEvent() << " events) — file closed.\n";
}
