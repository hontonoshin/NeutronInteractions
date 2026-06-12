/*===========================================================
  PhysicsList.cc
  
  Physics choices for neutron + Bi/GAGG study (10–100 GeV):
  
  EM:
    G4EmStandardPhysics — ionisation, bremsstrahlung, photoelectric
    for secondaries (protons, alphas, electrons from fragmentation).
    No optical physics registered (no scintillation/Cherenkov).
  
  Hadronic elastic:
    G4HadronElasticPhysicsHP — uses NeutronHP data for n < 20 MeV.
    At 10–100 GeV the standard parametric elastic takes over.
  
  Hadronic inelastic:
    G4HadronPhysicsINCLXX(true, true, false)
      arg1 = quasiElastic ON (uses INCL++ quasi-elastic above 20 MeV)
      arg2 = enable NeutronHP for inelastic n < 20 MeV
      arg3 = Bertini OFF (INCL++ handles all energies here)
    INCL++ valid range: ~1 MeV – ~3 GeV nucleon.
    Above ~3 GeV/nucleon INCL++ smoothly hands off to FTFP (internal
    to G4HadronPhysicsINCLXX).  At 10–100 GeV neutrons FTFP+BERT
    is used for projectile; INCL++ still handles target remnant 
    de-excitation via ABLA07 (fission, evaporation).
    
  Ion physics:
    G4IonINCLXXPhysics — transports all heavy fragments produced
    (Tl, Pb recoils, alphas, Li, Be, ...) correctly.
  
  Decay:
    G4DecayPhysics + G4RadioactiveDecayPhysics — handles gamma
    emission from excited nuclear states and β decays of
    fragmentation products.
  
  Neutron kill cut:
    G4NeutronTrackingCut at 1 ms — stops thermal neutrons that
    would otherwise diffuse for enormous times.
===========================================================*/

#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

// EM
#include "G4EmStandardPhysics.hh"

// Decays
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

// Hadronic elastic (with NeutronHP)
#include "G4HadronElasticPhysicsHP.hh"

// Hadronic inelastic — INCL++ backbone + ABLA07 + FTFP at high E
#include "G4HadronPhysicsINCLXX.hh"
#include "G4INCLXXInterfaceStore.hh"

// Ion transport
#include "G4IonINCLXXPhysics.hh"

// Neutron kill cut
#include "G4NeutronTrackingCut.hh"

// Nuclear de-excitation parameters (ABLA07 fission channel)
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"

// Particle constructors
#include "G4BaryonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4ProductionCutsTable.hh"
#include "G4DeexPrecoParameters.hh"

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
    SetVerboseLevel(1);

    // ── Electromagnetic (standard, handles secondaries from fragmentation) ────
    // No G4EmExtraPhysics, no G4OpticalPhysics — intentionally excluded.
    RegisterPhysics(new G4EmStandardPhysics(1));

    // ── Decays ────────────────────────────────────────────────────────────────
    RegisterPhysics(new G4DecayPhysics(1));
    RegisterPhysics(new G4RadioactiveDecayPhysics(1));

    // ── Hadronic elastic — NeutronHP data for n < 20 MeV ─────────────────────
    RegisterPhysics(new G4HadronElasticPhysicsHP(1));

    // ── Hadronic inelastic — INCL++/ABLA07 + FTFP at high energy ────────────
    //   Constructor signature: (name, quasiElastic, neutronHP, bert)
    //   quasiElastic = true  → INCL++ quasi-elastic model active
    //   neutronHP    = true  → HP data for n < 20 MeV inelastic
    //   bert         = false → no Bertini, INCL++ handles spallation range
    RegisterPhysics(new G4HadronPhysicsINCLXX("inclxx", true, true, false));

    // ── Ion physics — INCL++ for light fragments, parametric for heavy ────────
    RegisterPhysics(new G4IonINCLXXPhysics(1));

    // ── Kill thermal neutrons after 1 ms (prevents infinite diffusion) ────────
    auto* ntc = new G4NeutronTrackingCut(1);
    ntc->SetTimeLimit(1.0*CLHEP::ms);
    RegisterPhysics(ntc);
}

//PhysicsList::~PhysicsList() = default;

void PhysicsList::ConstructParticle()
{
    G4BaryonConstructor().ConstructParticle();
    G4BosonConstructor().ConstructParticle();
    G4LeptonConstructor().ConstructParticle();
    G4MesonConstructor().ConstructParticle();
    G4IonConstructor().ConstructParticle();
    G4ShortLivedConstructor().ConstructParticle();
    G4VModularPhysicsList::ConstructParticle();
}

void PhysicsList::ConstructProcess()
{
    G4VModularPhysicsList::ConstructProcess();

    // ── ABLA07 fission channel (Bi(n,f) enabled) ─────────────────────────────
    // fDeexChannelType = 4 → fAbla (ABLA07)
    // This sets the de-excitation model used by INCL++ after spallation.
    auto* deexParams = G4NuclearLevelData::GetInstance()->GetParameters();
    deexParams->SetDeexChannelsType(G4DeexChannelType(4));

    // ── INCL++ cluster emission up to A = 8 ───────────────────────────────────
    // Allows emission of He-3, He-4, Li-6/7, Be-8, ... from Bi remnant.
    G4INCLXXInterfaceStore::GetInstance()->SetMaxClusterMass(8);

    G4cout << "\n[PhysicsList] Configuration summary:"
           << "\n  EM:       G4EmStandardPhysics (NO optical)"
           << "\n  Elastic:  G4HadronElasticPhysicsHP (NeutronHP < 20 MeV)"
           << "\n  Inelast:  G4HadronPhysicsINCLXX (INCL++/ABLA07, FTFP>3GeV)"
           << "\n  Ions:     G4IonINCLXXPhysics"
           << "\n  Deex:     ABLA07 fission channel"
           << "\n  Clusters: up to A=8 from INCL++"
           << "\n  Neutron kill: 1 ms time cut\n"
           << G4endl;
}

void PhysicsList::SetCuts()
{
    // Production cuts — 1 mm for photons/e±, 0.1 mm for protons.
    // In a study focused on Bi fragmentation and GAGG Edep these are
    // reasonable: they avoid tracking micro-delta-rays while keeping
    // all physically meaningful secondary tracks.
    SetCutValue(1.0*mm, "gamma");
    SetCutValue(1.0*mm, "e-");
    SetCutValue(1.0*mm, "e+");
    SetCutValue(0.1*mm, "proton");

    // Energy range for production cuts table
    G4ProductionCutsTable::GetProductionCutsTable()
        ->SetEnergyRange(100.*eV, 100.*TeV);
}
