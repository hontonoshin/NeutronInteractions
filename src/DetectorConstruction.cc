/*===========================================================
  DetectorConstruction.cc
  
  Geometry:
    4 × 4 × 2 array of GAGG crystals (2.5 × 2.5 × 3.0 cm each)
    Each crystal wrapped in a 5-µm Bi-209 shell on all faces.
    
  No optical physics: no MPT, no scintillation, no Cherenkov.
  Focus: hadronic Edep + secondary-particle tracking.
  
  Two sensitive volumes:
    "Crystal"   — GAGG bulk    → CrystalSD
    "BiCoating" — Bi-209 shell → BiCoatingSD
===========================================================*/

#include "DetectorConstruction.hh"
#include "CrystalSD.hh"
#include "BiCoatingSD.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"

DetectorConstruction::DetectorConstruction()
    : G4VUserDetectorConstruction() {}

// ─────────────────────────────────────────────────────────
void DetectorConstruction::DefineMaterials()
{
    G4NistManager* nist = G4NistManager::Instance();

    // ── GAGG: Gd₃Al₂Ga₃O₁₂  (ρ = 6.63 g/cm³) ──────────────────────────────
    G4Element* Gd = nist->FindOrBuildElement("Gd");
    G4Element* Al = nist->FindOrBuildElement("Al");
    G4Element* Ga = nist->FindOrBuildElement("Ga");
    G4Element* O  = nist->FindOrBuildElement("O");

    auto* GAGG = new G4Material("GAGG", 6.63*g/cm3, 4, kStateSolid);
    GAGG->AddElement(Gd,  3);
    GAGG->AddElement(Al,  2);
    GAGG->AddElement(Ga,  3);
    GAGG->AddElement(O,  12);
    // No optical MPT — optical physics not registered.

    // ── Bi-209 (pure isotope, ρ = 9.747 g/cm³) ──────────────────────────────
    // Defined as a pure isotope so NeutronHP uses the correct Bi-209
    // cross-section data (ENDF/B-VII Bi-209 file).
    auto* bi209  = new G4Isotope("Bi209", 83, 209, 208.980399*g/mole);
    auto* Bi_el  = new G4Element("Bismuth_pure", "Bi", 1);
    Bi_el->AddIsotope(bi209, 1.0);

    auto* BiMat = new G4Material("Bismuth", 9.747*g/cm3, 1, kStateSolid);
    BiMat->AddElement(Bi_el, 1);

    // Air (world fill)
    nist->FindOrBuildMaterial("G4_AIR");
}

// ─────────────────────────────────────────────────────────
G4VPhysicalVolume* DetectorConstruction::Construct()
{
    DefineMaterials();

    auto* air   = G4Material::GetMaterial("G4_AIR");
    auto* GAGG  = G4Material::GetMaterial("GAGG");
    auto* BiMat = G4Material::GetMaterial("Bismuth");

    // ── Dimensions ───────────────────────────────────────────────────────────
    const G4double crystalX = 2.5*cm / 2.0;   // half-lengths
    const G4double crystalY = 2.5*cm / 2.0;
    const G4double crystalZ = 3.0*cm / 2.0;
    const G4double biThick  = 5.0*um;           // 5 µm per face

    // Bi outer half-lengths (shell encloses crystal on all 6 faces)
    const G4double biX = crystalX + biThick;
    const G4double biY = crystalY + biThick;
    const G4double biZ = crystalZ + biThick;

    // Small air gap between adjacent Bi shells (avoids overlapping volumes)
    const G4double gap    = 0.1*mm;
    const G4double pitchXY = 2.0*biX + gap;
    const G4double pitchZ  = 2.0*biZ + gap;

    // ── World ─────────────────────────────────────────────────────────────────
    const G4double caloHX = NX * pitchXY / 2.0;
    const G4double caloHY = NY * pitchXY / 2.0;
    const G4double caloHZ = NZ * pitchZ  / 2.0;

    auto* worldBox = new G4Box("World",
                               3.0*caloHX + 10.*cm,
                               3.0*caloHY + 10.*cm,
                               3.0*caloHZ + 30.*cm);
    fWorldLogical = new G4LogicalVolume(worldBox, air, "World");
    auto* worldPhys = new G4PVPlacement(nullptr, {}, fWorldLogical,
                                        "World", nullptr, false, 0, true);

    // ── Calorimeter envelope ──────────────────────────────────────────────────
    auto* caloBox = new G4Box("Calorimeter", caloHX, caloHY, caloHZ);
    auto* caloLog = new G4LogicalVolume(caloBox, air, "Calorimeter");
    new G4PVPlacement(nullptr, {}, caloLog, "Calorimeter",
                      fWorldLogical, false, 0, true);

    // ── Bi shell (mother of crystal) ─────────────────────────────────────────
    auto* biBox = new G4Box("BiCoating", biX, biY, biZ);
    fBiLogical  = new G4LogicalVolume(biBox, BiMat, "BiCoating");

    // ── GAGG crystal (daughter inside Bi shell) ───────────────────────────────
    auto* crystalBox = new G4Box("Crystal", crystalX, crystalY, crystalZ);
    fCrystalLogical  = new G4LogicalVolume(crystalBox, GAGG, "Crystal");

    // Crystal centred inside Bi shell
    new G4PVPlacement(nullptr, {}, fCrystalLogical, "Crystal",
                      fBiLogical, false, 0, true);

    // ── Place NX × NY × NZ Bi-wrapped crystals ───────────────────────────────
    G4int copyNum = 0;
    for (G4int iz = 0; iz < NZ; ++iz) {
        for (G4int iy = 0; iy < NY; ++iy) {
            for (G4int ix = 0; ix < NX; ++ix) {
                G4double xPos = (ix - (NX-1)/2.0) * pitchXY;
                G4double yPos = (iy - (NY-1)/2.0) * pitchXY;
                G4double zPos = (iz - (NZ-1)/2.0) * pitchZ;
                new G4PVPlacement(nullptr, {xPos, yPos, zPos},
                                  fBiLogical, "BiCoating",
                                  caloLog, false, copyNum++, true);
            }
        }
    }

    // ── Visualisation ─────────────────────────────────────────────────────────
    fWorldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
    caloLog->SetVisAttributes(G4VisAttributes::GetInvisible());

    auto* biVis = new G4VisAttributes(G4Colour(0.7, 0.7, 0.9, 0.3));
    biVis->SetForceSolid(true);
    fBiLogical->SetVisAttributes(biVis);

    auto* crystalVis = new G4VisAttributes(G4Colour(0.0, 0.8, 0.4, 0.6));
    crystalVis->SetForceWireframe(true);
    fCrystalLogical->SetVisAttributes(crystalVis);

    return worldPhys;
}

// ─────────────────────────────────────────────────────────
void DetectorConstruction::ConstructSDandField()
{
    auto* sdMgr = G4SDManager::GetSDMpointer();

    auto* crystalSD = new CrystalSD("CrystalSD", "CrystalHitsCollection");
    sdMgr->AddNewDetector(crystalSD);
    SetSensitiveDetector("Crystal", crystalSD);

    auto* biSD = new BiCoatingSD("BiCoatingSD", "BiHitsCollection");
    sdMgr->AddNewDetector(biSD);
    SetSensitiveDetector("BiCoating", biSD);
}
