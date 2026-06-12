#pragma once
// DetectorConstruction.hh
// GAGG crystal array with 5-µm Bi-209 coating.
// Optical physics intentionally excluded.

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    // Array layout: 4 × 4 × 2
    static constexpr G4int NX = 4;
    static constexpr G4int NY = 4;
    static constexpr G4int NZ = 2;

private:
    void DefineMaterials();

    G4LogicalVolume* fCrystalLogical = nullptr;
    G4LogicalVolume* fBiLogical      = nullptr;
    G4LogicalVolume* fWorldLogical   = nullptr;
};
