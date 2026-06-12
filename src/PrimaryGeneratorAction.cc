#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include <cmath>

G4ThreadLocal G4double PrimaryGeneratorAction::fgCurrentEnergy = 0.0;
G4ThreadLocal G4double PrimaryGeneratorAction::fgCurrentWeight = 1.0;

// ── True Gordon 2004 component probabilities ─────────────────────────────────
namespace True {
    constexpr double w_th  = 0.06;
    constexpr double w_ep  = 0.20;
    constexpr double w_ev  = 0.55;
    constexpr double w_hi  = 0.19;
}
// ── Biased (importance-sampled) probabilities ─────────────────────────────────
namespace Bias {
    constexpr double w_th  = 0.01;   // 6x reduced
    constexpr double w_ep  = 0.04;   // 5x reduced
    constexpr double w_ev  = 0.15;   // 3.7x reduced
    constexpr double w_hi  = 0.80;   // 4.2x enhanced
}
// ── Energy boundaries (MeV) ──────────────────────────────────────────────────
namespace E {
    constexpr double th_lo  = 1.0e-9;
    constexpr double th_hi  = 5.0e-7;   // ~0.5 eV
    constexpr double ep_hi  = 0.1;      // 100 keV
    constexpr double ev_hi  = 20.0;     // 20 MeV
    constexpr double hi_hi  = 1.0e5;    // 100 GeV
    constexpr double kT     = 0.0253e-3;// thermal kT in MeV
    constexpr double T_nuc  = 1.0;      // nuclear temperature (MeV)
    constexpr double gamma  = 1.7;      // power-law index above 20 MeV
}

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(),
      fParticleGun(new G4ParticleGun(1))
{
    auto* n = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
    fParticleGun->SetParticleDefinition(n);
    fParticleGun->SetParticleEnergy(1.0*MeV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,-1));
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{ delete fParticleGun; }

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
    G4double weight = 1.0;
    const G4double energy = SampleEnergy(weight);

    fgCurrentEnergy = energy;
    fgCurrentWeight = weight;

    fParticleGun->SetParticleEnergy(energy);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));

    // ── Position: uniform over the full calorimeter front face ───────────────
    // Crystal: 2.5×2.5 cm, Bi shell: +5µm per face, gap: 0.1 mm
    // pitchXY = 2*(1.25 cm + 5e-4 cm) + 0.01 cm = 2.51001 cm
    // caloHX  = NX * pitchXY / 2 = 4 * 2.51001 / 2 = 5.02002 cm
    // Same for Y.  Use exact arithmetic to match DetectorConstruction.
    constexpr G4double crystalHX = 2.5*cm / 2.0;
    constexpr G4double biThick   = 5.0*um;
    constexpr G4double gap       = 0.1*mm;
    constexpr G4double pitchXY   = 2.0*(crystalHX + biThick) + gap;
    constexpr G4double caloHX    = DetectorConstruction::NX * pitchXY / 2.0;
    constexpr G4double caloHY    = DetectorConstruction::NY * pitchXY / 2.0;

    const G4double x0 = (2.0*G4UniformRand() - 1.0) * caloHX;
    const G4double y0 = (2.0*G4UniformRand() - 1.0) * caloHY;

    fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, fSourceZ));
    fParticleGun->GeneratePrimaryVertex(event);
}

G4double PrimaryGeneratorAction::SampleEnergy(G4double& weight) const
{
    const double r = G4UniformRand();
    G4double energy = 0.;

    // ── Select component using BIASED weights ────────────────────────────────
    if (r < Bias::w_th) {
        // Thermal: Maxwell-Boltzmann ∝ E·exp(-E/kT)
        G4double E; int tries=0;
        do {
            double u1=std::max(G4UniformRand(),1e-30);
            double u2=std::max(G4UniformRand(),1e-30);
            E = E::kT*(-std::log(u1)-std::log(u2));
            ++tries;
        } while((E<E::th_lo||E>E::th_hi)&&tries<1000);
        energy = std::max(E::th_lo,std::min(E,E::th_hi)) * MeV;
        // Weight = (true component prob × biased total) / (biased component prob × true total)
        // Since total = 1 for both: weight = true_w / biased_w
        weight = True::w_th / Bias::w_th;
    }
    else if (r < Bias::w_th + Bias::w_ep) {
        // Epithermal: 1/E → log-uniform
        energy = E::th_hi * std::pow(E::ep_hi/E::th_hi, G4UniformRand()) * MeV;
        weight = True::w_ep / Bias::w_ep;
    }
    else if (r < Bias::w_th + Bias::w_ep + Bias::w_ev) {
        // Evaporation: Maxwell spectrum ∝ E·exp(-E/T_nuc)
        G4double E; int tries=0;
        do {
            double u1=std::max(G4UniformRand(),1e-30);
            double u2=std::max(G4UniformRand(),1e-30);
            E = E::T_nuc*(-std::log(u1)-std::log(u2));
            ++tries;
        } while((E<E::ep_hi||E>E::ev_hi)&&tries<1000);
        energy = std::max(E::ep_hi,std::min(E,E::ev_hi)) * MeV;
        weight = True::w_ev / Bias::w_ev;
    }
    else {
        // High-energy: power law E^{-γ}, inverse CDF
        const double alpha = 1.0 - E::gamma;   // = -0.7
        const double Ea = std::pow(E::ev_hi, alpha);
        const double Eb = std::pow(E::hi_hi, alpha);
        double Ep = Ea + G4UniformRand()*(Eb-Ea);
        if(Ep<=0.) Ep=Ea;
        energy = std::pow(Ep, 1./alpha) * MeV;
        weight = True::w_hi / Bias::w_hi;
    }
    return energy;
}

G4ThreeVector PrimaryGeneratorAction::SampleHemisphere() const
{
    // cos²θ weighting: cos θ = (rand)^{1/3}
    const double cosTheta = std::cbrt(std::max(G4UniformRand(), 1e-10));
    const double sinTheta = std::sqrt(std::max(0., 1.-cosTheta*cosTheta));
    const double phi = CLHEP::twopi * G4UniformRand();
    return G4ThreeVector(sinTheta*std::cos(phi),
                         sinTheta*std::sin(phi),
                        -cosTheta);
}
