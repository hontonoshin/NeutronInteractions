#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TColor.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TMath.h"
#include "TGaxis.h"

#include <cstdio>
#include <algorithm>

// -----------------------------------------------------------------------------
// Global style
// -----------------------------------------------------------------------------
void ApplyPubStyle()
{
    TStyle* s = new TStyle("pub_individual", "pub_individual");

    s->SetCanvasColor(0);
    s->SetFrameFillColor(0);
    s->SetFrameBorderMode(0);
    s->SetCanvasBorderMode(0);
    s->SetPadBorderMode(0);
    s->SetPadColor(0);
    s->SetStatColor(0);

    s->SetTextFont(42);
    s->SetLabelFont(42, "xyz");
    s->SetTitleFont(42, "xyz");
    s->SetLegendFont(42);

    s->SetOptStat(0);
    s->SetOptTitle(0);

    s->SetLabelSize(0.040, "xyz");
    s->SetTitleSize(0.048, "xyz");
    s->SetTitleOffset(1.12, "x");
    s->SetTitleOffset(1.35, "y");

    s->SetPadTickX(1);
    s->SetPadTickY(1);
    s->SetHistLineWidth(2);
    s->SetLegendBorderSize(0);
    s->SetLegendFillStyle(0);
    s->SetGridStyle(3);
    s->SetGridColor(kGray + 1);
    s->SetGridWidth(1);

    s->SetNumberContours(255);
    TGaxis::SetMaxDigits(4);

    gROOT->SetStyle("pub_individual");
    gROOT->ForceStyle();
}

// -----------------------------------------------------------------------------
// Small helpers
// -----------------------------------------------------------------------------
Bool_t HasBranch(TTree* t, const char* bname)
{
    return t && t->GetBranch(bname);
}

TString Sel(TTree* t, const char* condition)
{
    TString c = condition && TString(condition).Length() > 0 ? TString(condition) : TString("1");
    if (HasBranch(t, "weight")) return TString::Format("weight*(%s)", c.Data());
    return c;
}

void LogBins(Int_t n, Double_t lo, Double_t hi, Double_t* edge)
{
    const Double_t llo = TMath::Log10(lo);
    const Double_t lhi = TMath::Log10(hi);
    for (Int_t i = 0; i <= n; ++i)
        edge[i] = TMath::Power(10.0, llo + (lhi - llo) * i / n);
}

void PreparePad(Bool_t logx=false, Bool_t logy=false, Bool_t gridx=true, Bool_t gridy=true,
                Double_t left=0.16, Double_t right=0.075,
                Double_t bottom=0.155, Double_t top=0.055)
{
    gPad->SetLeftMargin(left);
    gPad->SetRightMargin(right);
    gPad->SetBottomMargin(bottom);
    gPad->SetTopMargin(top);
    gPad->SetTicks(1, 1);
    gPad->SetGrid(gridx, gridy);
    gPad->SetLogx(logx);
    gPad->SetLogy(logy);
}

void Format1D(TH1* h, const char* xtitle, const char* ytitle, Int_t color,
              Int_t fillColor=-1, Double_t yTitleOffset=1.32)
{
    h->SetLineColor(color);
    h->SetLineWidth(2);
    if (fillColor >= 0) {
        h->SetFillColor(fillColor);
        h->SetFillStyle(1001);
    } else {
        h->SetFillStyle(0);
    }

    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle(ytitle);
    h->GetXaxis()->SetTitleOffset(1.08);
    h->GetYaxis()->SetTitleOffset(yTitleOffset);
    h->GetXaxis()->SetLabelSize(0.040);
    h->GetYaxis()->SetLabelSize(0.040);
    h->GetXaxis()->SetTitleSize(0.048);
    h->GetYaxis()->SetTitleSize(0.048);
    h->GetXaxis()->CenterTitle(false);
    h->GetYaxis()->CenterTitle(false);
    h->GetXaxis()->SetNdivisions(510);
    h->GetYaxis()->SetNdivisions(510);
    h->GetXaxis()->SetMaxDigits(4);
    h->GetYaxis()->SetMaxDigits(4);
}

void SetLogYRange(TH1* h, Double_t minPositive=0.5, Double_t topFactor=6.0)
{
    if (!h) return;
    Double_t ymax = h->GetMaximum();
    if (ymax <= 0.0) ymax = 1.0;
    h->SetMinimum(minPositive);
    h->SetMaximum(ymax * topFactor);
}

void SaveCanvas(TCanvas* c, const char* outdir, const char* basename)
{

    c->Modified();
    c->Update();
    gSystem->ProcessEvents();

    const TString pdf = TString::Format("%s/%s.pdf", outdir, basename);
    c->Print(pdf, "pdf");
}

TCanvas* NewCanvas(const char* name, Int_t w=920, Int_t h=680)
{
    return new TCanvas(name, "", w, h);
}

// -----------------------------------------------------------------------------
// Main macro
// -----------------------------------------------------------------------------
void analyse(const char* fname="neutron_study_run0.root",
                                const char* outdir="plots")
{
    gROOT->SetBatch(kTRUE);
    ApplyPubStyle();
    TStopwatch timer;
    timer.Start();

    gSystem->mkdir(outdir, kTRUE);

    TFile* f = TFile::Open(fname, "READ");
    if (!f || f->IsZombie()) {
        printf("ERROR: Cannot open ROOT file: %s\n", fname);
        return;
    }

    TTree* t    = (TTree*)f->Get("NeutronStudy");
    TTree* tsec = (TTree*)f->Get("Secondaries");
    if (!t) {
        printf("ERROR: TTree 'NeutronStudy' is missing.\n");
        return;
    }
    if (!tsec) {
        printf("ERROR: TTree 'Secondaries' is missing. Re-run the simulation with secondary ntuple output enabled.\n");
        return;
    }

    const Long64_t Nevt = t->GetEntries();
    const Long64_t Nsec = tsec->GetEntries();
    printf("Events: %lld   Secondaries: %lld\n", Nevt, Nsec);
    printf("Weight branch: %s\n", HasBranch(t, "weight") ? "yes" : "no");

    const Int_t cBlack = kBlack;
    const Int_t cBlue  = TColor::GetColor("#1f77b4");
    const Int_t cRed   = TColor::GetColor("#d62728");
    const Int_t cGreen = TColor::GetColor("#2ca02c");
    const Int_t cOrange= TColor::GetColor("#ff7f0e");
    const Int_t cOrangeFill = TColor::GetColor("#ffbb78");
    const Int_t cPurple= TColor::GetColor("#9467bd");
    const Int_t cPurpleFill = TColor::GetColor("#c5b0d5");

    // -------------------------------------------------------------------------
    // 01. GAGG energy deposition. Use GeV on x-axis to avoid clipped x10^3.
    // -------------------------------------------------------------------------
    {
        TCanvas* c = NewCanvas("c_01_gagg_edep");
        PreparePad(false, true, true, true, 0.165, 0.085, 0.160, 0.050);
        TH1D* h = new TH1D("h_gagg_edep", "", 100, 0.0, 5.0);
        t->Project("h_gagg_edep", "edepGAGG_MeV/1000.0", Sel(t, "edepGAGG_MeV>0"));
        Format1D(h, "Energy deposition in GAGG, E_{dep}^{GAGG} (GeV)",
                 "Weighted counts / 50 MeV", cBlue, -1, 1.35);
        SetLogYRange(h, 0.5, 5.0);
        h->Draw("HIST");
        SaveCanvas(c, outdir, "01_gagg_edep");
    }

    // -------------------------------------------------------------------------
    // 02. Bi coating energy deposition.
    // -------------------------------------------------------------------------
    {
        TCanvas* c = NewCanvas("c_02_bi_coating_edep");
        PreparePad(false, true, true, true, 0.165, 0.085, 0.160, 0.050);
        TH1D* h = new TH1D("h_bi_edep", "", 100, 0.0, 200.0);
        t->Project("h_bi_edep", "edepBi_MeV", Sel(t, "edepBi_MeV>0"));
        Format1D(h, "Energy deposition in Bi coating, E_{dep}^{Bi} (MeV)",
                 "Weighted counts / 2 MeV", cRed, -1, 1.35);
        SetLogYRange(h, 0.5, 6.0);
        h->Draw("HIST");
        SaveCanvas(c, outdir, "02_bi_coating_edep");
    }

    // -------------------------------------------------------------------------
    // 03. Bi/total energy fraction. Scale x by 10^3 to remove ROOT x10^-3.
    // -------------------------------------------------------------------------
    {
        TCanvas* c = NewCanvas("c_03_bi_total_fraction");
        PreparePad(false, true, true, true, 0.165, 0.085, 0.160, 0.050);
        TH1D* h = new TH1D("h_bi_frac", "", 50, 0.0, 1.0);
        t->Project("h_bi_frac",
                   "1000.0*edepBi_MeV/(edepBi_MeV+edepGAGG_MeV)",
                   Sel(t, "edepBi_MeV>0 && (edepBi_MeV+edepGAGG_MeV)>0"));
        Format1D(h, "10^{3} E_{dep}^{Bi}/E_{dep}^{total}",
                 "Weighted counts / 0.02", cGreen, -1, 1.35);
        h->GetXaxis()->SetNoExponent(kTRUE);
        SetLogYRange(h, 0.5, 2.5);
        h->Draw("HIST");
        SaveCanvas(c, outdir, "03_bi_total_edep_fraction");
    }

    // -------------------------------------------------------------------------
    // 04. Primary neutron energy spectrum.
    // -------------------------------------------------------------------------
    {
        TCanvas* c = NewCanvas("c_04_primary_neutron_spectrum");
        PreparePad(true, true, true, true, 0.165, 0.090, 0.160, 0.050);
        const Int_t nb = 50;
        Double_t bins[nb+1];
        LogBins(nb, 9.0, 115.0, bins);
        TH1D* h = new TH1D("h_primary_E", "", nb, bins);
        t->Project("h_primary_E", "primaryE_GeV", Sel(t, "primaryE_GeV>0"));
        Format1D(h, "Primary neutron energy, E_{n} (GeV)",
                 "Weighted events", cBlue, -1, 1.35);
        SetLogYRange(h, 0.5, 2.5);
        h->Draw("HIST");
        SaveCanvas(c, outdir, "04_primary_neutron_energy_spectrum");
    }

    // -------------------------------------------------------------------------
    // 05. Bi secondary kinetic energy. Use GeV.
    // -------------------------------------------------------------------------
    {
        TCanvas* c = NewCanvas("c_05_bi_secondary_ke");
        PreparePad(false, true, true, true, 0.165, 0.085, 0.160, 0.050);
        TH1D* h = new TH1D("h_bi_sec_ke", "", 100, 0.0, 10.0);
        tsec->Project("h_bi_sec_ke", "ke_MeV/1000.0", "volume==0");
        Format1D(h, "Secondary kinetic energy from Bi, E_{k} (GeV)",
                 "Secondaries / 0.1 GeV", cRed, -1, 1.35);
        SetLogYRange(h, 0.5, 6.0);
        h->Draw("HIST");
        SaveCanvas(c, outdir, "05_bi_secondary_kinetic_energy");
    }

    // -------------------------------------------------------------------------
    // 06. GAGG secondary kinetic energy. Use GeV to avoid clipped x10^3.
    // -------------------------------------------------------------------------
    {
        TCanvas* c = NewCanvas("c_06_gagg_secondary_ke");
        PreparePad(false, true, true, true, 0.165, 0.085, 0.160, 0.050);
        TH1D* h = new TH1D("h_gagg_sec_ke", "", 100, 0.0, 10.0);
        tsec->Project("h_gagg_sec_ke", "ke_MeV/1000.0", "volume==1");
        Format1D(h, "Secondary kinetic energy from GAGG, E_{k} (GeV)",
                 "Secondaries / 0.1 GeV", cBlue, -1, 1.35);
        SetLogYRange(h, 0.5, 4.0);
        h->Draw("HIST");
        SaveCanvas(c, outdir, "06_gagg_secondary_kinetic_energy");
    }

    // -------------------------------------------------------------------------
    // 07. GAGG Edep composition. Use GeV on x-axis.
    // -------------------------------------------------------------------------
    {
        TCanvas* c = NewCanvas("c_07_gagg_edep_composition", 980, 720);
        PreparePad(false, true, true, true, 0.165, 0.085, 0.160, 0.050);
        TH1D* hTot = new TH1D("h_comp_tot", "", 80, 0.0, 4.0);
        TH1D* hNP  = new TH1D("h_comp_np",  "", 80, 0.0, 4.0);
        TH1D* hHI  = new TH1D("h_comp_hi",  "", 80, 0.0, 4.0);

        t->Project("h_comp_tot", "edepGAGG_MeV/1000.0", Sel(t, "edepGAGG_MeV>0"));
        t->Project("h_comp_np",  "edepNP_MeV/1000.0",   Sel(t, "edepNP_MeV>0"));
        t->Project("h_comp_hi",  "edepHeavyIon_MeV/1000.0", Sel(t, "edepHeavyIon_MeV>0"));

        Format1D(hTot, "Energy deposition in GAGG, E_{dep} (GeV)",
                 "Weighted counts / 50 MeV", cBlack, -1, 1.35);
        Format1D(hNP,  "", "", cBlue);
        Format1D(hHI,  "", "", cRed);
        SetLogYRange(hTot, 0.5, 6.0);
        hTot->Draw("HIST");
        hNP->Draw("HIST SAME");
        hHI->Draw("HIST SAME");

        TLegend* leg = new TLegend(0.56, 0.74, 0.92, 0.91);
        leg->SetTextSize(0.036);
        leg->AddEntry(hTot, "Total GAGG E_{dep}", "l");
        leg->AddEntry(hNP,  "n + p recoil component", "l");
        leg->AddEntry(hHI,  "Heavy ions, A > 4", "l");
        leg->Draw();

        SaveCanvas(c, outdir, "07_gagg_edep_composition");
    }

    // -------------------------------------------------------------------------
    // 08. EM energy deposition in GAGG.
    // -------------------------------------------------------------------------
    {
        TCanvas* c = NewCanvas("c_08_gagg_em_edep");
        PreparePad(false, true, true, true, 0.165, 0.085, 0.160, 0.050);
        TH1D* h = new TH1D("h_em_edep", "", 80, 0.0, 1000.0);
        t->Project("h_em_edep", "edepEM_MeV", Sel(t, "edepEM_MeV>0"));
        Format1D(h, "Electromagnetic energy deposition, E_{dep}^{EM} (MeV)",
                 "Weighted counts / 12.5 MeV", cGreen, -1, 1.35);
        SetLogYRange(h, 0.5, 5.0);
        h->Draw("HIST");
        SaveCanvas(c, outdir, "08_gagg_em_edep");
    }

    // -------------------------------------------------------------------------
    // 09. Heavy-ion Edep fraction. Log-y avoids clipped x10^6 y multiplier.
    // -------------------------------------------------------------------------
    {
        TCanvas* c = NewCanvas("c_09_heavy_ion_fraction");
        PreparePad(false, true, true, true, 0.165, 0.085, 0.160, 0.050);
        TH1D* h = new TH1D("h_hi_frac", "", 50, 0.0, 0.5);
        t->Project("h_hi_frac",
                   "edepHeavyIon_MeV/edepGAGG_MeV",
                   Sel(t, "edepGAGG_MeV>1 && edepHeavyIon_MeV>=0"));
        Format1D(h, "Heavy-ion energy fraction, E_{dep}^{HI}/E_{dep}^{GAGG}",
                 "Weighted events / 0.01", cOrange, cOrangeFill, 1.35);
        SetLogYRange(h, 0.5, 4.0);
        h->Draw("HIST");
        SaveCanvas(c, outdir, "09_gagg_heavy_ion_edep_fraction");
    }

    // -------------------------------------------------------------------------
    // 10. Standard secondaries from Bi by PDG code.
    // -------------------------------------------------------------------------
    {
        TCanvas* c = NewCanvas("c_10_bi_standard_secondaries");
        PreparePad(false, true, true, true, 0.165, 0.085, 0.160, 0.050);
        TH1D* h = new TH1D("h_std_pdg", "", 110, -250.0, 2500.0);
        tsec->Project("h_std_pdg", "pdg", "volume==0 && pdg>-250 && pdg<2500");
        Format1D(h, "PDG code", "Secondaries per code", cBlue, -1, 1.35);
        SetLogYRange(h, 0.5, 6.0);
        h->Draw("HIST");

        TLatex lab;
        lab.SetTextFont(42);
        lab.SetTextSize(0.032);
        lab.SetTextAlign(21);
        struct PdgLabel { Double_t x; const char* txt; } labels[] = {
            {22, "#gamma"}, {2112, "n"}, {2212, "p"}
        };
        for (auto& p : labels) {
            Double_t y = h->GetBinContent(h->FindBin(p.x));
            if (y > 0.5) lab.DrawLatex(p.x, y * 1.7, p.txt);
        }
        SaveCanvas(c, outdir, "10_bi_standard_secondaries_pdg");
    }

    // -------------------------------------------------------------------------
    // 11. Light ion secondaries from Bi.
    // -------------------------------------------------------------------------
    {
        TCanvas* c = NewCanvas("c_11_bi_light_ions");
        PreparePad(false, true, true, true, 0.165, 0.085, 0.160, 0.050);
        TH1D* h = new TH1D("h_light_ions", "", 7, 0.5, 7.5);
        const char* ionLab[] = {"d", "t", "^{3}He", "#alpha", "^{6}Li", "^{7}Li", "other"};
        for (Int_t i = 1; i <= 7; ++i) h->GetXaxis()->SetBinLabel(i, ionLab[i-1]);

        h->SetBinContent(1, tsec->GetEntries("volume==0 && pdg==1000010020"));
        h->SetBinContent(2, tsec->GetEntries("volume==0 && pdg==1000010030"));
        h->SetBinContent(3, tsec->GetEntries("volume==0 && pdg==1000020030"));
        h->SetBinContent(4, tsec->GetEntries("volume==0 && pdg==1000020040"));
        h->SetBinContent(5, tsec->GetEntries("volume==0 && pdg==1000030060"));
        h->SetBinContent(6, tsec->GetEntries("volume==0 && pdg==1000030070"));
        h->SetBinContent(7, tsec->GetEntries(
            "volume==0 && pdg>1000010000 && pdg<1000050000"
            " && pdg!=1000010020 && pdg!=1000010030"
            " && pdg!=1000020030 && pdg!=1000020040"
            " && pdg!=1000030060 && pdg!=1000030070"));

        Format1D(h, "Ion type", "Count", cOrange, cOrangeFill, 1.35);
        h->GetXaxis()->SetLabelSize(0.045);
        SetLogYRange(h, 0.5, 6.0);
        h->Draw("HIST");
        SaveCanvas(c, outdir, "11_bi_light_ion_secondaries");
    }

    // -------------------------------------------------------------------------
    // 12. Spallation fragment Z distribution from Bi.
    // -------------------------------------------------------------------------
    {
        TCanvas* c = NewCanvas("c_12_bi_fragment_Z");
        PreparePad(false, true, true, true, 0.165, 0.085, 0.160, 0.050);
        TH1D* h = new TH1D("h_fragment_Z", "", 14, 69.5, 83.5);
        tsec->Project("h_fragment_Z",
                      "fmod(floor(pdg/10000),1000)",
                      "volume==0 && pdg>1000050000"
                      " && fmod(floor(pdg/10000),1000)>=70"
                      " && fmod(floor(pdg/10000),1000)<=83");
        Format1D(h, "Atomic number Z", "Fragments", cPurple, cPurpleFill, 1.35);
        h->GetXaxis()->SetNdivisions(514);
        SetLogYRange(h, 0.5, 4.0);
        h->Draw("HIST");

        TLatex lab;
        lab.SetTextFont(42);
        lab.SetTextSize(0.030);
        lab.SetTextAlign(21);
        struct ElLabel { Int_t Z; const char* txt; } els[] = {{80,"Hg"},{81,"Tl"},{82,"Pb"},{83,"Bi"}};
        for (auto& e : els) {
            Double_t y = h->GetBinContent(h->FindBin(e.Z));
            if (y > 0.5) lab.DrawLatex((Double_t)e.Z, y * 1.45, e.txt);
        }
        SaveCanvas(c, outdir, "12_bi_spallation_fragment_z");
    }

    // -------------------------------------------------------------------------
    // 13. Spallation fragment A-Z map. Axis labels are bin labels, no overlapping
    // external TLatex labels.
    // -------------------------------------------------------------------------
    {
        TCanvas* c = NewCanvas("c_13_bi_fragment_AZ", 1080, 760);
        PreparePad(false, false, false, false, 0.155, 0.170, 0.140, 0.055);

        TH2D* h = new TH2D("h_fragment_AZ", "",
                           41, 169.5, 210.5,
                           14,  69.5,  83.5);

        tsec->Project("h_fragment_AZ",
                      "fmod(floor(pdg/10000),1000):fmod(floor(pdg/10),1000)",
                      "volume==0 && pdg>1000050000"
                      " && fmod(floor(pdg/10000),1000)>=70"
                      " && fmod(floor(pdg/10000),1000)<=83"
                      " && fmod(floor(pdg/10),1000)>=170"
                      " && fmod(floor(pdg/10),1000)<=210");

        gStyle->SetPalette(kViridis);
        if (h->GetMaximum() > 10.0) gPad->SetLogz();

        h->GetXaxis()->SetTitle("Mass number A");
        h->GetYaxis()->SetTitle("Atomic number Z");
        h->GetZaxis()->SetTitle("Fragment yield");
        h->GetXaxis()->SetTitleSize(0.046);
        h->GetYaxis()->SetTitleSize(0.046);
        h->GetZaxis()->SetTitleSize(0.046);
        h->GetXaxis()->SetLabelSize(0.036);
        h->GetYaxis()->SetLabelSize(0.032);
        h->GetZaxis()->SetLabelSize(0.034);
        h->GetXaxis()->SetTitleOffset(1.05);
        h->GetYaxis()->SetTitleOffset(1.22);
        h->GetZaxis()->SetTitleOffset(1.25);
        h->GetXaxis()->SetNdivisions(508);
        h->GetYaxis()->SetNdivisions(514);
        h->GetZaxis()->SetMaxDigits(3);

        const char* yLabels[] = {
            "70 Yb", "71 Lu", "72 Hf", "73 Ta", "74 W", "75 Re", "76 Os",
            "77 Ir", "78 Pt", "79 Au", "80 Hg", "81 Tl", "82 Pb", "83 Bi"
        };
        for (Int_t i = 1; i <= 14; ++i)
            h->GetYaxis()->SetBinLabel(i, yLabels[i-1]);

        if (h->GetMaximum() < 0.5) {
            h->Draw("AXIS");
            TLatex msg;
            msg.SetNDC();
            msg.SetTextFont(42);
            msg.SetTextSize(0.038);
            msg.SetTextAlign(22);
            msg.DrawLatex(0.50, 0.55, "No heavy fragments in the selected A-Z range");
        } else {
            h->Draw("COLZ");
        }

        SaveCanvas(c, outdir, "13_bi_spallation_fragment_AZ_chart");
    }

    // -------------------------------------------------------------------------
    // Text summary.
    // -------------------------------------------------------------------------
    TH1D* hMeanG  = new TH1D("h_mean_gagg", "", 100, 0.0, 1.0e5);
    TH1D* hMeanBi = new TH1D("h_mean_bi",   "", 100, 0.0, 1.0e4);
    TH1D* hMeanHI = new TH1D("h_mean_hi",   "", 100, 0.0, 1.0e5);
    TH1D* hMeanNP = new TH1D("h_mean_np",   "", 100, 0.0, 1.0e5);
    TH1D* hMeanEM = new TH1D("h_mean_em",   "", 100, 0.0, 1.0e5);

    t->Project("h_mean_gagg", "edepGAGG_MeV", "");
    t->Project("h_mean_bi",   "edepBi_MeV", "");
    t->Project("h_mean_hi",   "edepHeavyIon_MeV", "");
    t->Project("h_mean_np",   "edepNP_MeV", "");
    t->Project("h_mean_em",   "edepEM_MeV", "");

    printf("\n==================================================\n");
    printf("  Events                 : %lld\n", Nevt);
    printf("  Mean Edep GAGG         : %.4f MeV\n", hMeanG->GetMean());
    printf("  Mean Edep Bi coating   : %.6f MeV\n", hMeanBi->GetMean());
    printf("  Mean Edep heavy ions   : %.4f MeV\n", hMeanHI->GetMean());
    printf("  Mean Edep n/p          : %.4f MeV\n", hMeanNP->GetMean());
    printf("  Mean Edep e/gamma      : %.4f MeV\n", hMeanEM->GetMean());
    if (HasBranch(t, "hadronicInBi"))
        printf("  Had. inel. in Bi   %%   : %.4f\n", 100.0 * t->GetEntries("hadronicInBi==1") / (Double_t)Nevt);
    if (HasBranch(t, "hadronicInGAGG"))
        printf("  Had. inel. in GAGG %%   : %.4f\n", 100.0 * t->GetEntries("hadronicInGAGG==1") / (Double_t)Nevt);
    printf("  Secondaries ntuple     : %lld\n", Nsec);
    printf("  Output directory       : %s\n", outdir);
    printf("==================================================\n");

    timer.Stop();
    printf("Total macro time: %.2f s\n", timer.RealTime());
}
