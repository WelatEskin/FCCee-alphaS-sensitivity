// plot_one_minus_thrust.C
//
// ROOT macro for plotting normalized 1 - thrust distributions
// for different alphaS template samples.
//
// Output:
//   figures/one_minus_thrust_alphaS_scan.pdf
//   figures/one_minus_thrust_alphaS_scan.png
//
// Usage:
//   root -l -b -q macros/plot_one_minus_thrust.C

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"

#include <iostream>
#include <vector>
#include <string>

#include "FCCStyle.C"

void plot_one_minus_thrust()
{
    SetFCCStyle();

    gSystem->mkdir("figures", true);

    // -------------------------------------------------------------------------
    // Input template files
    // -------------------------------------------------------------------------

    std::vector<std::string> fileNames = {
        "data/template_alphaS_0110.root",
        "data/template_alphaS_0114.root",
        "data/template_alphaS_0118.root",
        "data/template_alphaS_0122.root",
        "data/template_alphaS_0126.root",
        "data/template_alphaS_0130.root"
    };

    std::vector<std::string> labels = {
        "#alpha_{s} = 0.110",
        "#alpha_{s} = 0.114",
        "#alpha_{s} = 0.118",
        "#alpha_{s} = 0.122",
        "#alpha_{s} = 0.126",
        "#alpha_{s} = 0.130"
    };

    // Color choice: ROOT standard colors, chosen to be visually distinct.
    std::vector<int> colors = {
        kBlue + 1,
        kAzure + 7,
        kBlack,
        kOrange + 7,
        kRed + 1,
        kMagenta + 2
    };

    std::vector<int> lineStyles = {
        1, 1, 1, 1, 1, 1
    };

    // -------------------------------------------------------------------------
    // Histogram configuration
    //
    // The 1 - T range is restricted to the phenomenologically useful region.
    // Very large values are rare at the Z pole and mostly populate the tail.
    // -------------------------------------------------------------------------

    const int nBins = 35;
    const double xMin = 0.0;
    const double xMax = 0.28;

    std::vector<TH1D*> hists;

    for (size_t i = 0; i < fileNames.size(); ++i) {

        TFile *file = TFile::Open(fileNames[i].c_str(), "READ");

        if (!file || file->IsZombie()) {
            std::cerr << "Error: could not open file: "
                      << fileNames[i] << std::endl;
            return;
        }

        TTree *tree = (TTree*)file->Get("Events");

        if (!tree) {
            std::cerr << "Error: could not find TTree 'Events' in file: "
                      << fileNames[i] << std::endl;
            return;
        }

        std::string histName = "h_oneMinusThrust_" + std::to_string(i);

        TH1D *hist = new TH1D(
            histName.c_str(),
            "",
            nBins,
            xMin,
            xMax
        );

        hist->Sumw2();

        // Fill only physically valid values in the plotted range.
        tree->Draw(
            ("oneMinusThrust >> " + histName).c_str(),
            "oneMinusThrust >= 0.0 && oneMinusThrust < 0.35",
            "goff"
        );

        NormalizeToUnitArea(hist);
        SetHistStyle(hist, colors[i], lineStyles[i]);

        hist->GetXaxis()->SetTitle("1 - T");
        hist->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{d(1 - T)}");

        hists.push_back(hist);
    }

    // Reference histogram for ratio panel: alphaS = 0.118
    TH1D *hRef = hists[2];

    // -------------------------------------------------------------------------
    // Canvas with ratio panel
    // -------------------------------------------------------------------------

    TCanvas *canvas = new TCanvas(
        "canvas",
        "one minus thrust alphaS scan",
        900,
        850
    );

    TPad *padTop = new TPad("padTop", "padTop", 0.0, 0.30, 1.0, 1.0);
    TPad *padBottom = new TPad("padBottom", "padBottom", 0.0, 0.0, 1.0, 0.30);

    padTop->SetBottomMargin(0.03);
    padTop->SetTopMargin(0.07);
    padTop->SetLeftMargin(0.18);
    padTop->SetRightMargin(0.05);

    padBottom->SetTopMargin(0.03);
    padBottom->SetBottomMargin(0.36);
    padBottom->SetLeftMargin(0.18);
    padBottom->SetRightMargin(0.05);

    padTop->Draw();
    padBottom->Draw();

    // -------------------------------------------------------------------------
    // Top panel: normalized distributions
    // -------------------------------------------------------------------------

    padTop->cd();

    double maxY = 0.0;

    for (auto *hist : hists) {
        if (hist->GetMaximum() > maxY) {
            maxY = hist->GetMaximum();
        }
    }

    hists[0]->SetMaximum(1.45 * maxY);
    hists[0]->SetMinimum(0.0);

    hists[0]->GetXaxis()->SetLabelSize(0.0);
    hists[0]->GetXaxis()->SetTitleSize(0.0);

    hists[0]->GetYaxis()->SetTitleSize(0.043);
    hists[0]->GetYaxis()->SetLabelSize(0.038);
    hists[0]->GetYaxis()->SetTitleOffset(1.75);

    hists[0]->Draw("hist");

    for (size_t i = 1; i < hists.size(); ++i) {
        hists[i]->Draw("hist same");
    }

    DrawFCCLabel(0.20, 0.88);

    TLegend *legend = new TLegend(0.66, 0.56, 0.92, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.032);

    for (size_t i = 0; i < hists.size(); ++i) {
        legend->AddEntry(hists[i], labels[i].c_str(), "l");
    }

    legend->Draw();

    // -------------------------------------------------------------------------
    // Bottom panel: ratios to alphaS = 0.118
    // -------------------------------------------------------------------------

    padBottom->cd();

    std::vector<TH1D*> ratios;

    for (size_t i = 0; i < hists.size(); ++i) {
        std::string ratioName = "ratio_" + std::to_string(i);

        TH1D *ratio = MakeRatio(
            hists[i],
            hRef,
            ratioName.c_str()
        );

        SetHistStyle(ratio, colors[i], lineStyles[i]);
        ratio->GetXaxis()->SetTitle("1 - T");
        ratio->GetYaxis()->SetTitle("Ratio");

        ratio->GetXaxis()->SetTitleSize(0.105);
        ratio->GetXaxis()->SetLabelSize(0.085);
        ratio->GetXaxis()->SetTitleOffset(1.05);

        ratio->GetYaxis()->SetTitleSize(0.085);
        ratio->GetYaxis()->SetLabelSize(0.075);
        ratio->GetYaxis()->SetTitleOffset(0.75);
        ratio->GetYaxis()->SetNdivisions(505);

        ratio->SetMinimum(0.70);
        ratio->SetMaximum(1.30);

        ratios.push_back(ratio);
    }

    ratios[0]->Draw("hist");

    for (size_t i = 1; i < ratios.size(); ++i) {
        ratios[i]->Draw("hist same");
    }

    DrawRatioUnityLine(xMin, xMax);

    // -------------------------------------------------------------------------
    // Save output
    // -------------------------------------------------------------------------

    canvas->cd();

    canvas->SaveAs("figures/one_minus_thrust_alphaS_scan.pdf");
    canvas->SaveAs("figures/one_minus_thrust_alphaS_scan.png");

    std::cout << "Output written:" << std::endl;
    std::cout << "  figures/one_minus_thrust_alphaS_scan.pdf" << std::endl;
    std::cout << "  figures/one_minus_thrust_alphaS_scan.png" << std::endl;
}