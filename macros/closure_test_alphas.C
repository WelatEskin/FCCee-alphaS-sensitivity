// closure_test_alphas.C
//
// Template-based pseudo-data closure test for alphaS extraction.
//
// Observable:
//   oneMinusThrust = 1 - T
//
// Pseudo-data:
//   alphaS_true = 0.118, independent random seed
//
// Templates:
//   alphaS = 0.110, 0.114, 0.118, 0.122, 0.126, 0.130
//
// Output:
//   figures/closure_test_alphas.pdf
//   figures/closure_test_alphas.png
//
// Usage:
//   root -l -b -q macros/closure_test_alphas.C

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>

#include "FCCStyle.C"

// -----------------------------------------------------------------------------
// Helper function: create normalized histogram from a ROOT TTree
// -----------------------------------------------------------------------------

TH1D* MakeNormalizedHist(
    const std::string& fileName,
    const std::string& histName,
    int nBins,
    double xMin,
    double xMax,
    const std::string& selection
)
{
    TFile *file = TFile::Open(fileName.c_str(), "READ");

    if (!file || file->IsZombie()) {
        std::cerr << "Error: could not open file: "
                  << fileName << std::endl;
        return nullptr;
    }

    TTree *tree = (TTree*)file->Get("Events");

    if (!tree) {
        std::cerr << "Error: could not find TTree 'Events' in file: "
                  << fileName << std::endl;
        return nullptr;
    }

    TH1D *hist = new TH1D(histName.c_str(), "", nBins, xMin, xMax);
    hist->Sumw2();

    tree->Draw(
        ("oneMinusThrust >> " + histName).c_str(),
        selection.c_str(),
        "goff"
    );

    NormalizeToUnitArea(hist);

    return hist;
}

// -----------------------------------------------------------------------------
// Binned chi-square using normalized histograms
// -----------------------------------------------------------------------------

double ComputeChi2(
    TH1D *data,
    TH1D *templ,
    double fitMin,
    double fitMax,
    int &nUsedBins
)
{
    if (!data || !templ) return -1.0;

    double chi2 = 0.0;
    nUsedBins = 0;

    for (int i = 1; i <= data->GetNbinsX(); ++i) {

        double x = data->GetBinCenter(i);

        if (x < fitMin || x > fitMax) continue;

        double d = data->GetBinContent(i);
        double t = templ->GetBinContent(i);

        double ed = data->GetBinError(i);
        double et = templ->GetBinError(i);

        double sigma2 = ed * ed + et * et;

        if (sigma2 <= 0.0) continue;

        chi2 += (d - t) * (d - t) / sigma2;
        ++nUsedBins;
    }

    return chi2;
}

// -----------------------------------------------------------------------------
// Main macro
// -----------------------------------------------------------------------------

void closure_test_alphas()
{
    SetFCCStyle();

    gSystem->mkdir("figures", true);

    // -------------------------------------------------------------------------
    // Input files
    // -------------------------------------------------------------------------

    std::string pseudoDataFile = "data/pseudo_data_alphaS_0118.root";

    std::vector<double> alphaSValues = {
        0.110,
        0.114,
        0.118,
        0.122,
        0.126,
        0.130
    };

    std::vector<std::string> templateFiles = {
        "data/template_alphaS_0110.root",
        "data/template_alphaS_0114.root",
        "data/template_alphaS_0118.root",
        "data/template_alphaS_0122.root",
        "data/template_alphaS_0126.root",
        "data/template_alphaS_0130.root"
    };

    // -------------------------------------------------------------------------
    // Histogram and fit-region configuration
    //
    // Histogram range:
    //   0 < 1 - T < 0.30
    //
    // Fit region:
    //   0.02 < 1 - T < 0.25
    //
    // This avoids the endpoint region and the sparse high-tail region.
    // -------------------------------------------------------------------------

    const int nBins = 30;
    const double xMin = 0.0;
    const double xMax = 0.30;

    const double fitMin = 0.02;
    const double fitMax = 0.25;

    std::string selection =
        "oneMinusThrust >= 0.0 && oneMinusThrust < 0.30";

    // -------------------------------------------------------------------------
    // Build pseudo-data histogram
    // -------------------------------------------------------------------------

    TH1D *hData = MakeNormalizedHist(
        pseudoDataFile,
        "h_pseudo_data",
        nBins,
        xMin,
        xMax,
        selection
    );

    if (!hData) {
        std::cerr << "Error: pseudo-data histogram could not be created."
                  << std::endl;
        return;
    }

    // -------------------------------------------------------------------------
    // Build template histograms and compute chi-square values
    // -------------------------------------------------------------------------

    std::vector<TH1D*> hTemplates;
    std::vector<double> chi2Values;
    std::vector<double> chi2Errors;
    std::vector<double> alphaSErrors;

    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Template-based closure test" << std::endl;
    std::cout << "Observable: 1 - T" << std::endl;
    std::cout << "Fit region: " << fitMin << " < 1-T < " << fitMax << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    int nUsedBinsLast = 0;

    for (size_t i = 0; i < templateFiles.size(); ++i) {

        std::string histName = "h_template_" + std::to_string(i);

        TH1D *hTemp = MakeNormalizedHist(
            templateFiles[i],
            histName,
            nBins,
            xMin,
            xMax,
            selection
        );

        if (!hTemp) {
            std::cerr << "Error: template histogram could not be created for "
                      << templateFiles[i] << std::endl;
            return;
        }

        int nUsedBins = 0;
        double chi2 = ComputeChi2(hData, hTemp, fitMin, fitMax, nUsedBins);

        hTemplates.push_back(hTemp);
        chi2Values.push_back(chi2);
        chi2Errors.push_back(0.0);
        alphaSErrors.push_back(0.0);

        nUsedBinsLast = nUsedBins;

        std::cout << std::fixed << std::setprecision(3)
                  << "alphaS = " << alphaSValues[i]
                  << "   chi2 = " << chi2
                  << "   used bins = " << nUsedBins
                  << std::endl;
    }

    // -------------------------------------------------------------------------
    // Find discrete minimum
    // -------------------------------------------------------------------------

    double minChi2 = chi2Values[0];
    double bestAlphaSDiscrete = alphaSValues[0];
    int minIndex = 0;

    for (size_t i = 1; i < chi2Values.size(); ++i) {
        if (chi2Values[i] < minChi2) {
            minChi2 = chi2Values[i];
            bestAlphaSDiscrete = alphaSValues[i];
            minIndex = static_cast<int>(i);
        }
    }

    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Discrete minimum:" << std::endl;
    std::cout << "alphaS_best = " << bestAlphaSDiscrete
              << ", chi2_min = " << minChi2 << std::endl;

    if (minIndex <= 0 || minIndex >= static_cast<int>(alphaSValues.size()) - 1) {
        std::cerr << "Error: chi2 minimum is at the boundary of the scan."
                  << std::endl;
        return;
    }

    // -------------------------------------------------------------------------
    // Use only the three points closest to the minimum for interpolation
    // -------------------------------------------------------------------------

    std::vector<double> alphaSFitValues = {
        alphaSValues[minIndex - 1],
        alphaSValues[minIndex],
        alphaSValues[minIndex + 1]
    };

    std::vector<double> chi2FitValues = {
        chi2Values[minIndex - 1],
        chi2Values[minIndex],
        chi2Values[minIndex + 1]
    };

    std::vector<double> fitErrors = {
        0.0,
        0.0,
        0.0
    };

    TGraphErrors *graphFit = new TGraphErrors(
        alphaSFitValues.size(),
        alphaSFitValues.data(),
        chi2FitValues.data(),
        fitErrors.data(),
        fitErrors.data()
    );

    // -------------------------------------------------------------------------
    // Parabolic fit:
    //
    //   chi2(alphaS) = curvature * (alphaS - alphaS_fit)^2 + chi2_min
    // -------------------------------------------------------------------------

    TF1 *parabola = new TF1(
        "parabola",
        "[0]*(x-[1])*(x-[1]) + [2]",
        alphaSFitValues.front() - 0.001,
        alphaSFitValues.back() + 0.001
    );

    parabola->SetParameter(0, 1.0e7);
    parabola->SetParameter(1, bestAlphaSDiscrete);
    parabola->SetParameter(2, minChi2);

    parabola->SetParName(0, "curvature");
    parabola->SetParName(1, "alphaS_fit");
    parabola->SetParName(2, "chi2_min");

    graphFit->Fit(parabola, "RQ");

    double curvature = parabola->GetParameter(0);
    double alphaSFit = parabola->GetParameter(1);
    double chi2FitMin = parabola->GetParameter(2);

    // Statistical uncertainty from Delta chi2 = 1
    double alphaSFitErr = -1.0;

    if (curvature > 0.0) {
        alphaSFitErr = 1.0 / std::sqrt(curvature);
    }

    std::cout << "Parabolic interpolation:" << std::endl;
    std::cout << "alphaS_fit = " << std::fixed << std::setprecision(5)
              << alphaSFit << " +/- " << alphaSFitErr
              << "  (Delta chi2 = 1)" << std::endl;
    std::cout << "chi2_min_fit = " << chi2FitMin << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    // -------------------------------------------------------------------------
    // Build zoomed Delta chi2 curve around the fitted minimum.
    //
    // This plot intentionally shows only the local uncertainty region.
    // The full template chi2 values are printed in the terminal output.
    // -------------------------------------------------------------------------

    const double xPlotMin = alphaSFit - 0.00120;
    const double xPlotMax = alphaSFit + 0.00120;

    const int nScanPoints = 400;

    std::vector<double> xZoom;
    std::vector<double> yZoom;

    for (int i = 0; i < nScanPoints; ++i) {
        double x = xPlotMin + (xPlotMax - xPlotMin) * i / (nScanPoints - 1.0);
        double y = curvature * (x - alphaSFit) * (x - alphaSFit);

        xZoom.push_back(x);
        yZoom.push_back(y);
    }

    TGraph *graphDelta = new TGraph(
        xZoom.size(),
        xZoom.data(),
        yZoom.data()
    );

    graphDelta->SetName("g_delta_chi2_zoom");
    graphDelta->SetTitle("");
    graphDelta->SetLineColor(kRed + 1);
    graphDelta->SetLineWidth(2);

    // Best-fit point
    TGraph *graphBest = new TGraph(1);
    graphBest->SetPoint(0, alphaSFit, 0.0);
    graphBest->SetMarkerStyle(20);
    graphBest->SetMarkerSize(1.2);
    graphBest->SetMarkerColor(kBlack);

    // Delta chi2 = 1 crossing points
    TGraph *graphErrPoints = new TGraph(2);
    graphErrPoints->SetPoint(0, alphaSFit - alphaSFitErr, 1.0);
    graphErrPoints->SetPoint(1, alphaSFit + alphaSFitErr, 1.0);
    graphErrPoints->SetMarkerStyle(24);
    graphErrPoints->SetMarkerSize(1.1);
    graphErrPoints->SetMarkerColor(kRed + 1);

    // -------------------------------------------------------------------------
    // Canvas
    // -------------------------------------------------------------------------

    TCanvas *canvas = new TCanvas(
        "canvas",
        "closure test alphaS",
        1050,
        800
    );

    canvas->SetLeftMargin(0.17);
    canvas->SetRightMargin(0.08);
    canvas->SetTopMargin(0.09);
    canvas->SetBottomMargin(0.14);

    graphDelta->GetXaxis()->SetTitle("#alpha_{s}");
    graphDelta->GetYaxis()->SetTitle("#Delta#chi^{2}");
    graphDelta->GetXaxis()->SetLimits(xPlotMin, xPlotMax);

    graphDelta->SetMinimum(0.0);
    graphDelta->SetMaximum(6.5);

    graphDelta->GetXaxis()->SetTitleSize(0.040);
    graphDelta->GetXaxis()->SetLabelSize(0.034);
    graphDelta->GetXaxis()->SetTitleOffset(1.05);

    graphDelta->GetYaxis()->SetTitleSize(0.040);
    graphDelta->GetYaxis()->SetLabelSize(0.034);
    graphDelta->GetYaxis()->SetTitleOffset(1.55);

    graphDelta->Draw("AL");

    // Horizontal Delta chi2 = 1 line
    TLine *lineDelta1 = new TLine(xPlotMin, 1.0, xPlotMax, 1.0);
    lineDelta1->SetLineColor(kGray + 2);
    lineDelta1->SetLineStyle(2);
    lineDelta1->SetLineWidth(2);
    lineDelta1->Draw("same");

    // Central fit line
    TLine *lineFit = new TLine(alphaSFit, 0.0, alphaSFit, 1.0);
    lineFit->SetLineColor(kRed + 1);
    lineFit->SetLineStyle(2);
    lineFit->SetLineWidth(2);
    lineFit->Draw("same");

    // Error interval lines
    TLine *lineFitLow = new TLine(
        alphaSFit - alphaSFitErr,
        0.0,
        alphaSFit - alphaSFitErr,
        1.0
    );

    lineFitLow->SetLineColor(kRed + 1);
    lineFitLow->SetLineStyle(3);
    lineFitLow->SetLineWidth(2);
    lineFitLow->Draw("same");

    TLine *lineFitHigh = new TLine(
        alphaSFit + alphaSFitErr,
        0.0,
        alphaSFit + alphaSFitErr,
        1.0
    );

    lineFitHigh->SetLineColor(kRed + 1);
    lineFitHigh->SetLineStyle(3);
    lineFitHigh->SetLineWidth(2);
    lineFitHigh->Draw("same");

    graphBest->Draw("P same");
    graphErrPoints->Draw("P same");

    // -------------------------------------------------------------------------
    // Compact FCC-ee label
    //
    // For the zoomed closure-test plot, the long process label is intentionally
    // omitted from the canvas to avoid overlap with the parabola. It should be
    // stated in the figure caption.
    // -------------------------------------------------------------------------

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);

    label.SetTextSize(0.032);
    label.DrawLatex(0.20, 0.86, "#bf{FCC-ee Simulation}");

    // -------------------------------------------------------------------------
    // Fit result text
    // -------------------------------------------------------------------------

    TLatex resultText;
    resultText.SetNDC();
    resultText.SetTextFont(42);
    resultText.SetTextSize(0.022);

    resultText.DrawLatex(
        0.69,
        0.86,
        "#alpha_{s}^{true} = 0.118"
    );

    resultText.DrawLatex(
        0.69,
        0.825,
        Form("#alpha_{s}^{fit} = %.5f #pm %.5f", alphaSFit, alphaSFitErr)
    );

    resultText.DrawLatex(
        0.69,
        0.790,
        "0.02 < 1-T < 0.25"
    );

    resultText.DrawLatex(
        0.69,
        0.755,
        "Zoomed #Delta#chi^{2} around minimum"
    );

    // -------------------------------------------------------------------------
    // Legend
    // -------------------------------------------------------------------------

    TLegend *legend = new TLegend(0.69, 0.57, 0.91, 0.70);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.022);
    legend->AddEntry(graphDelta, "Parabolic #Delta#chi^{2}", "l");
    legend->AddEntry(lineDelta1, "#Delta#chi^{2}=1", "l");
    legend->AddEntry(graphBest, "Best fit", "p");
    legend->Draw();

    // -------------------------------------------------------------------------
    // Save output
    // -------------------------------------------------------------------------

    canvas->SaveAs("figures/closure_test_alphas.pdf");
    canvas->SaveAs("figures/closure_test_alphas.png");

    std::cout << "Output written:" << std::endl;
    std::cout << "  figures/closure_test_alphas.pdf" << std::endl;
    std::cout << "  figures/closure_test_alphas.png" << std::endl;
}