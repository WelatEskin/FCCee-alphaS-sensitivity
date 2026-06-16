// generate_fccee_alphas_root.cc
//
// Generator-level FCC-ee alphaS sensitivity study
//
// Process:
//   e+ e- -> Z/gamma* -> q qbar -> hadrons
//
// Workflow:
//   Pythia8 event generation
//   FastJet jet reconstruction
//   ROOT TTree output
//
// Usage:
//   ./generate_fccee_alphas_root <alphaS> <nEvents> <seed> <output.root>
//
// Example:
//   ./generate_fccee_alphas_root 0.118 10000 12345 data/validation_alphaS_0118.root

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include "Pythia8/Pythia.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include "TFile.h"
#include "TTree.h"

using namespace Pythia8;

// -----------------------------------------------------------------------------
// Simple thrust calculation
//
// Thrust is defined as:
//   T = max_n sum_i |p_i . n| / sum_i |p_i|
//
// A fully precise thrust calculation requires an optimization over all possible
// axes. For this generator-level study, we estimate the thrust axis by scanning
// over the directions of the visible final-state particles. This approximation
// is sufficient for a first closure-test framework and can later be replaced by
// a more specialized thrust minimizer if needed.
// -----------------------------------------------------------------------------

double computeThrust(const std::vector<fastjet::PseudoJet>& particles)
{
    if (particles.empty()) return -1.0;

    double denominator = 0.0;

    for (const auto& p : particles) {
        const double pAbs = std::sqrt(
            p.px() * p.px() +
            p.py() * p.py() +
            p.pz() * p.pz()
        );
        denominator += pAbs;
    }

    if (denominator <= 0.0) return -1.0;

    double bestThrust = 0.0;

    for (const auto& axisCandidate : particles) {
        const double axisNorm = std::sqrt(
            axisCandidate.px() * axisCandidate.px() +
            axisCandidate.py() * axisCandidate.py() +
            axisCandidate.pz() * axisCandidate.pz()
        );

        if (axisNorm <= 0.0) continue;

        const double nx = axisCandidate.px() / axisNorm;
        const double ny = axisCandidate.py() / axisNorm;
        const double nz = axisCandidate.pz() / axisNorm;

        double numerator = 0.0;

        for (const auto& p : particles) {
            const double projection = p.px() * nx + p.py() * ny + p.pz() * nz;
            numerator += std::abs(projection);
        }

        const double thrust = numerator / denominator;

        if (thrust > bestThrust) {
            bestThrust = thrust;
        }
    }

    return bestThrust;
}

// -----------------------------------------------------------------------------
// Helper function: get leading jet energy
// -----------------------------------------------------------------------------

double getLeadingJetEnergy(const std::vector<fastjet::PseudoJet>& jets)
{
    if (jets.empty()) return -1.0;

    double maxEnergy = -1.0;

    for (const auto& jet : jets) {
        if (jet.E() > maxEnergy) {
            maxEnergy = jet.E();
        }
    }

    return maxEnergy;
}


// -----------------------------------------------------------------------------
// Main program
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    // -------------------------------------------------------------------------
    // Read command-line arguments
    // -------------------------------------------------------------------------

    if (argc != 5) {
        std::cerr << "Usage:\n"
                  << "  " << argv[0]
                  << " <alphaS> <nEvents> <seed> <output.root>\n\n"
                  << "Example:\n"
                  << "  " << argv[0]
                  << " 0.118 10000 12345 data/validation_alphaS_0118.root\n";
        return 1;
    }

    const double inputAlphaS = std::stod(argv[1]);
    const int nEvents = std::stoi(argv[2]);
    const int seed = std::stoi(argv[3]);
    const std::string outputFileName = argv[4];

    std::cout << "----------------------------------------\n";
    std::cout << "FCC-ee alphaS generator-level production\n";
    std::cout << "alphaS input : " << inputAlphaS << "\n";
    std::cout << "nEvents      : " << nEvents << "\n";
    std::cout << "seed         : " << seed << "\n";
    std::cout << "output file  : " << outputFileName << "\n";
    std::cout << "----------------------------------------\n";

    // -------------------------------------------------------------------------
    // ROOT output setup
    // -------------------------------------------------------------------------

    TFile* outFile = new TFile(outputFileName.c_str(), "RECREATE");

    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: could not create output file "
                  << outputFileName << "\n";
        return 1;
    }

    TTree* tree = new TTree("Events", "FCC-ee alphaS event-level observables");

    int event_id = -1;
    double alphaS = inputAlphaS;

    double thrust = -1.0;
    double oneMinusThrust = -1.0;

    int nJets_durham = -1;
    int nJets_antikt = -1;

    double leadingJetEnergy_durham = -1.0;
    double leadingJetEnergy_antikt = -1.0;

    double y23 = -1.0;

    tree->Branch("event_id", &event_id, "event_id/I");
    tree->Branch("alphaS", &alphaS, "alphaS/D");

    tree->Branch("thrust", &thrust, "thrust/D");
    tree->Branch("oneMinusThrust", &oneMinusThrust, "oneMinusThrust/D");

    tree->Branch("nJets_durham", &nJets_durham, "nJets_durham/I");
    tree->Branch("nJets_antikt", &nJets_antikt, "nJets_antikt/I");

    tree->Branch("leadingJetEnergy_durham",
                 &leadingJetEnergy_durham,
                 "leadingJetEnergy_durham/D");

    tree->Branch("leadingJetEnergy_antikt",
                 &leadingJetEnergy_antikt,
                 "leadingJetEnergy_antikt/D");

    tree->Branch("y23", &y23, "y23/D");

    // -------------------------------------------------------------------------
    // Pythia8 setup
    // -------------------------------------------------------------------------

    Pythia pythia;

    // Beam setup:
    // idA = 11  -> electron
    // idB = -11 -> positron
    // eCM = 91.2 GeV corresponds approximately to the Z pole.
    pythia.readString("Beams:idA = 11");
    pythia.readString("Beams:idB = -11");
    pythia.readString("Beams:eCM = 91.2");

    // Hard process:
    // e+ e- -> gamma*/Z -> f fbar
    // At the Z pole, this generates the desired hadronic process once the
    // Z decay channels are restricted to quarks.
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");

    // Restrict Z decays to hadronic final states:
    // 1 = d, 2 = u, 3 = s, 4 = c, 5 = b
    // Top quarks are not included because sqrt(s)=91.2 GeV is below ttbar threshold.
    pythia.readString("23:onMode = off");
    pythia.readString("23:onIfAny = 1 2 3 4 5");

    // Strong coupling setting for final-state shower.
    // This is the main parameter varied in the template scan.
    pythia.readString("TimeShower:alphaSvalue = " + std::to_string(inputAlphaS));    pythia.readString("TimeShower:alphaSorder = 2");

    // Nominal shower renormalization scale.
    // Scale variations will be studied separately later.
    pythia.readString("TimeShower:renormMultFac = 1.0");

    // Random seed for statistically independent samples.
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(seed));

    // Initialize Pythia.
    if (!pythia.init()) {
        std::cerr << "Error: Pythia initialization failed.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // FastJet setup
    // -------------------------------------------------------------------------

    // Durham e+e- kT algorithm:
    // This is the main jet algorithm for e+e- event-shape and jet-rate studies.
    fastjet::JetDefinition durham_def(fastjet::ee_kt_algorithm);

    // e+e- generalized anti-kT algorithm:
    // Used as an algorithm-dependence cross-check.
    // The parameter p = -1 corresponds to anti-kT-like clustering.
    fastjet::JetDefinition antikt_def(fastjet::ee_genkt_algorithm, 0.4, -1.0);

    // Energy threshold used only for the anti-kT-like inclusive jet collection.
    // Durham jets are defined separately through a y_cut resolution parameter.
    const double jetEnergyMin = 1.0; // GeV

    // -------------------------------------------------------------------------
    // Event loop
    // -------------------------------------------------------------------------

    int nAccepted = 0;

    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {

        if (!pythia.next()) continue;

        event_id = iEvent;

        std::vector<fastjet::PseudoJet> visibleParticles;

        // Select visible final-state particles.
        for (int i = 0; i < pythia.event.size(); ++i) {

            if (!pythia.event[i].isFinal()) continue;
            if (!pythia.event[i].isVisible()) continue;

            const double px = pythia.event[i].px();
            const double py = pythia.event[i].py();
            const double pz = pythia.event[i].pz();
            const double e  = pythia.event[i].e();

            // Avoid pathological zero-energy entries.
            if (e <= 0.0) continue;

            fastjet::PseudoJet particle(px, py, pz, e);
            visibleParticles.push_back(particle);
        }

        // Require at least two visible particles to define event-shape observables.
        if (visibleParticles.size() < 2) continue;

        // ---------------------------------------------------------------------
        // Thrust and 1 - thrust
        // ---------------------------------------------------------------------

        thrust = computeThrust(visibleParticles);
        oneMinusThrust = 1.0 - thrust;

        if (thrust <= 0.0 || thrust > 1.0001) continue;

        // ---------------------------------------------------------------------
        // Durham jet reconstruction
        // ---------------------------------------------------------------------

        fastjet::ClusterSequence cs_durham(visibleParticles, durham_def);

        // Durham e+e- jet multiplicity should be defined through a y_cut
        // resolution parameter rather than through an energy threshold.
        // The value y_cut = 0.002 is used as a nominal working point.
        // It can later be varied as a jet-definition systematic check.
        const double yCutDurham = 0.002;

        nJets_durham = cs_durham.n_exclusive_jets_ycut(yCutDurham);

        std::vector<fastjet::PseudoJet> jets_durham =
            fastjet::sorted_by_E(cs_durham.exclusive_jets_ycut(yCutDurham));

        leadingJetEnergy_durham = getLeadingJetEnergy(jets_durham);

        // Durham y23:
        // exclusive_ymerge(3) gives the transition scale associated with
        //merging from 3 jets to 2 jets in the Durham clustering sequence.
        if (visibleParticles.size() >= 3) {
             y23 = cs_durham.exclusive_ymerge(3);
        } else {
            y23 = -1.0;
}

        // ---------------------------------------------------------------------
        // e+e- anti-kT-like reconstruction
        // ---------------------------------------------------------------------

        fastjet::ClusterSequence cs_antikt(visibleParticles, antikt_def);

        // For the e+e- generalized anti-kT algorithm, we first retrieve the inclusive
        // jet collection with a very loose internal cut and then apply the physical
        // jet-energy threshold manually. This avoids possible ambiguities in the
        // interpretation of the inclusive_jets argument for e+e- clustering.
        std::vector<fastjet::PseudoJet> jets_antikt_all =
            fastjet::sorted_by_E(cs_antikt.inclusive_jets(0.0));

        std::vector<fastjet::PseudoJet> jets_antikt;

        for (const auto& jet : jets_antikt_all) {
             if (jet.E() > jetEnergyMin) {
                jets_antikt.push_back(jet);
            }
}

jets_antikt = fastjet::sorted_by_E(jets_antikt);

nJets_antikt = static_cast<int>(jets_antikt.size());
leadingJetEnergy_antikt = getLeadingJetEnergy(jets_antikt);

        tree->Fill();
        ++nAccepted;

        // Progress printout for long productions.
        if ((iEvent + 1) % 10000 == 0) {
            std::cout << "Processed " << (iEvent + 1)
                      << " / " << nEvents
                      << " events, accepted " << nAccepted << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Finalize
    // -------------------------------------------------------------------------

    pythia.stat();

    outFile->cd();
    tree->Write();
    outFile->Close();

    std::cout << "----------------------------------------\n";
    std::cout << "Production completed.\n";
    std::cout << "Generated events : " << nEvents << "\n";
    std::cout << "Accepted events  : " << nAccepted << "\n";
    std::cout << "Output written   : " << outputFileName << "\n";
    std::cout << "----------------------------------------\n";

    return 0;
}