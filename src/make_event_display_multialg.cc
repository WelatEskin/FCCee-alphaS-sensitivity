#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"

#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

using namespace Pythia8;
using namespace fastjet;

// ============================================================
// Helper: assign jet id to each input particle from FastJet jets
// ============================================================

std::vector<int> assign_jet_ids_from_constituents(
    const std::vector<PseudoJet>& jets,
    int n_particles
) {
    std::vector<int> jet_id(n_particles, -1);

    for (int j = 0; j < (int)jets.size(); ++j) {
        std::vector<PseudoJet> constituents = jets[j].constituents();

        for (const auto& c : constituents) {
            int idx = c.user_index();

            if (idx >= 0 && idx < n_particles) {
                jet_id[idx] = j;
            }
        }
    }

    return jet_id;
}


// ============================================================
// Main
// ============================================================

int main(int argc, char* argv[]) {

    // ------------------------------------------------------------
    // Configuration
    // ------------------------------------------------------------

    const int nEvents = 1000000;

    const double eCM = 91.2;
    const double alphaS = 0.118;

    const double y_cut = 0.01;

    // Used for ee_genkt_algorithm.
    // This is the same R-like parameter used in your anti-kT-like setup.
    const double R = 0.4;

    const char* output_name = "data/event_display_alphaS_0118_multialg.root";

    if (argc > 1) {
        output_name = argv[1];
    }

    std::cout << "----------------------------------------" << std::endl;
    std::cout << "FCC-ee multi-algorithm event-display ROOT producer" << std::endl;
    std::cout << "Output file: " << output_name << std::endl;
    std::cout << "sqrt(s)    : " << eCM << " GeV" << std::endl;
    std::cout << "alphaS     : " << alphaS << std::endl;
    std::cout << "nEvents    : " << nEvents << std::endl;
    std::cout << "y_cut      : " << y_cut << std::endl;
    std::cout << "R          : " << R << std::endl;
    std::cout << "----------------------------------------" << std::endl;


    // ------------------------------------------------------------
    // ROOT output
    // ------------------------------------------------------------

    TFile* outFile = new TFile(output_name, "RECREATE");
    TTree* tree = new TTree("events", "FCC-ee particle-level event display tree");

    int event_id;
    int particle_id;
    int pdg_id;
    int charge3;

    double px;
    double py;
    double pz;
    double E;
    double pt;
    double eta;
    double rapidity;
    double phi;

    int durham_jet_id;
    int antikt_jet_id;
    int kt_jet_id;
    int cambridge_jet_id;

    tree->Branch("event_id", &event_id, "event_id/I");
    tree->Branch("particle_id", &particle_id, "particle_id/I");
    tree->Branch("pdg_id", &pdg_id, "pdg_id/I");
    tree->Branch("charge3", &charge3, "charge3/I");

    tree->Branch("px", &px, "px/D");
    tree->Branch("py", &py, "py/D");
    tree->Branch("pz", &pz, "pz/D");
    tree->Branch("E", &E, "E/D");
    tree->Branch("pt", &pt, "pt/D");
    tree->Branch("eta", &eta, "eta/D");
    tree->Branch("rapidity", &rapidity, "rapidity/D");
    tree->Branch("phi", &phi, "phi/D");

    tree->Branch("durham_jet_id", &durham_jet_id, "durham_jet_id/I");
    tree->Branch("antikt_jet_id", &antikt_jet_id, "antikt_jet_id/I");
    tree->Branch("kt_jet_id", &kt_jet_id, "kt_jet_id/I");
    tree->Branch("cambridge_jet_id", &cambridge_jet_id, "cambridge_jet_id/I");


    // ------------------------------------------------------------
    // Pythia8 setup
    // ------------------------------------------------------------

    Pythia pythia;
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = 118"); 
    
    pythia.readString("Beams:idA = 11");
    pythia.readString("Beams:idB = -11");
    pythia.readString("Beams:eCM = 91.2");

    // e+e- -> gamma*/Z -> q qbar
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");

    // Z -> hadrons only
    pythia.readString("23:onMode = off");
    pythia.readString("23:onIfAny = 1 2 3 4 5");

    // Hadron-level final state
    pythia.readString("HadronLevel:all = on");

    // Alpha_s reference setting.
    // These two settings are commonly used for shower alpha_s control.
    pythia.readString("TimeShower:alphaSvalue = 0.118");
    pythia.readString("SpaceShower:alphaSvalue = 0.118");


    pythia.init();


    // ------------------------------------------------------------
    // FastJet definitions
    // ------------------------------------------------------------

    // Durham e+e- kT
    JetDefinition durham_def(ee_kt_algorithm);

    // Generalized e+e- kT family:
    // p = -1 : anti-kT-like
    // p =  1 : kT-like
    // p =  0 : Cambridge/Aachen-like
    JetDefinition antikt_def(ee_genkt_algorithm, R, -1.0);
    JetDefinition kt_def(ee_genkt_algorithm, R, 1.0);
    JetDefinition cambridge_def(ee_genkt_algorithm, R, 0.0);


    // ------------------------------------------------------------
    // Event loop
    // ------------------------------------------------------------

    int written_events = 0;
    int written_particles = 0;

    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {

        if (!pythia.next()) {
            continue;
        }

        std::vector<PseudoJet> particles;

        std::vector<int> particle_pdg_id;
        std::vector<int> particle_charge3;

        std::vector<double> particle_px;
        std::vector<double> particle_py;
        std::vector<double> particle_pz;
        std::vector<double> particle_E;

        // --------------------------------------------------------
        // Collect visible final-state particles
        // --------------------------------------------------------

        for (int i = 0; i < pythia.event.size(); ++i) {

            Particle& p = pythia.event[i];

            if (!p.isFinal()) continue;
            if (!p.isVisible()) continue;

            double px_i = p.px();
            double py_i = p.py();
            double pz_i = p.pz();
            double E_i  = p.e();

            PseudoJet pj(px_i, py_i, pz_i, E_i);
            pj.set_user_index((int)particles.size());

            particles.push_back(pj);

            particle_pdg_id.push_back(p.id());
            particle_charge3.push_back(p.chargeType());

            particle_px.push_back(px_i);
            particle_py.push_back(py_i);
            particle_pz.push_back(pz_i);
            particle_E.push_back(E_i);
        }

        if (particles.size() < 2) {
            continue;
        }


        // --------------------------------------------------------
        // Cluster with all algorithms
        // --------------------------------------------------------

        ClusterSequence cs_durham(particles, durham_def);
        ClusterSequence cs_antikt(particles, antikt_def);
        ClusterSequence cs_kt(particles, kt_def);
        ClusterSequence cs_cambridge(particles, cambridge_def);

        std::vector<PseudoJet> jets_durham =
            sorted_by_E(cs_durham.exclusive_jets_ycut(y_cut));

        std::vector<PseudoJet> jets_antikt =
            sorted_by_E(cs_antikt.exclusive_jets_ycut(y_cut));

        std::vector<PseudoJet> jets_kt =
            sorted_by_E(cs_kt.exclusive_jets_ycut(y_cut));

        std::vector<PseudoJet> jets_cambridge =
            sorted_by_E(cs_cambridge.exclusive_jets_ycut(y_cut));


        // --------------------------------------------------------
        // Assign particle -> jet labels
        // --------------------------------------------------------

        std::vector<int> durham_ids =
            assign_jet_ids_from_constituents(jets_durham, particles.size());

        std::vector<int> antikt_ids =
            assign_jet_ids_from_constituents(jets_antikt, particles.size());

        std::vector<int> kt_ids =
            assign_jet_ids_from_constituents(jets_kt, particles.size());

        std::vector<int> cambridge_ids =
            assign_jet_ids_from_constituents(jets_cambridge, particles.size());


        // --------------------------------------------------------
        // Fill flat particle-level tree
        // --------------------------------------------------------

        event_id = iEvent;

        for (int ip = 0; ip < (int)particles.size(); ++ip) {

            particle_id = ip;

            pdg_id = particle_pdg_id[ip];
            charge3 = particle_charge3[ip];

            px = particle_px[ip];
            py = particle_py[ip];
            pz = particle_pz[ip];
            E  = particle_E[ip];

            pt = std::sqrt(px * px + py * py);
            phi = std::atan2(py, px);

            double p_abs = std::sqrt(px * px + py * py + pz * pz);

            if (p_abs != std::abs(pz)) {
                eta = 0.5 * std::log((p_abs + pz) / (p_abs - pz));
            } else {
                eta = (pz >= 0.0) ? 1.0e10 : -1.0e10;
            }

            if (E != std::abs(pz)) {
                rapidity = 0.5 * std::log((E + pz) / (E - pz));
            } else {
                rapidity = (pz >= 0.0) ? 1.0e10 : -1.0e10;
            }

            durham_jet_id = durham_ids[ip];
            antikt_jet_id = antikt_ids[ip];
            kt_jet_id = kt_ids[ip];
            cambridge_jet_id = cambridge_ids[ip];

            tree->Fill();

            written_particles++;
        }

        written_events++;

        if (iEvent % 1000 == 0) {
            std::cout
                << "Processed event " << iEvent
                << " | visible particles = " << particles.size()
                << " | Njet Durham = " << jets_durham.size()
                << " | anti-kT-like = " << jets_antikt.size()
                << " | kT-like = " << jets_kt.size()
                << " | Cambridge-like = " << jets_cambridge.size()
                << std::endl;
        }
    }


    // ------------------------------------------------------------
    // Write output
    // ------------------------------------------------------------

    outFile->cd();
    tree->Write();
    outFile->Close();

    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Done." << std::endl;
    std::cout << "Written events    : " << written_events << std::endl;
    std::cout << "Written particles : " << written_particles << std::endl;
    std::cout << "Output file       : " << output_name << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    return 0;
}