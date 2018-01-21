// An analyzer to produce key selection criteria plots for Z' -> ditau processes.
// arXiv:1707.07016v1
// USAGE:
// gSystem->Load("<path_to_delphes>/Delphes/libDelphes.so");
// .x DiTauAnalyzer.C("input_file.root", "Z' 500 GeV", nbins);

#ifdef __CLING__
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "lester_mt2_bisect.h"
#include <utility> // std::pair, std::make_pair
#else
class ExRootTreeReader;
class ExRootResult;
#endif


// Declare global variables.
Int_t accepted_events = 0;


void DiTauAnalyzer(const char *file_name, const char *sample_desc, int nbins,
                   bool apply_cuts = false) {
  gSystem->Load("libDelphes.so");

  // Create a chain of root trees.
  TChain chain("Delphes");
  chain.Add(file_name);

  // Create object of class ExRootTreeReader.
  ExRootTreeReader *tree_reader = new ExRootTreeReader(&chain);
  Long64_t number_of_entries = tree_reader->GetEntries();

  // Get pointers to branches used in this analysis.
  TClonesArray *branch_jet = tree_reader->UseBranch("Jet");
  TClonesArray *branch_met = tree_reader->UseBranch("MissingET");

  // Book histograms.
  TH1F *hist_pair_mass = new TH1F("pair_mass", "Max{M(tau,j)}", nbins, 0., 500.);
  TH1F *hist_MET = new TH1F("met", "Normalized Missing ET", nbins, 0., 2.);
  TH1F *hist_MT2 = new TH1F("tmass", "MT2 (Ditau + MET)", nbins, 0., 150.);
  TH1F *hist_HT_LT = new TH1F("HT_LT", "HT - LT", nbins, -500., 500.);
  TH1F *hist_ditau_mass = new TH1F("ditau_mass", "M(tau+,tau-)", nbins, 0., 500.);
  TH1F *hist_pt_btag = new TH1F("pt_btag", "BTag Pt", nbins, 0., 300.);
  TH1F *hist_pt_jet = new TH1F("pt_jet", "Secondary Jet Pt", nbins, 0., 300.);
  TH1F *hist_pt_tau_p = new TH1F("pt_tau_p", "Tau+ Pt", nbins, 0., 300.);
  TH1F *hist_pt_tau_m = new TH1F("pt_tau_m", "Tau- Pt", nbins, 0., 300.);
  TH1F *hist_pt_ditau = new TH1F("pt_ditau", "Ditau Pt", nbins, 0., 300.);
  TH1F *hist_e_btag = new TH1F("e_btag", "BTag Jet E", nbins, 0., 300.);
  TH1F *hist_e_jet = new TH1F("e_jet", "Secondary Jet E", nbins, 0., 300.);
  TH1F *hist_e_tau_p = new TH1F("e_tau_p", "E(#tau^{+}", nbins, 0., 300.);
  TH1F *hist_e_tau_m = new TH1F("e_tau_m", "E(#tau^{-})", nbins, 0., 300.);
  TH1F *hist_e_ratio = new TH1F("e_ratio", "", nbins, 0., 1.);
  TH1F *hist_unboosted_MT2 = new TH1F("unboost_mt2", "", nbins, 0., 150.);
  TH1F *hist_denis = new TH1F("hist_denis", "", nbins, -.3, .3);
  TH1F *hist_topology = new TH1F("hist_topo", "", nbins, -3.14, 3.14);
  TH1F *hist_j_topology = new TH1F("hist_unboosted_topo", "", nbins, -3.14, 3.14);
  TH1F *hist_deltaR_jets = new TH1F("hist_deltaR_1", "", nbins, 0., 6.);
  TH1F *hist_deltaR_taus = new TH1F("hist_deltaR_2", "", nbins, 0., 6.);
  TH1F *hist_deltaR_tau_jet = new TH1F("hist_deltaR_3", "", nbins, 0., 6.);
  TH1F *hist_deltaR_tau_met = new TH1F("hist_deltaR_4", "", nbins, 0., 6.);
  TH1F *hist_deltaR_topology = new TH1F("hist_deltaR_topology", "", nbins, -6., 6.);

  // Book 2D histograms for exploratory analysis.
  TH2F *hist2d_topology = new TH2F("topo", "WJets(ditau) Topology", nbins,
                                 0., 50., nbins, -3.14, 3.14);
  hist2d_topology->GetXaxis()->SetTitle("M_{T2}");
  hist2d_topology->GetYaxis()->SetTitle("max[#delta#phi] - #delta#phi(#tau^{+}, #tau^{-})");

  // E(tau+) vs. E(tau-).
  TH2F *hist2d_e_tau = new TH2F("e_tautau", "E(#tau^{+}) vs. E(#tau^{-})",
                                nbins, 0., 500., nbins, 0., 500.);
  hist2d_e_tau->GetXaxis()->SetTitle("E(#tau^{+}) [GeV]");
  hist2d_e_tau->GetYaxis()->SetTitle("E(#tau^{-}) [GeV]");

  // E vs Pt for Tau+
  TH2F *hist2d_e_pt_tau = new TH2F("e_pt_tau", "E(#tau) vs. P_{T}(#tau)",
                                nbins, 0., 500., nbins, 0., 300.);
  hist2d_e_pt_tau->GetXaxis()->SetTitle("E(#tau) [GeV]");
  hist2d_e_pt_tau->GetYaxis()->SetTitle("P_{T}(#tau) [GeV]");

  // Pt vs. Pt for tau+/tau-
  TH2F *hist2d_pt_tau = new TH2F("pt_tautau", "P_{T}(#tau^{+}) vs. P_{T}(#tau^{-})",
                                 nbins, 0., 300., nbins, 0., 300.);
  hist2d_pt_tau->GetXaxis()->SetTitle("P_{T}(#tau^{+}) [GeV]");
  hist2d_pt_tau->GetYaxis()->SetTitle("P_{T}(#tau^{-}) [GeV]");

  // MT2 vs. PT(tau+)
  TH2F *hist2d_e_pt_btag = new TH2F("mt2_pt_tau", "MT2 vs. P_{T}(#tau+)",
                                    nbins, 0., 100., nbins, 0., 300.);
  hist2d_e_pt_btag->GetXaxis()->SetTitle("MT2");
  hist2d_e_pt_btag->GetYaxis()->SetTitle("P_{T}(#tau) [GeV]");

  // Pt vs. Pt for MATCHED BTag and Tau
  TH2F *hist2d_pt_btag_tau = new TH2F("pt_btag_tau", "P_{T}(b) vs. P_{T}(#tau)",
                                      nbins, 0., 300., nbins, 0., 300.);
  hist2d_pt_btag_tau->GetXaxis()->SetTitle("P_{T}(b) [GeV]");
  hist2d_pt_btag_tau->GetYaxis()->SetTitle("P_{T}(#tau) [GeV]");

  // Eta vs. Eta for MATCHED BTag and Tau
  TH2F *hist2d_eta_btag_tau = new TH2F("eta_btag_tau", "#eta(b) vs. #eta(#tau)",
                                       nbins, 0., 4., nbins, 0., 4.);
  hist2d_eta_btag_tau->GetXaxis()->SetTitle("#eta(b)");
  hist2d_eta_btag_tau->GetYaxis()->SetTitle("#eta(#tau)");


  // Main event loop.
  for(Int_t entry = 0; entry < number_of_entries; ++entry) {

    tree_reader->ReadEntry(entry);

    // Declare physics objects.
    Muon *muon;
    Jet *jet;
    MissingET *met;
    TLorentzVector btag_jet;
    TLorentzVector second_jet;
    TLorentzVector tau_plus;
    TLorentzVector tau_minus;

    // Declare running indices.
    Int_t jet_size = branch_jet->GetEntries();
    Int_t met_size = branch_met->GetEntries();

    // Require 2 jets, at least one of them BTagged, and 2 OS Tau's
    // (which also count as Delphes jets)..
    if (jet_size < 4) continue;
    bool found_btag = false;
    bool found_tau_plus = false;
    bool found_tau_minus = false;

    Int_t leading_btag_id = -1;
    // Loop over jets and find the highest pt BTag, the second highest pt, non-TauTag jet,
    // and two opposite sign TauTag jets.
    for (Int_t ii = 0; ii < jet_size; ii++) {
      jet = (Jet*) branch_jet->At(ii);
      // Find the leading b-tagged jet.
      if (jet->BTag == 1 && jet->PT > btag_jet.Pt()) {
        btag_jet.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
        leading_btag_id = ii;
        found_btag = true;
      }
      // Find the leading OS Tau pair.
      if (jet->TauTag == 1) {
        if (jet->Charge == 1 && jet->PT > tau_plus.Pt()) {
          tau_plus.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
          found_tau_plus = true;
        }
        if (jet->Charge == -1 && jet->PT > tau_minus.Pt()) {
          tau_minus.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
          found_tau_minus = true;
        }
      }
    }

    // Loop to find the other highest-pt jet with ID different than the leading BTag.
    for (Int_t ii = 0; ii < jet_size; ii++) {
      jet = (Jet*) branch_jet->At(ii);
      if (jet->TauTag == 0 && jet->PT > second_jet.Pt() && ii != leading_btag_id) {
        second_jet.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
      }
    }

    // Make preselection cuts.
    if (!found_tau_plus || !found_tau_minus) continue;
    if (!found_btag) continue;
    if (btag_jet.Pt() < 30.) continue;
    if (second_jet.Pt() < 30.) continue;
    if (tau_plus.Pt() < 70.) continue;
    if (tau_minus.Pt() < 70.) continue;
    if (abs(btag_jet.Eta()) > 2.4 || abs(second_jet.Eta()) > 2.4) continue;

    accepted_events++;

    // Now that we have ID'd our 4 particles, calculate kinematics.
    // p = +, m = -
    // b = btag jet, j = other leading jet (btag or non-btag)
    TLorentzVector tau_p_b = tau_plus + btag_jet;
    TLorentzVector tau_m_j = tau_minus + second_jet;
    TLorentzVector tau_p_j = tau_plus + second_jet;
    TLorentzVector tau_m_b = tau_minus + btag_jet;
    TLorentzVector ditau = tau_plus + tau_minus;
    TLorentzVector dijet = btag_jet + second_jet;
    TLorentzVector met_P4;
    std::pair <TLorentzVector, TLorentzVector> b_tau_pair; // Matched pair.

    // Find the right tau-jet pairing by taking the permutation of tau-jet
    // pairs (where jet = b or non-b) with the smallest mass difference, then
    // picking the highest pair mass of that permutation.
    Double_t choice_pair_mass;
    if (abs(tau_p_b.M() - tau_m_j.M()) < abs(tau_p_j.M() - tau_m_b.M())) {
      b_tau_pair = std::make_pair(btag_jet, tau_plus);
      choice_pair_mass = std::max(tau_p_b.M(), tau_m_j.M());
    } else {
      b_tau_pair = std::make_pair(btag_jet, tau_minus);
      choice_pair_mass = std::max(tau_p_j.M(), tau_m_b.M());
    }
    hist_pair_mass->Fill(choice_pair_mass);

    // Calculate deltaR(j,j), deltaR(tau,tau), deltaR(j,tau)_matched, deltaR(tau,MET)
    Double_t deltaR_jets = btag_jet.DeltaR(second_jet);
    Double_t deltaR_taus = tau_plus.DeltaR(tau_minus);
    Double_t deltaR_tau_jet = b_tau_pair.first.DeltaR(b_tau_pair.second);
    Double_t deltaR_tau_met = ditau.DeltaR(met_P4);
    Double_t deltaR_tau_met_max = std::max(tau_plus.DeltaR(met_P4), tau_minus.DeltaR(met_P4));
    Double_t deltaR_topology = deltaR_tau_met_max - deltaR_taus;
    hist_deltaR_jets->Fill(deltaR_jets);
    hist_deltaR_taus->Fill(deltaR_taus);
    hist_deltaR_tau_jet->Fill(deltaR_tau_jet);
    hist_deltaR_tau_met->Fill(deltaR_tau_met);
    hist_deltaR_topology->Fill(deltaR_topology);

    //////////// Begin MT2 related calculations.

    // Calculate and fill MET for the event.
    met = (MissingET*) branch_met->At(0);
    hist_MET->Fill(met->MET / ditau.M());

    // Calculate MT2 (http://www.hep.phy.cam.ac.uk/~lester/mt2/).
    met_P4.SetPtEtaPhiE(met->MET, 0, met->Phi, met->MET);
    Double_t mVisA = tau_plus.M(); // Mass of visible object on side A.
    Double_t pxA = tau_plus.Px(); // x momentum of visible object on side A.
    Double_t pyA = tau_plus.Py(); // y momentum of visible object on side A.

    Double_t mVisB = tau_minus.M(); // Mass of visible object on side B.
    Double_t pxB = tau_minus.Px(); // x momentum of visible object on side B.
    Double_t pyB = tau_minus.Py(); // y momentum of visible object on side B.

    Double_t pxMiss = met_P4.Px(); // x component of missing transverse momentum.
    Double_t pyMiss = met_P4.Py(); // y component of missing transverse momentum.

    Double_t chiA = 0.0; // Hypothesised mass of invisible on side A.
    Double_t chiB = 0.0; // Hypothesised mass of invisible on side B.

    Double_t desiredPrecisionOnMt2 = 0; // Algo aims for machine precision.

    Double_t MT2 = asymm_mt2_lester_bisect::get_mT2(
                   mVisA, pxA, pyA,
                   mVisB, pxB, pyB,
                   pxMiss, pyMiss,
                   chiA, chiB,
                   desiredPrecisionOnMt2);
    hist_MT2->Fill(MT2);

    // UNBOOSTED system.
    // First calculate recoil in the transverse plane.
    TLorentzVector jet_recoil = btag_jet + second_jet;
    TLorentzVector unboosted_tau_plus = tau_plus + jet_recoil;
    TLorentzVector unboosted_tau_minus = tau_minus + jet_recoil;
    TLorentzVector unboosted_met = met_P4 + jet_recoil;

    // Calculate UNBOOSTED MT2.
    Double_t unboost_pxA = unboosted_tau_plus.Px();
    Double_t unboost_pyA = unboosted_tau_plus.Py();
    Double_t unboost_pxB = unboosted_tau_minus.Px();
    Double_t unboost_pyB = unboosted_tau_minus.Py();
    Double_t unboost_pxMiss = unboosted_met.Px();
    Double_t unboost_pyMiss = unboosted_met.Py();

    Double_t unboosted_MT2 = asymm_mt2_lester_bisect::get_mT2(
                             mVisA, unboost_pxA, unboost_pyA,
                             mVisB, unboost_pxB, unboost_pyB,
                             unboost_pxMiss, unboost_pyMiss, chiA,
                             chiB, desiredPrecisionOnMt2);
    hist_unboosted_MT2->Fill(unboosted_MT2 * deltaR_taus);

    hist_denis->Fill((MT2 - unboosted_MT2) / ditau.M());

    // Calculate and fill topology plots.

    // Calculate max{dPhi(MET, mu)} - dPhi(mu,mu)
    // positive --> Topo 1
    // negative --> Topo 2 (trivial zero)
    Double_t max_dphi = std::max(abs(tau_plus.DeltaPhi(met_P4)),
                                 abs(tau_minus.DeltaPhi(met_P4)));
    Double_t dphi_taus = abs(tau_plus.DeltaPhi(tau_minus));
    hist2d_topology->Fill(MT2, max_dphi - dphi_taus);
    hist_topology->Fill(max_dphi - dphi_taus);

    // Calculate unboosted topology.
    Double_t unboost_max_dphi = std::max(abs(unboosted_tau_plus.DeltaPhi(unboosted_met)),
                                         abs(unboosted_tau_minus.DeltaPhi(unboosted_met)));
    Double_t unboost_dphi_taus = abs(unboosted_tau_plus.DeltaPhi(unboosted_tau_minus));
    Double_t tau_met_phi = abs(ditau.DeltaPhi(met_P4));
    hist_j_topology->Fill(tau_met_phi - dphi_taus);

    //////////// End MT2 calculations.

    // Calculate and fill HT - LT.
    Double_t H_T = abs(btag_jet.Pt()) + abs(second_jet.Pt());
    Double_t L_T = abs(tau_plus.Pt()) + abs(tau_minus.Pt());
    hist_HT_LT->Fill(H_T - L_T);

    // Fill ditau pair mass.
    hist_ditau_mass->Fill(ditau.M());

    // Fill Pt spectrums.
    hist_pt_btag->Fill(btag_jet.Pt());
    hist_pt_jet->Fill(second_jet.Pt());
    hist_pt_tau_p->Fill(tau_plus.Pt());
    hist_pt_tau_m->Fill(tau_minus.Pt());

    // Fill Ditau pt.
    hist_pt_ditau->Fill(ditau.Pt());

    // Fill energies.
    hist_e_btag->Fill(btag_jet.E());
    hist_e_jet->Fill(second_jet.E());
    hist_e_tau_p->Fill(tau_plus.E());
    hist_e_tau_m->Fill(tau_minus.E());

    // Fill E ratio.
    hist_e_ratio->Fill(tau_plus.E() / (tau_plus.E() + tau_minus.E()));

    // Fill 2D histograms.
    hist2d_e_tau->Fill(tau_plus.E(), tau_minus.E());
    hist2d_e_pt_tau->Fill(tau_plus.E(), tau_plus.Pt());
    hist2d_pt_tau->Fill(tau_plus.Pt(), tau_minus.Pt());
    hist2d_e_pt_btag->Fill(MT2, tau_plus.Pt());
    hist2d_pt_btag_tau->Fill(b_tau_pair.first.Pt(), b_tau_pair.second.Pt());
    hist2d_eta_btag_tau->Fill(b_tau_pair.first.Eta(), b_tau_pair.second.Eta());

  } // End event loop.

  printf("Out of 40k events %d were accepted. \n", accepted_events);

  // Draw histograms and save them in a .root format.
  hist_pair_mass->Scale(1/hist_pair_mass->Integral());
  hist_MET->Scale(1/hist_MET->Integral());
  hist_HT_LT->Scale(1/hist_HT_LT->Integral());
  hist_ditau_mass->Scale(1/hist_ditau_mass->Integral());
  hist_MT2->Scale(1/hist_MT2->Integral());
  hist_unboosted_MT2->Scale(1/hist_unboosted_MT2->Integral());
  hist_denis->Scale(1/hist_denis->Integral());
  hist_topology->Scale(1/hist_topology->Integral());
  hist_j_topology->Scale(1/hist_j_topology->Integral());
  hist_deltaR_jets->Scale(1/hist_deltaR_jets->Integral());
  hist_deltaR_taus->Scale(1/hist_deltaR_taus->Integral());
  hist_deltaR_tau_jet->Scale(1/hist_deltaR_tau_jet->Integral());
  hist_deltaR_tau_met->Scale(1/hist_deltaR_tau_met->Integral());
  hist_deltaR_topology->Scale(1/hist_deltaR_topology->Integral());
  hist_pt_btag->Scale(1/hist_pt_btag->Integral());
  hist_pt_jet->Scale(1/hist_pt_jet->Integral());
  hist_pt_tau_p->Scale(1/hist_pt_tau_p->Integral());
  hist_pt_tau_m->Scale(1/hist_pt_tau_m->Integral());
  hist_pt_ditau->Scale(1/hist_pt_ditau->Integral());
  hist_e_tau_p->Scale(1/hist_e_tau_p->Integral());
  hist_e_tau_m->Scale(1/hist_e_tau_m->Integral());
  hist_e_btag->Scale(1/hist_e_btag->Integral());
  hist_e_jet->Scale(1/hist_e_jet->Integral());
  hist_e_ratio->Scale(1/hist_e_ratio->Integral());

  // Draw 2D hists.
  TCanvas *c1 = new TCanvas();
  hist2d_e_tau->Draw("COL2Z");

  TCanvas *c2 = new TCanvas();
  hist2d_e_pt_tau->Draw("COL2Z");

  TCanvas *c3 = new TCanvas();
  hist2d_pt_tau->Draw("COL2Z");

  TCanvas *c4 = new TCanvas();
  hist2d_e_pt_btag->Draw("COL2Z");

  TCanvas *c5 = new TCanvas();
  hist2d_pt_btag_tau->Draw("COL2Z");

  TCanvas *c6 = new TCanvas();
  hist2d_eta_btag_tau->Draw("COL2Z");

  TCanvas *c7 = new TCanvas();
  hist2d_topology->Draw("COL2Z");


  TFile *f1 = new TFile("hist_taujet_mass.root", "UPDATE");
  hist_pair_mass->Write(sample_desc, TObject::kOverwrite);

  TFile *f2 = new TFile("hist_MET.root", "UPDATE");
  hist_MET->Write(sample_desc, TObject::kOverwrite);

  TFile *f3 = new TFile("hist_HT_LT.root", "UPDATE");
  hist_HT_LT->Write(sample_desc, TObject::kOverwrite);

  TFile *f4 = new TFile("hist_ditau_mass.root", "UPDATE");
  hist_ditau_mass->Write(sample_desc, TObject::kOverwrite);

  TFile *f5 = new TFile("hist_MT2.root" ,"UPDATE");
  hist_MT2->Write(sample_desc, TObject::kOverwrite);

  TFile *f6 = new TFile("hist_pt_btag.root", "UPDATE");
  hist_pt_btag->Write(sample_desc, TObject::kOverwrite);

  TFile *f7 = new TFile("hist_pt_jet.root", "UPDATE");
  hist_pt_jet->Write(sample_desc, TObject::kOverwrite);

  TFile *f8 = new TFile("hist_pt_taup.root", "UPDATE");
  hist_pt_tau_p->Write(sample_desc, TObject::kOverwrite);

  TFile *f9 = new TFile("hist_pt_taum.root", "UPDATE");
  hist_pt_tau_m->Write(sample_desc, TObject::kOverwrite);

  TFile *f10 = new TFile("hist_pt_ditau.root", "UPDATE");
  hist_pt_ditau->Write(sample_desc, TObject::kOverwrite);

  TFile *f11 = new TFile("hist_e_btag.root", "UPDATE");
  hist_e_btag->Write(sample_desc, TObject::kOverwrite);

  TFile *f12 = new TFile("hist_e_jet.root", "UPDATE");
  hist_e_jet->Write(sample_desc, TObject::kOverwrite);

  TFile *f13 = new TFile("hist_e_taup.root", "UPDATE");
  hist_e_tau_p->Write(sample_desc, TObject::kOverwrite);

  TFile *f14 = new TFile("hist_e_taum.root", "UPDATE");
  hist_e_tau_m->Write(sample_desc, TObject::kOverwrite);

  TFile *f15 = new TFile("hist2d_tautau.root", "UPDATE");
  hist2d_e_tau->Write(sample_desc, TObject::kOverwrite);

  TFile *f16 = new TFile("hist_e_ratio.root", "UPDATE");
  hist_e_ratio->Write(sample_desc, TObject::kOverwrite);

  TFile *f17 = new TFile("hist_unboost_MT2.root", "UPDATE");
  hist_unboosted_MT2->Write(sample_desc, TObject::kOverwrite);

  TFile *f18 = new TFile("hist_topology.root", "UPDATE");
  hist_topology->Write(sample_desc, TObject::kOverwrite);

  TFile *f19 = new TFile("hist_denis_ditau.root", "UPDATE");
  hist_denis->Write(sample_desc, TObject::kOverwrite);

  TFile *f20 = new TFile("hist_deltaR_jets.root", "UPDATE");
  hist_deltaR_jets->Write(sample_desc, TObject::kOverwrite);

  TFile *f21 = new TFile("hist_deltaR_taus.root", "UPDATE");
  hist_deltaR_taus->Write(sample_desc, TObject::kOverwrite);

  TFile *f22 = new TFile("hist_deltaR_tau_jet.root", "UPDATE");
  hist_deltaR_tau_jet->Write(sample_desc, TObject::kOverwrite);

  TFile *f23 = new TFile("hist_deltaR_tau_met.root", "UPDATE");
  hist_deltaR_tau_met->Write(sample_desc, TObject::kOverwrite);

  TFile *f24 = new TFile("hist_deltaR_topology.root","UPDATE");
  hist_deltaR_topology->Write(sample_desc, TObject::kOverwrite);

  TFile *f25 = new TFile("hist_j_topology.root", "UPDATE");
  hist_j_topology->Write(sample_desc, TObject::kOverwrite);
}
