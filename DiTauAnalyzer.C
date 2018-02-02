// An analyzer to produce key selection criteria plots for Z' -> ditau processes.
// arXiv:1707.07016v1
// USAGE:
// gSystem->Load("<path_to_delphes>/Delphes/libDelphes.so");
// .x DiTauAnalyzer.C("input_file.root", "Z' 500 GeV", nbins);

#ifdef __CLING__
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TEfficiency.h>
#include <TCanvas.h>

#include "lester_mt2_bisect.h"
#include <vector>
#include <string>
#include <utility>
#include <iostream>

// Set constants.
double kTauMass = 1.77682; // GeV

// MT2 Calculator.
double mt2(TLorentzVector vis1, TLorentzVector vis2, TLorentzVector miss) {
  double mVisA = vis1.M(); //vis1.M(); // Mass of visible object on side A.
  double pxA = vis1.Px(); // x momentum of visible object on side A.
  double pyA = vis1.Py(); // y momentum of visible object on side A.

  double mVisB = vis1.M(); //vis2.M(); // Mass of visible object on side B.
  double pxB = vis2.Px(); // x momentum of visible object on side B.
  double pyB = vis2.Py(); // y momentum of visible object on side B.

  double pxMiss = miss.Px(); // x component of missing transverse momentum.
  double pyMiss = miss.Py(); // y component of missing transverse momentum.

  double chiA = 0.0; // Hypothesised mass of invisible on side A.
  double chiB = 0.0; // Hypothesised mass of invisible on side B.

  double desiredPrecisionOnMt2 = 0.0; // Algo aims for machine precision.

  asymm_mt2_lester_bisect lester;
  return lester.get_mT2(mVisA, pxA, pyA, mVisB, pxB, pyB,
                        pxMiss, pyMiss, chiA, chiB,
                        desiredPrecisionOnMt2);
}


// Helper function to express a vector V in a new basis defined by two vectors (b1, b2).
TLorentzVector change_basis(TLorentzVector V, TLorentzVector b1,
                            TLorentzVector b2) {
  TVector3 normal = (b1.Vect()).Cross(b2.Vect());

  // Define orthonormal basis vectors (e1, e2, e3) in terms of global coords.
  TVector3 e1 = (1 / (b1.Vect()).Mag()) * (b1.Vect());
  TVector3 e3 = (1 / normal.Mag()) * normal;
  TVector3 e2 = e3.Cross(e1);

  // Define local coords.
  TVector3 w1(1., 0., 0.);
  TVector3 w2(0., 1., 0.);
  TVector3 w3(0., 0., 1.);

  TVector3 V_prime = ((V.Vect()).Dot(e1) * w1)
                   + ((V.Vect()).Dot(e2) * w2)
                   + ((V.Vect()).Dot(e3) * w3);
  TLorentzVector W;
  W.SetPxPyPzE(V_prime.Px(), V_prime.Py(), V_prime.Pz(), V.E());
  return W;
}


// Helper function to project a vector into the plane defined by two vectors.
TLorentzVector plane_projection(TLorentzVector a, TLorentzVector plane_vec1,
                                TLorentzVector plane_vec2) {
  TVector3 normal = (plane_vec1.Vect()).Cross(plane_vec2.Vect());
  Double_t scalar_proj_on_normal = (a.Vect()).Dot(normal) / normal.Mag2();
  TVector3 a_prime = a.Vect() - scalar_proj_on_normal * normal;
  TLorentzVector P;
  P.SetPxPyPzE(a_prime.Px(), a_prime.Py(), a_prime.Pz(), a.E());
  return P;
}



void DiTauAnalyzer(const char *file_name, const char *sample_desc, int nbins,
                   bool apply_cuts = false) {
  gSystem->Load("libDelphes.so");

  // Create a chain of root trees.
  TChain chain("Delphes");
  chain.Add(file_name);

  // Create object of class ExRootTreeReader.
  ExRootTreeReader *tree_reader = new ExRootTreeReader(&chain);
  Long64_t number_of_entries = tree_reader->GetEntries();
  Int_t accepted_events = 0;

  // Get pointers to branches used in this analysis.
  TClonesArray *branch_jet = tree_reader->UseBranch("Jet");
  TClonesArray *branch_met = tree_reader->UseBranch("MissingET");

  // Book histograms.
  TH1F *hist_pair_mass = new TH1F("pair_mass", "Max{M(tau,j)}", nbins, 0., 500.);
  TH1F *hist_MET = new TH1F("met", "Normalized Missing ET", nbins, 0., 2.);
  TH1F *hist_MT2 = new TH1F("tmass", "MT2 (Ditau + MET)", nbins, 0., 150.);
  TH1F *hist_HT_LT = new TH1F("HT_LT", "HT - LT", nbins, -1000., 1000.);
  TH1F *hist_ditau_mass = new TH1F("ditau_mass", "M(tau+,tau-)", nbins, 0., 500.);
  TH1F *hist_unboosted_MT2 = new TH1F("unboost_mt2", "", nbins, 0., 150.);
  TH1F *hist_primed_htlt = new TH1F("hist_unbhtlt", "", nbins, -1600, 600);
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
  TH1F *hist_topology = new TH1F("hist_topo", "", nbins, -3.14, 3.14);
  TH1F *hist_j_topology = new TH1F("hist_unboosted_topo", "", nbins, -3.14, 3.14);
  TH1F *hist_deltaR_jets = new TH1F("hist_deltaR_1", "", nbins, 0., 6.);
  TH1F *hist_deltaR_taus = new TH1F("hist_deltaR_2", "", nbins, 0., 3.14);
  TH1F *hist_deltaR_tau_jet = new TH1F("hist_deltaR_3", "", nbins, 0., 6.);
  TH1F *hist_deltaR_tau_met = new TH1F("hist_deltaR_4", "", nbins, 0., 6.);
  TH1F *hist_deltaR_topology = new TH1F("hist_deltaR_topology", "", nbins, -6., 6.);
  TH1F *hist_deltaPhi_tau_met = new TH1F("hist_deltaPhi_1", "", nbins, -6., 6.);

  // Book 2D histograms for exploratory analysis.
  TH2F *hist2d_topology = new TH2F("topo", "Topology", nbins,
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

    if (entry % 10000 == 0) printf("On event %d / %lld \n", entry, number_of_entries);
    tree_reader->ReadEntry(entry);
    vector<int> b_pos;
    vector<int> tau_pos;
    int nonb_pos = -1;
    float nonb_pt = 0.;
    bool dibottom = false;
    bool bottom_and_jet = false;
    bool os_taus = false;

    // If the event contains at least one jet, push back the vector of b's and
    // record the branch position of the highest-pT non-b jet.
    // Also push back a vector of Tau's.
    for (unsigned j = 0; j < branch_jet->GetEntries(); ++j) {
      Jet *jet = (Jet*) branch_jet->At(j);

      if (jet->PT < 30.) continue;  // Implicitly cut on jets pT < 30 GeV.
      if (jet->TauTag && jet->PT < 70.) continue;  // For taus, cut pT < 70 GeV.

      if (jet->BTag) {
        b_pos.push_back(j);
      } else if (nonb_pt < jet->PT && !(jet->TauTag)) {
        nonb_pos = j;
        nonb_pt = jet->PT;
      } else if (jet->TauTag) {
        tau_pos.push_back(j);
      }
    }

    if (b_pos.size() >= 2) dibottom = true;
    if (b_pos.size() == 1 && nonb_pt > 0.) bottom_and_jet = true;

    Jet *tau1, *tau2;
    // If the event contains at least two muons, check if they are OS and
    // keep them.
    if (tau_pos.size() > 1) {
      tau1 = (Jet *) branch_jet->At(tau_pos[0]);
      tau2 = (Jet *) branch_jet->At(tau_pos[1]);

      if((tau1->Charge) != (tau2->Charge)) os_taus = true;
    }

    // Make preselection cut.
    if (!os_taus) continue;
    if (!(dibottom || bottom_and_jet)) continue;
    accepted_events++;

    // Make jet assignments.
    Jet *b1, *b2;
    if (dibottom) {
      b1 = (Jet *) branch_jet->At(b_pos[0]);
      b2 = (Jet *) branch_jet->At(b_pos[1]);
    } else if (bottom_and_jet) {
      b1 = (Jet *) branch_jet->At(b_pos[0]);
      b2 = (Jet *) branch_jet->At(nonb_pos);
    }

    // End of selection. ////////////////////////////////////////

    bool PassMissingET = false;
    bool PassHTLT = false;
    bool PassMuJetMass = false;
    bool PassUnboostedMT2 = false;

    // Declare physics objects.
    MissingET *MPT = (MissingET *) branch_met->At(0);
    TLorentzVector tau1_p4 = tau1->P4();
    TLorentzVector tau2_p4 = tau2->P4();
    TLorentzVector met_p4 = MPT->P4();
    TLorentzVector b1_p4 = b1->P4();
    TLorentzVector b2_p4 = b2->P4();


    // Now that we have ID'd our 4 particles, calculate kinematics.
    // p = +, m = -
    // b = btag jet, j = other leading jet (btag or non-btag)
    TLorentzVector pm_11 = tau1_p4 + b1_p4;
    TLorentzVector pm_22 = tau2_p4 + b2_p4;
    TLorentzVector pm_12 = tau1_p4 + b2_p4;
    TLorentzVector pm_21 = tau2_p4 + b1_p4;
    TLorentzVector ditau = tau1_p4 + tau2_p4;
    TLorentzVector dijet = b1_p4 + b2_p4;
    std::pair <TLorentzVector, TLorentzVector> b_tau_pair; // Matched pair.


    // (1) max{M(tau,j)}
    // Find the right tau-jet pairing by taking the permutation of tau-jet
    // pairs (where jet = b or non-b) with the smallest mass difference, then
    // picking the highest pair mass of that permutation.
    Double_t choice_pair_mass;
    if (abs(pm_11.M() - pm_22.M()) < abs(pm_12.M() - pm_21.M())) {
      b_tau_pair = std::make_pair(b1_p4, tau1_p4);
      choice_pair_mass = std::max(pm_11.M(), pm_22.M());
    } else {
      b_tau_pair = std::make_pair(b1_p4, tau2_p4);
      choice_pair_mass = std::max(pm_12.M(), pm_21.M());
    }
    hist_pair_mass->Fill(choice_pair_mass);

    // (2) MET / M(tau,tau)
    hist_MET->Fill(MPT->MET / ditau.M());

    // (3) HT - LT
    Double_t H_T = b1_p4.Pt() + b2_p4.Pt();
    Double_t L_T = tau1_p4.Pt() + tau2_p4.Pt();
    hist_HT_LT->Fill(H_T - L_T);

    //////////// Begin MT2 related calculations.

    // Calculate MT2 (http://www.hep.phy.cam.ac.uk/~lester/mt2/).
    double MT2 = mt2(tau1_p4, tau2_p4, met_p4);
    hist_MT2->Fill(MT2);

    // UNBOOSTED system.
    TLorentzVector jet_recoil = plane_projection(dijet, tau1_p4, tau2_p4);
    TLorentzVector proj_met = plane_projection(met_p4, tau1_p4, tau2_p4);

    // Re-represent Taus, recoil, and met in an orthnormal basis defined by the ditau plane.
    TLorentzVector tp_prime = change_basis(tau1_p4, tau1_p4, tau2_p4);
    TLorentzVector tm_prime = change_basis(tau2_p4, tau1_p4, tau2_p4);
    TLorentzVector met_prime = change_basis(proj_met, tau1_p4, tau2_p4);
    TLorentzVector recoil_prime = change_basis(jet_recoil, tau1_p4, tau2_p4);
    TLorentzVector bjet_prime = change_basis(b1_p4, tau1_p4, tau2_p4);
    TLorentzVector jet_prime = change_basis(b2_p4, tau1_p4, tau2_p4);

    TLorentzVector unb_tp_prime = tp_prime + recoil_prime;
    TLorentzVector unb_tm_prime = tm_prime + recoil_prime;
    TLorentzVector unb_met_prime = met_prime + recoil_prime;


    // Calculate MT2 in the ditau plane.
    // double MT2_prime = mt2(unb_tp_prime, unb_tm_prime, unb_met_prime);
    // Regular unboost.
    // double MT2_prime = mt2(tau1_p4 + dijet, tau2_p4 + dijet, met_p4 + dijet);
    // Manual.
    double MT2_prime = asymm_mt2_lester_bisect::get_mT2(tau1_p4.M(),
                        unb_tp_prime.Px(), unb_tp_prime.Py(), tau2_p4.M(),
                        unb_tm_prime.Px(), unb_met_prime.Px(),
                        unb_met_prime.Py(), 0., 0., 0.);
    hist_unboosted_MT2->Fill(MT2_prime);
    hist_primed_htlt->Fill(bjet_prime.Pt() + jet_prime.Pt() - tp_prime.Pt()
                           - tm_prime.Pt());

    //////////// End MT2 calculations.


    // Calculate max{dPhi(MET, mu)} - dPhi(mu,mu)
    // positive --> Topo 1
    // negative --> Topo 2 (trivial zero)
    Double_t max_dphi = std::max(abs(unb_tp_prime.DeltaPhi(unb_met_prime)),
                                 abs(unb_tm_prime.DeltaPhi(unb_met_prime)));
    Double_t dphi_taus = abs(unb_tp_prime.DeltaPhi(unb_tm_prime));
    hist2d_topology->Fill(MT2_prime, max_dphi - dphi_taus);
    hist_topology->Fill(max_dphi - dphi_taus);

    // Fill ditau pair mass.
    hist_ditau_mass->Fill(ditau.M());

    // Fill Pt spectrums.
    hist_pt_btag->Fill(b1_p4.Pt());
    hist_pt_jet->Fill(b2_p4.Pt());
    hist_pt_tau_p->Fill(tau1_p4.Pt());
    hist_pt_tau_m->Fill(tau2_p4.Pt());

    // Fill Ditau pt.
    hist_pt_ditau->Fill(ditau.Pt());

    // Fill energies.
    hist_e_btag->Fill(b1_p4.E());
    hist_e_jet->Fill(b2_p4.E());
    hist_e_tau_p->Fill(tau1_p4.E());
    hist_e_tau_m->Fill(tau2_p4.E());

    // Fill E ratio.
    hist_e_ratio->Fill(tau1_p4.E() / (tau1_p4.E() + tau2_p4.E()));

    // Fill angular quantities.
    Double_t deltaR_jets = b1_p4.DeltaR(b2_p4);
    Double_t deltaR_taus = tau1_p4.DeltaR(tau2_p4);
    Double_t deltaR_tau_jet = b_tau_pair.first.DeltaR(b_tau_pair.second);
    Double_t deltaR_tau_met = ditau.DeltaR(met_p4);
    Double_t deltaR_tau_met_max = std::max(tau1_p4.DeltaR(met_p4), tau2_p4.DeltaR(met_p4));
    Double_t deltaR_topology = deltaR_tau_met_max - deltaR_taus;
    hist_deltaR_jets->Fill(deltaR_jets);
    hist_deltaR_tau_jet->Fill(deltaR_tau_jet);
    hist_deltaR_tau_met->Fill(deltaR_tau_met);
    hist_deltaR_topology->Fill(deltaR_topology);
    hist_deltaR_taus->Fill(tau1_p4.DeltaR(tau2_p4));

    if (tau1_p4.Pt() > tau2_p4.Pt()) {
      hist_deltaPhi_tau_met->Fill(fabs(tau1_p4.DeltaPhi(met_p4)));
    } else {
      hist_deltaPhi_tau_met->Fill(fabs(tau2_p4.DeltaPhi(met_p4)));
    }

    // Fill 2D histograms.
    hist2d_e_tau->Fill(tau1_p4.E(), tau2_p4.E());
    hist2d_e_pt_tau->Fill(tau1_p4.E(), tau1_p4.Pt());
    hist2d_pt_tau->Fill(tau1_p4.Pt(), tau2_p4.Pt());
    hist2d_e_pt_btag->Fill(MT2, tau1_p4.Pt());
    hist2d_pt_btag_tau->Fill(b_tau_pair.first.Pt(), b_tau_pair.second.Pt());
    hist2d_eta_btag_tau->Fill(b_tau_pair.first.Eta(), b_tau_pair.second.Eta());

  } // End event loop.

  printf("%d / %lld accepted \n", accepted_events, number_of_entries);

  // Draw histograms and save them in a .root format.
  hist_pair_mass->Scale(1/hist_pair_mass->Integral());
  hist_MET->Scale(1/hist_MET->Integral());
  hist_HT_LT->Scale(1/hist_HT_LT->Integral());
  hist_ditau_mass->Scale(1/hist_ditau_mass->Integral());
  hist_MT2->Scale(1/hist_MT2->Integral());
  hist_unboosted_MT2->Scale(1/hist_unboosted_MT2->Integral());
  hist_primed_htlt->Scale(1/hist_primed_htlt->Integral());
  hist_topology->Scale(1/hist_topology->Integral());
  hist_j_topology->Scale(1/hist_j_topology->Integral());
  hist_deltaR_jets->Scale(1/hist_deltaR_jets->Integral());
  hist_deltaR_taus->Scale(1/hist_deltaR_taus->Integral());
  hist_deltaR_tau_jet->Scale(1/hist_deltaR_tau_jet->Integral());
  hist_deltaR_tau_met->Scale(1/hist_deltaR_tau_met->Integral());
  hist_deltaR_topology->Scale(1/hist_deltaR_topology->Integral());
  hist_deltaPhi_tau_met->Scale(1/hist_deltaPhi_tau_met->Integral());
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

  TFile *f19 = new TFile("hist_primed_htlt_ditau.root", "UPDATE");
  hist_primed_htlt->Write(sample_desc, TObject::kOverwrite);

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

  TFile *f26 = new TFile("hist_deltaPhi_tau_met.root", "UPDATE");
  hist_deltaPhi_tau_met->Write(sample_desc, TObject::kOverwrite);
}
