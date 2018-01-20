// An analyzer to produce key selection criteria plots for dimuon processes.
// Run with:
// gSystem->Load("<path_to_delphes>/Delphes/libDelphes.so");

#ifdef __CLING__
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif

#include <algorithm>
#include <iostream>
#include <string>
#include <lester_mt2_bisect.h>

Double_t kMuonMass = 0.105658389;  // GeV
Int_t count_pre = 0;
Int_t count_top_bound = 0;
Int_t count_met = 0;
Int_t count_ht_lt = 0;
Int_t count_UnbMT2 = 0;
Int_t event_count = 0;

void DiMuonAnalyzer(const char *sample_desc, int nbins) {
  gSystem->Load("libDelphes.so");

  // Create a chain of root trees.
  TChain chain("Delphes");

  // ttbar samples.
  chain.Add("/fdata/scratch/mexanick/ttbar/tt_b.root");
  chain.Add("/fdata/scratch/mexanick/ttbar/tt_c.root");
  chain.Add("/fdata/scratch/mexanick/ttbar/tt_g.root");
  chain.Add("/fdata/scratch/mexanick/ttbar/tt_j.root");
  chain.Add("/fdata/scratch/mexanick/ttbar/tt_k.root");
  chain.Add("/fdata/scratch/mexanick/ttbar/tt_l.root");

  // Z' (500 GeV) samples.
/*  chain.Add("/fdata/scratch/mexanick/z_prime_500GeV_cs1/Events/run_05/tag_2_delphes_events.root");
  chain.Add("/fdata/scratch/mexanick/z_prime_500GeV_cs1/Events/run_06/tag_2_delphes_events.root");
  chain.Add("/fdata/scratch/mexanick/z_prime_500GeV_cs1/Events/run_07/tag_2_delphes_events.root");
  chain.Add("/fdata/scratch/mexanick/z_prime_500GeV_cs1/Events/run_08/tag_1_delphes_events.root");
  chain.Add("/fdata/scratch/mexanick/z_prime_500GeV_cs2/Events/run_05/tag_2_delphes_events.root");
  chain.Add("/fdata/scratch/mexanick/z_prime_500GeV_cs2/Events/run_06/tag_2_delphes_events.root");
  chain.Add("/fdata/scratch/mexanick/z_prime_500GeV_cs2/Events/run_07/tag_2_delphes_events.root");
  chain.Add("/fdata/scratch/mexanick/z_prime_500GeV_cs2/Events/run_08/tag_1_delphes_events.root");
  chain.Add("/fdata/scratch/mexanick/z_prime_500GeV_cs5/Events/run_01/tag_6_delphes_events.root");
  chain.Add("/fdata/scratch/mexanick/z_prime_500GeV_cs5/Events/run_02/tag_5_delphes_events.root");
  chain.Add("/fdata/scratch/mexanick/z_prime_500GeV_cs5/Events/run_03/tag_5_delphes_events.root");
*/
  // Create object of class ExRootTreeReader.
  ExRootTreeReader *tree_reader = new ExRootTreeReader(&chain);
  Long64_t number_of_entries = tree_reader->GetEntries();

  // Get pointers to branches used in this analysis.
  TClonesArray *branch_jet = tree_reader->UseBranch("Jet");
  TClonesArray *branch_met = tree_reader->UseBranch("MissingET");
  TClonesArray *branch_muon = tree_reader->UseBranch("Muon");

  // Book histograms.
  TH1F *hist_pair_mass = new TH1F("pair_mass", "Max{M(j,mu)}", nbins, 0., 300.);
  TH1F *hist_MET = new TH1F("met", "Normalized Missing ET", nbins, 0., 1.);
  TH1F *hist_HT_LT = new TH1F("HT_LT", "HT - LT", nbins, -500., 500.);
  TH1F *hist_MT2 = new TH1F("MT2", "MT2", nbins, 0., 1.);
  TH1F *hist_unboosted_MT2 = new TH1F("unboost_mt2", "", nbins, 0., 40.);

  // Event loop.
  for(Int_t entry = 0; entry < number_of_entries; ++entry) {

    tree_reader->ReadEntry(entry);

    event_count++;
    if (event_count % 10000 == 0) {
      printf("Processed %d / %lld events \n", event_count, number_of_entries);
      printf("Preselection: %d / %d events \n", count_pre, event_count);
      printf("Top Mass Bound: %d / %d events \n", count_top_bound, event_count);
      printf("Normalized MET: %d / %d events \n", count_met, event_count);
      printf("HT - LT: %d / %d events \n", count_ht_lt, event_count);
      printf("Unboosted MT2: %d / %d events \n", count_UnbMT2, event_count);
    }

    // Declare physics objects.
    Muon *muon;
    Jet *jet;
    MissingET *met;
    TLorentzVector btag_jet;
    TLorentzVector second_jet;
    TLorentzVector mu_plus;
    TLorentzVector mu_minus;

    Int_t jet_size = branch_jet->GetEntries();
    Int_t met_size = branch_met->GetEntries();
    Int_t muon_size = branch_muon->GetEntries();

    if (jet_size < 2) continue;
    bool found_btag = false;
    bool found_jet = false;
    bool found_mu_plus = false;
    bool found_mu_minus = false;

    // Find the leading OS Muon pair.
    for (Int_t ii = 0; ii < muon_size; ii++) {
      muon = (Muon*) branch_muon->At(ii);
      // Looking for highest-pT mu+.
      if (muon->Charge == 1 && muon->PT > mu_plus.Pt()) {
        mu_plus.SetPtEtaPhiM(muon->PT, muon->Eta, muon->Phi, kMuonMass);
        found_mu_plus = true;
      }
      // Looking for highest-pT mu-.
      if (muon->Charge == -1 && muon->PT > mu_minus.Pt()) {
        mu_minus.SetPtEtaPhiM(muon->PT, muon->Eta, muon->Phi, kMuonMass);
        found_mu_minus = true;
      }
    }

    // Loop over jets and find the highest-pT BTag jet.
    Int_t leading_btag_id = -1;
    for (Int_t ii = 0; ii < jet_size; ii++) {
      jet = (Jet*) branch_jet->At(ii);
      if (jet->BTag == 1 && jet->PT > btag_jet.Pt()) {
        btag_jet.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
        leading_btag_id = ii;
        found_btag = true;
      }
    }

    // Find the other highest-pT jet with ID different than the leading BTag.
    for (Int_t ii = 0; ii < jet_size; ii++) {
      jet = (Jet*) branch_jet->At(ii);
      // Must not be a Tau jet with ID different than the leading b-tag.
      if (jet->TauTag == 0 && jet->PT > second_jet.Pt() && ii != leading_btag_id) {
        second_jet.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
        found_jet = true;
      }
    }
    // Place preselection cuts.
    if (!found_btag || !found_jet) continue;
    if (!found_mu_plus || !found_mu_minus) continue;
    if (btag_jet.Pt() < 30.) continue;
    if (second_jet.Pt() < 30.) continue;
    if (abs(btag_jet.Eta()) > 2.4 || abs(second_jet.Eta()) > 2.4) continue;

    // Now that we have ID'd our 4 particles, construct more 4-vectors.
    // p = +, m = -
    // b = btag jet, j = other leading jet (btag or non-btag)
    TLorentzVector mu_p_b = mu_plus + btag_jet;
    TLorentzVector mu_m_j = mu_minus + second_jet;
    TLorentzVector mu_p_j = mu_plus + second_jet;
    TLorentzVector mu_m_b = mu_minus + btag_jet;
    TLorentzVector dimuon = mu_plus + mu_minus;
    TLorentzVector met_P4;

    // Fill histograms.

    // (1) Calculate and fill muon/jet pair mass.
    Double_t top_mass_bound;
    // If (mu+, b) - (mu-, j) has the smallest mass difference...
    if (abs(mu_p_b.M() - mu_m_j.M()) < abs(mu_p_j.M() - mu_m_b.M())) {
      top_mass_bound = std::max(mu_p_b.M(), mu_m_j.M()); // ... take the max pair.
    } else { // otherwise, take the max of the other set of pairs.
      top_mass_bound = std::max(mu_p_j.M(), mu_m_b.M());
    }
    hist_pair_mass->Fill(top_mass_bound);

    // (2) Calculate and fill MET for the event.
    met = (MissingET*) branch_met->At(0);
    Double_t normalized_met = (met->MET) / (dimuon.M());
    hist_MET->Fill(normalized_met);

    // (3) Calculate and fill HT - LT.
    Double_t H_T = abs(btag_jet.Pt()) + abs(second_jet.Pt());
    Double_t L_T = abs(mu_plus.Pt()) + abs(mu_minus.Pt());
    hist_HT_LT->Fill(H_T - L_T);

    // (4) Calculate MT2 (http://www.hep.phy.cam.ac.uk/~lester/mt2/).
    met_P4.SetPtEtaPhiM(met->MET, met->Eta, met->Phi, 0);
    double mVisA = kMuonMass; // Mass of visible object on side A.
    double pxA = mu_plus.Px(); // x momentum of visible object on side A.
    double pyA = mu_plus.Py(); // y momentum of visible object on side A.

    double mVisB = kMuonMass; // Mass of visible object on side B.
    double pxB = mu_minus.Px(); // x momentum of visible object on side B.
    double pyB = mu_minus.Py(); // y momentum of visible object on side B.

    double pxMiss = met_P4.Px(); // x component of missing transverse momentum.
    double pyMiss = met_P4.Py(); // y component of missing transverse momentum.

    double chiA = 0.0; // Hypothesised mass of invisible on side A.
    double chiB = 0.0; // Hypothesised mass of invisible on side B.

    double desiredPrecisionOnMt2 = 0.0; // Algo aims for machine precision.

    double MT2 = asymm_mt2_lester_bisect::get_mT2(mVisA, pxA, pyA,
                                                mVisB, pxB, pyB,
                                                pxMiss, pyMiss,
                                                chiA, chiB,
                                                desiredPrecisionOnMt2);
    hist_MT2->Fill(MT2);

    // UNBOOSTED system.
    // First calculate recoil in the transverse plane.
    TLorentzVector jet_recoil = btag_jet + second_jet;
    TLorentzVector unboosted_mu_plus = jet_recoil + mu_plus;
    TLorentzVector unboosted_mu_minus = jet_recoil + mu_minus;
    TLorentzVector unboosted_met = jet_recoil + met_P4;

    // (5) Calculate UNBOOSTED MT2.
    double unboost_pxA = unboosted_mu_plus.Px();
    double unboost_pyA = unboosted_mu_plus.Py();

    double unboost_pxB = unboosted_mu_minus.Px();
    double unboost_pyB = unboosted_mu_minus.Py();
    double unboost_pxMiss = unboosted_met.Px();
    double unboost_pyMiss = unboosted_met.Py();

    double unboosted_MT2 = asymm_mt2_lester_bisect::get_mT2(
                             mVisA, unboost_pxA, unboost_pyA,
                             mVisB, unboost_pxB, unboost_pyB,
                             unboost_pxMiss, unboost_pyMiss, chiA, chiB,
                             desiredPrecisionOnMt2);
    hist_unboosted_MT2->Fill(unboosted_MT2);


    // Calculate things for acceptance.
    count_pre++;
    if (top_mass_bound > 170) {
      count_top_bound++;
      if (normalized_met < 0.2) {
        count_met++;
        if (H_T - L_T < 0) {
          count_ht_lt++;
        }
      }
    }

  } // End event loop.

  printf("Preselection: %d / %d events \n", count_pre, event_count);
  printf("Top Mass Bound: %d / %d events \n", count_top_bound, event_count);
  printf("Normalized MET: %d / %d events \n", count_met, event_count);
  printf("HT - LT: %d / %d events \n", count_ht_lt, event_count);
  printf("Unboosted MT2: %d / %d events \n", count_UnbMT2, event_count);



  // Draw histograms and save them in a .root format.
  hist_pair_mass->Scale(1/hist_pair_mass->Integral());
  hist_MET->Scale(1/hist_MET->Integral());
  hist_HT_LT->Scale(1/hist_HT_LT->Integral());
  hist_MT2->Scale(1/hist_MT2->Integral());
  hist_unboosted_MT2->Scale(1/hist_unboosted_MT2->Integral());

  TFile *f1 = new TFile("hist_dimuon_mujet_mass.root", "UPDATE");
  hist_pair_mass->Write(sample_desc, TObject::kOverwrite);

  TFile *f2 = new TFile("hist_dimuon_MET.root", "UPDATE");
  hist_MET->Write(sample_desc, TObject::kOverwrite);

  TFile *f3 = new TFile("hist_dimuon_HT_LT.root", "UPDATE");
  hist_HT_LT->Write(sample_desc, TObject::kOverwrite);

  TFile *f4 = new TFile("hist_dimuon_MT2.root" ,"UPDATE");
  hist_MT2->Write(sample_desc, TObject::kOverwrite);

  TFile *f5 = new TFile("hist_dimuon_unboost_MT2.root", "UPDATE");
  hist_unboosted_MT2->Write(sample_desc, TObject::kOverwrite);

}
