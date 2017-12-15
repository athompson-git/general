// An analyzer to produce key selection criteria plots for Z' -> ditau processes.
// USAGE:
// gSystem->Load("<path_to_delphes>/Delphes/libDelphes.so");
// .x ditau_analyzer.C("input_file.root", "saved_histogram_name.root", nbins);

#ifdef __CLING__
//R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif

void ditau_analyzer(const char *file_name, const char *sample_desc, int nbins) {
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
  TH1F *hist_pair_mass = new TH1F("pair_mass", "Max{M(j,tau)}", nbins, 0., 300.);
  hist_pair_mass->GetXaxis()->SetTitle("M(tau,j)");
  hist_pair_mass->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_MET = new TH1F("met", "Normalized Missing ET", nbins, 0., 1.);
  hist_MET->GetXaxis()->SetTitle("MET/M(tautau)");
  hist_MET->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_HT_LT = new TH1F("HT_LT", "HT - LT", nbins, -500., 500.);
  hist_HT_LT->GetXaxis()->SetTitle("HT - LT");
  hist_HT_LT->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_eta = new TH1F("eta", "BTag Eta", nbins, -6., 6.);
  hist_eta->GetXaxis()->SetTitle("Eta");
  hist_eta->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_pt = new TH1F("pt", "BTag Pt", nbins, 0., 300.);
  hist_pt->GetXaxis()->SetTitle("Pt (GeV)");
  hist_pt->GetYaxis()->SetTitle("a.u.");

  // Event loop.
  for(Int_t entry = 0; entry < number_of_entries; ++entry) {

    tree_reader->ReadEntry(entry);

    // Declare physics objects.
    Muon *muon;
    Jet *jet;
    MissingET *met;

    // Declare running indices.
    Int_t jet_size = branch_jet->GetEntries();
    Int_t met_size = branch_met->GetEntries();

    // Skip events if they do not have 4 or more jets.
    if (jet_size < 4) {
      printf("Event #%d: did not find at least 4 jets, skipping... \n", entry);
      continue;
    }

    // Skip events if they do not have a BTag jet.
    bool found_btag = false;
    // Skip events if they do not have a Tau+ and Tau-.
    bool found_tau_plus = false;
    bool found_tau_minus = false;

    // Declare 4-vectors.
    TLorentzVector btag_jet;
    TLorentzVector second_jet;
    TLorentzVector tau_plus;
    TLorentzVector tau_minus;

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

    if (!found_tau_plus || !found_tau_minus) {
      printf("Event # %d: did not find all Tau's, skipping... \n", entry);
      continue;
    }
    if (!found_btag) {
      printf("Event # %d: did not find a b-tagged jet, skipping... \n", entry);
      continue;
    }

    // Now that we have ID'd our 4 particles, construct more 4-vectors.
    // p = +, m = -
    // b = btag jet, j = other leading jet (btag or non-btag)
    TLorentzVector tau_p_b = tau_plus + btag_jet;
    TLorentzVector tau_m_j = tau_minus + second_jet;
    TLorentzVector tau_p_j = tau_plus + second_jet;
    TLorentzVector tau_m_b = tau_minus + btag_jet;
    TLorentzVector ditau = tau_plus + tau_minus;

    // Fill histograms.

    // Calculate and fill tau-jet pair mass.
    Double_t choice_pair_mass;
    if (abs(tau_p_b.M() - tau_m_j.M()) < abs(tau_p_j.M() - tau_m_b.M())) {
      choice_pair_mass = std::max(tau_p_b.M(), tau_m_j.M());
    } else {
      choice_pair_mass = std::max(tau_p_j.M(), tau_m_b.M());
    }
    hist_pair_mass->Fill(choice_pair_mass);

    // Calculate and fill normalized Missing ET.
    Double_t met_integral = 0.0;
    for (Int_t ii = 0; ii < met_size; ii++) {
      met = (MissingET*) branch_met->At(ii);
      met_integral += met->MET;
    }
    hist_MET->Fill(met_integral / ditau.M());

    // Calculate and fill HT - LT.
    Double_t H_T = abs(btag_jet.Pt()) + abs(second_jet.Pt());
    Double_t L_T = abs(tau_plus.Pt()) + abs(tau_minus.Pt());
    hist_HT_LT->Fill(H_T - L_T);

    // Fill eta spectrum for BTag.
    hist_eta->Fill(btag_jet.Eta());

    // Fill Pt spectrum for the BTag
    hist_pt->Fill(btag_jet.Pt());

  } // End event loop.

  // Draw histograms and save them in a .root format.
  TCanvas *c1 = new TCanvas();
  hist_pair_mass->Scale(1/hist_pair_mass->Integral());
  hist_pair_mass->Draw("HIST");

  TCanvas *c2 = new TCanvas();
  hist_MET->Scale(1/hist_MET->Integral());
  hist_MET->Draw("HIST");

  TCanvas *c3 = new TCanvas();
  hist_HT_LT->Scale(1/hist_HT_LT->Integral());
  hist_HT_LT->Draw("HIST");

  TCanvas *c4 = new TCanvas();
  hist_eta->Scale(1/hist_eta->Integral());
  hist_eta->Draw("HIST");

  TCanvas *c5 = new TCanvas();
  hist_pt->Scale(1/hist_eta->Integral());
  hist_pt->Draw("HIST");

  std::string mass_name = string(sample_desc) + "_pair_mass";
  std::string met_name = string(sample_desc) + "_met";
  std::string HT_LT_name = string(sample_desc) + "_HT_LT";
  std::string eta_name = string(sample_desc) + "_eta";
  std::string pt_name = string(sample_desc) + "_pt";

  TFile *f1 = new TFile("pair_mass.root", "UPDATE");
  hist_pair_mass->Write(mass_name.c_str(), TObject::kOverwrite);

  TFile *f2 = new TFile("hist_MET.root", "UPDATE");
  hist_MET->Write(met_name.c_str(), TObject::kOverwrite);

  TFile *f3 = new TFile("hist_HT_LT.root", "UPDATE");
  hist_HT_LT->Write(HT_LT_name.c_str(), TObject::kOverwrite);

  TFile *f4 = new TFile("hist_eta.root" ,"UPDATE");
  hist_eta->Write(eta_name.c_str(), TObject::kOverwrite);

  TFile *f5 = new TFile("hist_pt.root", "UPDATE");
  hist_pt->Write(pt_name.c_str(), TObject::kOverwrite);

}
