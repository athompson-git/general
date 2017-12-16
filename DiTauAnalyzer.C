// An analyzer to produce key selection criteria plots for Z' -> ditau processes.
// USAGE:
// gSystem->Load("<path_to_delphes>/Delphes/libDelphes.so");
// .x DiTauAnalyzer.C("input_file.root", "saved_histogram_name.root", nbins);

#ifdef __CLING__
//R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif

void DiTauAnalyzer(const char *file_name, const char *sample_desc, int nbins) {
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
  TH1F *hist_pair_mass = new TH1F("pair_mass", "Max{M(tau,j)}", nbins, 0., 300.);
  hist_pair_mass->GetXaxis()->SetTitle("M(tau,j)");
  hist_pair_mass->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_MET = new TH1F("met", "Normalized Missing ET", nbins, 0., 1.);
  hist_MET->GetXaxis()->SetTitle("MET/M(tau+,tau-)");
  hist_MET->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_HT_LT = new TH1F("HT_LT", "HT - LT", nbins, -500., 500.);
  hist_HT_LT->GetXaxis()->SetTitle("HT - LT");
  hist_HT_LT->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_ditau_mass = new TH1F("ditau_mass", "M(tau+,tau-)", nbins, 0., 300.);
  hist_ditau_mass->GetXaxis()->SetTitle("M(tau+, tau-)");
  hist_ditau_mass->GetYaxis()->SetTitle("a.u.");

  // Eta and Pt.
  TH1F *hist_eta_btag = new TH1F("eta_btag", "BTag Eta", nbins, -6., 6.);
  hist_eta_btag->GetXaxis()->SetTitle("Eta");
  hist_eta_btag->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_eta_jet = new TH1F("eta_btag", "Secondary Jet Eta", nbins, -6., 6.);
  hist_eta_jet->GetXaxis()->SetTitle("Eta");
  hist_eta_jet->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_eta_tau_p = new TH1F("eta_btag", "Tau+ Eta", nbins, -6., 6.);
  hist_eta_tau_p->GetXaxis()->SetTitle("Eta");
  hist_eta_tau_p->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_eta_tau_m = new TH1F("eta_btag", "Tau- Eta", nbins, -6., 6.);
  hist_eta_tau_m->GetXaxis()->SetTitle("Eta");
  hist_eta_tau_m->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_pt_btag = new TH1F("pt_btag", "BTag Pt", nbins, 0., 300.);
  hist_pt_btag->GetXaxis()->SetTitle("Pt (GeV)");
  hist_pt_btag->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_pt_jet = new TH1F("pt_btag", "Secondary Jet Pt", nbins, 0., 300.);
  hist_pt_jet->GetXaxis()->SetTitle("Pt (GeV)");
  hist_pt_jet->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_pt_tau_p = new TH1F("pt_tau_p", "Tau+ Pt", nbins, 0., 300.);
  hist_pt_tau_p->GetXaxis()->SetTitle("Pt (GeV)");
  hist_pt_tau_p->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_pt_tau_m = new TH1F("pt_tau_m", "Tau- Pt", nbins, 0., 300.);
  hist_pt_tau_m->GetXaxis()->SetTitle("Pt (GeV)");
  hist_pt_tau_m->GetYaxis()->SetTitle("a.u.");

  // Main event loop.
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
      printf("Event # %d: did not find at least 4 jets, skipping... \n", entry);
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

    // Now that we have ID'd our 4 particles, calculate kinematics.
    // p = +, m = -
    // b = btag jet, j = other leading jet (btag or non-btag)
    TLorentzVector tau_p_b = tau_plus + btag_jet;
    TLorentzVector tau_m_j = tau_minus + second_jet;
    TLorentzVector tau_p_j = tau_plus + second_jet;
    TLorentzVector tau_m_b = tau_minus + btag_jet;
    TLorentzVector ditau = tau_plus + tau_minus;

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

    // Fill ditau pair mass.
    hist_ditau_mass->Fill(ditau.M());

    // Fill eta spectrums.
    hist_eta_btag->Fill(btag_jet.Eta());
    hist_eta_jet->Fill(second_jet.Eta());
    hist_eta_tau_p->Fill(tau_plus.Eta());
    hist_eta_tau_m->Fill(tau_minus.Eta());

    // Fill Pt spectrums.
    hist_pt_btag->Fill(btag_jet.Pt());
    hist_pt_jet->Fill(second_jet.Pt());
    hist_pt_tau_p->Fill(tau_plus.Pt());
    hist_pt_tau_m->Fill(tau_minus.Pt());

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
  hist_ditau_mass->Scale(1/hist_ditau_mass->Integral());
  hist_ditau_mass->Draw("HIST");

  TCanvas *c5 = new TCanvas();
  hist_eta_btag->Scale(1/hist_eta_btag->Integral());
  hist_eta_btag->Draw("HIST");

  TCanvas *c6 = new TCanvas();
  hist_eta_jet->Scale(1/hist_eta_jet->Integral());
  hist_eta_jet->Draw("HIST");

  TCanvas *c7 = new TCanvas();
  hist_eta_tau_p->Scale(1/hist_eta_tau_p->Integral());
  hist_eta_tau_p->Draw("HIST");

  TCanvas *c8 = new TCanvas();
  hist_eta_tau_m->Scale(1/hist_eta_tau_m->Integral());
  hist_eta_tau_m->Draw("HIST");

  TCanvas *c9 = new TCanvas();
  hist_pt_btag->Scale(1/hist_pt_btag->Integral());
  hist_pt_btag->Draw("HIST");

  TCanvas *c10 = new TCanvas();
  hist_pt_jet->Scale(1/hist_pt_jet->Integral());
  hist_pt_jet->Draw("HIST");

  TCanvas *c11 = new TCanvas();
  hist_pt_tau_p->Scale(1/hist_pt_tau_p->Integral());
  hist_pt_tau_p->Draw("HIST");

  TCanvas *c12 = new TCanvas();
  hist_pt_tau_m->Scale(1/hist_pt_tau_m->Integral());
  hist_pt_tau_m->Draw("HIST");

  std::string mass_name = string(sample_desc) + "_taujet_mass";
  std::string met_name = string(sample_desc) + "_met";
  std::string HT_LT_name = string(sample_desc) + "_HT_LT";
  std::string ditau_mass_name = string(sample_desc) + "_ditau_mass";

  std::string eta_btag_name = string(sample_desc) + "_eta_btag";
  std::string eta_jet_name = string(sample_desc) + "_eta_jet";
  std::string eta_taup_name = string(sample_desc) + "_eta_taup";
  std::string eta_taum_name = string(sample_desc) + "_eta_taum";

  std::string pt_btag_name = string(sample_desc) + "_pt_btag";
  std::string pt_jet_name = string(sample_desc) + "_pt_jet";
  std::string pt_taup_name = string(sample_desc) + "_pt_taup";
  std::string pt_taum_name = string(sample_desc) + "_pt_taum";


  TFile *f1 = new TFile("hist_taujet_mass.root", "UPDATE");
  hist_pair_mass->Write(mass_name.c_str(), TObject::kOverwrite);

  TFile *f2 = new TFile("hist_MET.root", "UPDATE");
  hist_MET->Write(met_name.c_str(), TObject::kOverwrite);

  TFile *f3 = new TFile("hist_HT_LT.root", "UPDATE");
  hist_HT_LT->Write(HT_LT_name.c_str(), TObject::kOverwrite);

  TFile *f4 = new TFile("hist_ditau_mass.root", "UPDATE");
  hist_ditau_mass->Write(ditau_mass_name.c_str(), TObject::kOverwrite);

  TFile *f5 = new TFile("hist_eta_btag.root" ,"UPDATE");
  hist_eta_btag->Write(eta_btag_name.c_str(), TObject::kOverwrite);

  TFile *f6 = new TFile("hist_eta_jet.root" ,"UPDATE");
  hist_eta_jet->Write(eta_jet_name.c_str(), TObject::kOverwrite);

  TFile *f7 = new TFile("hist_eta_taup.root" ,"UPDATE");
  hist_eta_tau_p->Write(eta_taup_name.c_str(), TObject::kOverwrite);

  TFile *f8 = new TFile("hist_eta_taum.root" ,"UPDATE");
  hist_eta_tau_m->Write(eta_taum_name.c_str(), TObject::kOverwrite);

  TFile *f9 = new TFile("hist_pt_btag.root", "UPDATE");
  hist_pt_btag->Write(pt_btag_name.c_str(), TObject::kOverwrite);

  TFile *f10 = new TFile("hist_pt_jet.root", "UPDATE");
  hist_pt_jet->Write(pt_jet_name.c_str(), TObject::kOverwrite);

  TFile *f11 = new TFile("hist_pt_taup.root", "UPDATE");
  hist_pt_tau_p->Write(pt_taup_name.c_str(), TObject::kOverwrite);

  TFile *f12 = new TFile("hist_pt_taum.root", "UPDATE");
  hist_pt_tau_m->Write(pt_taum_name.c_str(), TObject::kOverwrite);

}
