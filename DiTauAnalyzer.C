// An analyzer to produce key selection criteria plots for Z' -> ditau processes.
// arXiv:1707.07016v1
// USAGE:
// gSystem->Load("<path_to_delphes>/Delphes/libDelphes.so");
// .x DiTauAnalyzer.C("input_file.root", "Z' 500 GeV", nbins);

#ifdef __CLING__
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

  TH1F *hist_tmass = new TH1F("tmass", "Transverse (Ditau + MET) mass", nbins, 0., 500.);
  hist_tmass->GetXaxis()->SetTitle("MT2(tau+,tau-,met) [GeV]");
  hist_tmass->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_HT_LT = new TH1F("HT_LT", "HT - LT", nbins, -500., 500.);
  hist_HT_LT->GetXaxis()->SetTitle("HT - LT");
  hist_HT_LT->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_ditau_mass = new TH1F("ditau_mass", "M(tau+,tau-)", nbins, 0., 500.);
  hist_ditau_mass->GetXaxis()->SetTitle("M(tau+, tau-)");
  hist_ditau_mass->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_pt_btag = new TH1F("pt_btag", "BTag Pt", nbins, 0., 300.);
  hist_pt_btag->GetXaxis()->SetTitle("Pt (GeV)");
  hist_pt_btag->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_pt_jet = new TH1F("pt_jet", "Secondary Jet Pt", nbins, 0., 300.);
  hist_pt_jet->GetXaxis()->SetTitle("Pt (GeV)");
  hist_pt_jet->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_pt_tau_p = new TH1F("pt_tau_p", "Tau+ Pt", nbins, 0., 300.);
  hist_pt_tau_p->GetXaxis()->SetTitle("Pt (GeV)");
  hist_pt_tau_p->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_pt_tau_m = new TH1F("pt_tau_m", "Tau- Pt", nbins, 0., 300.);
  hist_pt_tau_m->GetXaxis()->SetTitle("Pt (GeV)");
  hist_pt_tau_m->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_pt_ditau = new TH1F("pt_ditau", "Ditau Pt", nbins, 0., 300.);
  hist_pt_ditau->GetXaxis()->SetTitle("Pt (GeV)");
  hist_pt_ditau->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_dzeta = new TH1F("dzeta", "DZeta, #alpha = 0.50", nbins, -300., 200.);
  hist_dzeta->GetXaxis()->SetTitle("DZeta [GeV]");
  hist_dzeta->GetYaxis()->SetTitle("a.u.");

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
    if (jet_size < 4) continue;

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

    // Make Pre-selection cuts.
    if (!found_tau_plus || !found_tau_minus) continue;
    if (!found_btag) continue;

    // Now that we have ID'd our 4 particles, calculate kinematics.
    // p = +, m = -
    // b = btag jet, j = other leading jet (btag or non-btag)
    TLorentzVector tau_p_b = tau_plus + btag_jet;
    TLorentzVector tau_m_j = tau_minus + second_jet;
    TLorentzVector tau_p_j = tau_plus + second_jet;
    TLorentzVector tau_m_b = tau_minus + btag_jet;
    TLorentzVector ditau = tau_plus + tau_minus;
    TLorentzVector met_p4, transverse_tau_plus, transverse_tau_minus;

    // Calculate and fill tau-jet pair mass by taking the permutation of tau-jet
    // pairs (where jet = b or non-b) with the smallest mass difference, then
    // picking the highest pair mass of that permutation.
    Double_t choice_pair_mass;
    if (abs(tau_p_b.M() - tau_m_j.M()) < abs(tau_p_j.M() - tau_m_b.M())) {
      choice_pair_mass = std::max(tau_p_b.M(), tau_m_j.M());
    } else {
      choice_pair_mass = std::max(tau_p_j.M(), tau_m_b.M());
    }
    hist_pair_mass->Fill(choice_pair_mass);

    // Calculate and fill Missing ET normalized to the ditau mass.
    met = (MissingET*) branch_met->At(0);
    printf("MET is %f \n", met->MET);
    hist_MET->Fill(met->MET / ditau.M());

    // Fill ditau + MET transverse mass.
    met_p4.SetPtEtaPhiE(met->MET, 0, met->Phi, met->MET);
    transverse_tau_plus.SetPtEtaPhiE(tau_plus.Pt(), 0., tau_plus.Phi(),
                                      tau_plus.Et());
    transverse_tau_minus.SetPtEtaPhiE(tau_minus.Pt(), 0., tau_minus.Phi(),
                                       tau_minus.Et());
    TLorentzVector transverse_p4 = met_p4 + transverse_tau_plus
                                 + transverse_tau_minus;
    hist_tmass->Fill(transverse_p4.M());

    // Calculate D-Zeta to distinguish W+jets backgrounds.
    // (these are actually 2-component vectors in the Eta=0 plane)
    TVector3 p_vis_tau_plus = transverse_tau_plus.Vect();
    TVector3 p_vis_tau_minus = transverse_tau_minus.Vect();
    TVector3 zeta = (p_vis_tau_minus.Mag() * p_vis_tau_plus
                     + p_vis_tau_plus.Mag() * p_vis_tau_minus);
    TVector3 zeta_hat = zeta * (1 / zeta.Mag());
    Double_t p_vis = p_vis_tau_plus.Dot(zeta_hat) + p_vis_tau_minus.Dot(zeta_hat);
    Double_t p_miss = (met_p4.Vect()).Dot(zeta_hat);

    hist_dzeta->Fill(p_miss - 0.50 * p_vis); // alpha = 0.50 (Tau_hadronic).

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
  hist_tmass->Scale(1/hist_tmass->Integral());
  hist_tmass->Draw("HIST");

  TCanvas *c6 = new TCanvas();
  hist_pt_btag->Scale(1/hist_pt_btag->Integral());
  hist_pt_btag->Draw("HIST");

  TCanvas *c7 = new TCanvas();
  hist_pt_jet->Scale(1/hist_pt_jet->Integral());
  hist_pt_jet->Draw("HIST");

  TCanvas *c8 = new TCanvas();
  hist_pt_tau_p->Scale(1/hist_pt_tau_p->Integral());
  hist_pt_tau_p->Draw("HIST");

  TCanvas *c9 = new TCanvas();
  hist_pt_tau_m->Scale(1/hist_pt_tau_m->Integral());
  hist_pt_tau_m->Draw("HIST");

  TCanvas *c10 = new TCanvas();
  hist_pt_ditau->Scale(1/hist_pt_ditau->Integral());
  hist_pt_ditau->Draw("HIST");

  TCanvas *c11 = new TCanvas();
  hist_dzeta->Scale(1/hist_dzeta->Integral());
  hist_dzeta->Draw("HIST");

  std::string mass_name = string(sample_desc) + " Tau-Jet Mass";
  std::string met_name = string(sample_desc) + " MET";
  std::string tmass_name = string(sample_desc) + " Transverse Mass";
  std::string HT_LT_name = string(sample_desc) + " HT - LT";
  std::string ditau_mass_name = string(sample_desc) + " ditau Mass";
  std::string pt_btag_name = string(sample_desc) + " BTag PT";
  std::string pt_jet_name = string(sample_desc) + " Secondary Jet PT";
  std::string pt_taup_name = string(sample_desc) + " Tau+ PT";
  std::string pt_taum_name = string(sample_desc) + " Tau- PT";
  std::string pt_ditau_name = string(sample_desc) + " Ditau PT";
  std::string dzeta_name = string(sample_desc) + " DZeta";

  TFile *f1 = new TFile("hist_taujet_mass.root", "UPDATE");
  hist_pair_mass->Write(mass_name.c_str(), TObject::kOverwrite);

  TFile *f2 = new TFile("hist_MET.root", "UPDATE");
  hist_MET->Write(met_name.c_str(), TObject::kOverwrite);

  TFile *f3 = new TFile("hist_HT_LT.root", "UPDATE");
  hist_HT_LT->Write(HT_LT_name.c_str(), TObject::kOverwrite);

  TFile *f4 = new TFile("hist_ditau_mass.root", "UPDATE");
  hist_ditau_mass->Write(ditau_mass_name.c_str(), TObject::kOverwrite);

  TFile *f5 = new TFile("hist_tmass.root" ,"UPDATE");
  hist_tmass->Write(tmass_name.c_str(), TObject::kOverwrite);

  TFile *f6 = new TFile("hist_pt_btag.root", "UPDATE");
  hist_pt_btag->Write(pt_btag_name.c_str(), TObject::kOverwrite);

  TFile *f7 = new TFile("hist_pt_jet.root", "UPDATE");
  hist_pt_jet->Write(pt_jet_name.c_str(), TObject::kOverwrite);

  TFile *f8 = new TFile("hist_pt_taup.root", "UPDATE");
  hist_pt_tau_p->Write(pt_taup_name.c_str(), TObject::kOverwrite);

  TFile *f9 = new TFile("hist_pt_taum.root", "UPDATE");
  hist_pt_tau_m->Write(pt_taum_name.c_str(), TObject::kOverwrite);

  TFile *f10 = new TFile("hist_pt_ditau.root", "UPDATE");
  hist_pt_ditau->Write(pt_ditau_name.c_str(), TObject::kOverwrite);

  TFile *f11 = new TFile("hist_dzeta.root", "UPDATE");
  hist_dzeta->Write(dzeta_name.c_str(), TObject::kOverwrite);
}
