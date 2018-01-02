// An analyzer to produce key selection criteria plots for dimuon processes.
// EXAMPLE USAGE:
// gSystem->Load("<path_to_delphes>/Delphes/libDelphes.so");
// .x DiMuonAnalyzer.C("input_file.root", "Zprime (500 GeV)", nbins);

#ifdef __CLING__
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "lester_mt2_bisect.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif

// Define global constants.
Double_t kMuonMass = 0.105658389;  // GeV
Double_t previous_mt2 = 999.0;

void DiMuonAnalyzer(const char *file_name, const char *sample_desc, int nbins) {
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
  TClonesArray *branch_muon = tree_reader->UseBranch("Muon");

  // Book histograms.
  TH1F *hist_pair_mass = new TH1F("pair_mass", "Max{M(j,mu)}", nbins, 0., 300.);
  hist_pair_mass->GetXaxis()->SetTitle("M(mu,j)");
  hist_pair_mass->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_MET = new TH1F("met", "Normalized Missing ET", nbins, 0., 1.);
  hist_MET->GetXaxis()->SetTitle("MET/Mmumu");
  hist_MET->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_HT_LT = new TH1F("HT_LT", "HT - LT", nbins, -500., 500.);
  hist_HT_LT->GetXaxis()->SetTitle("HT - LT");
  hist_HT_LT->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_MT2 = new TH1F("MT2", "MT2", nbins, 0., 100.);
  hist_MT2->GetXaxis()->SetTitle("M_{T2}");
  hist_MT2->GetYaxis()->SetTitle("a.u.");

  TH2F *hist_MT2_dPhi = new TH2F("mt2_vs_dPhi", "Max #delta #phi vs. M_{T2}", nbins,
                                0., 100., nbins, 0., 3.14);
  hist_MT2_dPhi->GetXaxis()->SetTitle("M_{T2}");
  hist_MT2_dPhi->GetYaxis()->SetTitle("max[dPhi(mu+, MET), dPhi(mu-, MET)]");

  TH1F *hist_pt_mu = new TH1F("mupt", "", nbins, 0., 300.);

  // Event loop.
  for(Int_t entry = 0; entry < number_of_entries; ++entry) {

    tree_reader->ReadEntry(entry);

    // Declare physics objects.
    Muon *muon;
    Jet *jet;
    MissingET *met;
    TLorentzVector btag_jet;
    TLorentzVector second_jet;
    TLorentzVector mu_plus;
    TLorentzVector mu_minus;

    // Declare running indices.
    Int_t jet_size = branch_jet->GetEntries();
    Int_t met_size = branch_met->GetEntries();
    Int_t muon_size = branch_muon->GetEntries();

    // Redeclare preselection states.
    if (jet_size < 2) continue;
    bool found_btag = false;
    bool found_mu_plus = false;
    bool found_mu_minus = false;

    // Find the leading OS Muon pair.
    for (Int_t ii = 0; ii < muon_size; ii++) {
      muon = (Muon*) branch_muon->At(ii);
      if (muon->Charge == 1 && muon->PT > mu_plus.Pt()) {
        mu_plus.SetPtEtaPhiM(muon->PT, muon->Eta, muon->Phi, kMuonMass);
        found_mu_plus = true;
      }
      if (muon->Charge == -1 && muon->PT > mu_minus.Pt()) {
        mu_minus.SetPtEtaPhiM(muon->PT, muon->Eta, muon->Phi, kMuonMass);
        found_mu_minus = true;
      }
    }

    // Loop over jets and find the highest pt BTag jet.
    Int_t leading_btag_id = -1;
    for (Int_t ii = 0; ii < jet_size; ii++) {
      jet = (Jet*) branch_jet->At(ii);
      // Find the leading b-tagged jet.
      if (jet->BTag == 1 && jet->PT > btag_jet.Pt()) {
        btag_jet.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
        leading_btag_id = ii;
        found_btag = true;
      }
    }

    // Find the other highest-pt jet with ID different than the leading BTag.
    for (Int_t ii = 0; ii < jet_size; ii++) {
      jet = (Jet*) branch_jet->At(ii);
      if (jet->TauTag == 0 && jet->PT > second_jet.Pt() && ii != leading_btag_id) {
        second_jet.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
      }
    }

    // Place preselection cuts.
    if (!found_btag) continue;
    if (!found_mu_plus || !found_mu_minus) continue;
    // Place kinematic cuts.
    if (mu_plus.Pt() < 30.) continue;
    if (mu_minus.Pt() < 30.) continue;
    if (btag_jet.Pt() == 0.) continue;
    if (second_jet.Pt() == 0.) continue;

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

    hist_pt_mu->Fill(mu_plus.Pt());

    // Calculate and fill muon/jet pair mass.
    Double_t choice_pair_mass;
    if (abs(mu_p_b.M() - mu_m_j.M()) < abs(mu_p_j.M() - mu_m_b.M())) {
      choice_pair_mass = std::max(mu_p_b.M(), mu_m_j.M());
    } else {
      choice_pair_mass = std::max(mu_p_j.M(), mu_m_b.M());
    }
    hist_pair_mass->Fill(choice_pair_mass);

    // Calculate and fill MET for the event.
    met = (MissingET*) branch_met->At(0);
    hist_MET->Fill(met->MET / dimuon.M());

    // Calculate MT2 (http://www.hep.phy.cam.ac.uk/~lester/mt2/).
    met_P4.SetPtEtaPhiE(met->MET, 0, met->Phi, met->MET);
    Double_t mVisA = mu_plus.M(); // Mass of visible object on side A.
    Double_t pxA = mu_plus.Px(); // x momentum of visible object on side A.
    Double_t pyA = mu_plus.Py(); // y momentum of visible object on side A.

    Double_t mVisB = mu_minus.M(); // Mass of visible object on side B.
    Double_t pxB = mu_minus.Px(); // x momentum of visible object on side B.
    Double_t pyB = mu_minus.Py(); // y momentum of visible object on side B.

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

    // Report lowest MT2 value.
    if (MT2 < previous_mt2) {
      previous_mt2 = MT2;
      printf("%f \n", MT2);
    }

    // Calculate and fill HT - LT.
    Double_t H_T = abs(btag_jet.Pt()) + abs(second_jet.Pt());
    Double_t L_T = abs(mu_plus.Pt()) + abs(mu_minus.Pt());
    hist_HT_LT->Fill(H_T - L_T);

    // Calculate dphi/mt2
    Double_t max_dphi = std::max(abs(mu_plus.DeltaPhi(met_P4)),
                                 abs(mu_minus.DeltaPhi(met_P4)));

    hist_MT2_dPhi->Fill(MT2, max_dphi);

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
  hist_MT2->Scale(1/hist_MT2->Integral());
  hist_MT2->Draw("HIST");

  gStyle->SetOptStat(0);
  TCanvas *c5 = new TCanvas();
  hist_MT2_dPhi->Draw("COL2Z");

  TCanvas *c6 = new TCanvas();
  hist_pt_mu->Scale(1/hist_pt_mu->Integral());
  hist_pt_mu->Draw("HIST");

  std::string mass_name = string(sample_desc) + " max(M(#mu,j))";
  std::string met_name = string(sample_desc) + " MET";
  std::string HT_LT_name = string(sample_desc) + " H_{T} - L_{T}";
  std::string mt2_name = string(sample_desc) + " MT2";

  TFile *f1 = new TFile("hist_dimuon_mujet_mass.root", "UPDATE");
  hist_pair_mass->Write(mass_name.c_str(), TObject::kOverwrite);

  TFile *f2 = new TFile("hist_dimuon_MET.root", "UPDATE");
  hist_MET->Write(met_name.c_str(), TObject::kOverwrite);

  TFile *f3 = new TFile("hist_dimuon_HT_LT.root", "UPDATE");
  hist_HT_LT->Write(HT_LT_name.c_str(), TObject::kOverwrite);

  TFile *f4 = new TFile("hist_dimuon_MT2.root" ,"UPDATE");
  hist_MT2->Write(mt2_name.c_str(), TObject::kOverwrite);

  TFile *f5 = new TFile("hist_dimuon_ptmu.root", "UPDATE");
  hist_pt_mu->Write(sample_desc, TObject::kOverwrite);

}
