// An analyzer to produce the b-tagged jet / b-tagged + mu PT ratio.
// USAGE:
// gSystem->Load("<path_to_delphes>/Delphes/libDelphes.so");
// .x PTRatioAnalyzer.C("input_file.root", "saved_histogram_name.root", "plot_name", nbins);

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif

void PTRatioAnalyzer(const char *file_name, const char *save_name, const char *plot_name, int nbins$
  gSystem->Load("/home/thompson/MG5_aMC_2_6_0/Delphes/libDelphes.so");

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
  TH1F *hist_MET = new TH1F(plot_name + ": Normalized Missing ET, "", nbins, 0., 1.);
  hist_MET->GetXaxis()->SetTitle("MET/Mmumu");
  hist_MET->GetYaxis()->SetTitle("Nb events");

  TH1F *hist_lj_mass = new TH1F(plot_name + ": M(lj)", "", nbins, 0., 1.);
  hist_lj_mass->GetXaxis()->SetTitle("M(mu,j)");
  hist_lj_mass->GetYaxis()->SetTitle("Nb events");

  TH1F *hist_HT_LT = new TH1F(plot_name + ": HT - LT", "", nbins, 0., 1.);
  hist_ht_lt->GetXaxis()->SetTitle("HT - LT");
  hist_ht_lt->GetYaxis()->SetTitle("Nb events");

  TH1F *hist_zeta = new TH1F(plot_name + ": Zeta", "", nbins, 0., 1.);
  hist_zeta->GetXaxis()->SetTitle("zeta");
  hist_zeta->GetYaxis()->SetTitle("Nb events");

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
    if (jet_size < 4) continue;
   
    // Skip events if they do not have a BTag jet.
    // Skip events if they do not have a Tau+ and Tau-.
    
    // Declare 4-vectors.
    TLorentzVector *btag_jet;
    TLorentzVector *second_jet;
    TLorentzVector *tau_plus;
    TLorentzVector *tau_minus;
    
    Int_t leading_btag_id = -1;
    // Top mass bound and HT - LT loop.
    // Loop over jets and find the highest pt BTag, the second highest pt, non-TauTag jet,
    // and two opposite sign TauTag jets.
    for (Int_t ii = 0; ii < jet_size; ii++) {
      jet = (Jet*) branch_jet->At(ii);
      // Find the leading b-tagged jet.
      if (jet->BTag == 1 && jet->PT > btag_jet->PT()) {
        btag_jet->SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->M);
        leading_btag_id = ii;
      }
      // Find the leading OS tau pair.
      if (jet->TauTag == 1) {
        if (jet->Charge == 1 && jet->PT > tau_plus->PT()) {
          tau_plus->SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->M);
        }
        if (jet->Charge == -1 && jet->PT > tau_minus->PT()) {
          tau_minus->SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->M);
        }
      }
    }
    
    // Loop to find the other highest-pt jet with ID different than the leading BTag. 
    for (Int_t ii = 0; ii < jet_size; ii++) {
      jet = (Jet*) branch_jet->At(ii);
      if (jet->TauTag == 0 && jet->PT > second_jet->PT() && ii != leading_btag_id) {
        second_jet->SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->M);
      }
    }
    
    // Determine tau-jet pair mass.
    
    // Normalized Missing ET loop.
    
    // Calculate HT - LT.
    
    

  }

  // Draw histograms and save them in a .root format.
  TCanvas *c = new TCanvas();

  hist_lj_mass->Scale(1/hist_lj_mass->Integral());
  hist_lj_mass->Draw();

  TFile *f = new TFile(save_name,"RECREATE");
  hist_lj_mass->Write("",TObject::kOverwrite);
}
