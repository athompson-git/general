#define zprime_BFF_dimuon_analysis_cxx

#include "zprime_BFF_dimuon_analysis.h"
#include <TH2.h>
#include <TStyle.h>

void zprime_BFF_dimuon_analysis::Begin(TTree * /*tree*/)
{
   TString option = GetOption();
}

void zprime_BFF_dimuon_analysis::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();
   // Define a TH1F 1D histogram to fill with the ditau invariant mass.
   hist = new TH1F("DimuonMass", "Dimuon Invariant Mass, g g > zp b b~ $$ b b~, (zp > mu- mu+), MZP=200 GeV", 100, 0.0, 1000.0);
   hist->SetLineColor(1);
   hist->GetXaxis()->SetTitle("Dimuon Invariant Mass M_mu-mu+ (GeV)");
   hist->GetYaxis()->SetTitle("Number of events");
}

Bool_t zprime_BFF_dimuon_analysis::Process(Long64_t entry)
{
   fReader.SetEntry(entry);
   GetEntry(entry);

   // Loop over particles in the event and find Tau pairs.
   // Calculate the invariant mass from the Tau pairs to reconstruct the Z prime.
   while(fReader.Next()) {
     int nParticles = Event_Nparticles[0];
     Int_t mu = 0;
     Int_t amu = 0;
     TLorentzVector dimuon;

     for(int p = 0; p < nParticles; ++p) {
       // Muon PID = +/- 13
       if(Particle_PID[p] == -13) {
         mu = p;
       }

       if(Particle_PID[p] == 13) {
         amu = p;
       }
     }

     // Build 4-vectors for each Muon.
     TLorentzVector fmomentum_mu(Particle_Px[mu], Particle_Py[mu], Particle_Pz[mu], Particle_E[mu]);
     TLorentzVector fmomentum_amu(Particle_Px[amu], Particle_Py[amu], Particle_Pz[amu], Particle_E[amu]);
     dimuon = fmomentum_mu + fmomentum_amu;

     hist->Fill(dimuon.M());
   }

   return kTRUE;
}

void zprime_BFF_dimuon_analysis::SlaveTerminate()
{
   hist->Draw();
}

void zprime_BFF_dimuon_analysis::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
