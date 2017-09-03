//derived from Delphes Example1 as stated directly below
/*
Simple macro showing how to access branches from the delphes output root file,
loop over events, and plot simple quantities such as the jet pt and the di-Muon invariant
mass.

root -l examples/Example1.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include <TROOT.h>
#include <TString.h>
#include <TVector2.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <TProfile.h>

#endif

//------------------------------------------------------------------------------

//example usage: root -l -b -q 'Cutflow.C(-3,1)' gives the result for requring at least 2 jets with at least 1 bottom-tag for the 200 GeV Z'+1b signal sample in the result file zpb_200.root
void Cutflow(int input, int BottomReq)
{
  bool verbose_=false;
  gSystem->Load("libDelphes");
  float cs=0.;

  //define input files here
  vector<TString> fNames;
  TString savename;
  if(input==-3)     {fNames.push_back("200oneb.root"); cs = 0.020;   savename="zpb_200.root";} //old 0.037
  else if(input==-2){fNames.push_back("350oneb.root"); cs = 0.002;   savename="zpb_350.root";} //old 0.008//0.001
  else if(input==-1){fNames.push_back("500oneb.root"); cs = 0.0012;  savename="zpb_500.root";} //old 0.0012//0.0001
  else if(input==0) {fNames.push_back("tt/tt_*.root"); cs = 36.;     savename="ttbar.root";   }//Peisi-space
  else if(input==1) {fNames.push_back("200GeV.root");  cs = 0.040;   savename="zpbb_200.root";}//old 0.08//200 GeV, Peisi-space
  else if(input==2) {fNames.push_back("350GeV.root");  cs = 0.004;   savename="zpbb_350.root";}//old 0.017//350 GeV, Peisi-space //0.001
  else if(input==3) {fNames.push_back("500GeV.root");  cs = 0.0024;  savename="zpbb_500.root";}//old 0.002//500 GeV, Peisi-space //0.0001
  else if(input==4) {fNames.push_back("mumu_bb_offshellZ/Events/delphes_root/run*_delphes_events.root");   cs = 3.46; savename="offZ_mumu_bb.root";}//Mysha-space
  else if(input==5) {fNames.push_back("mumu_jets_offshellZ/Events/delphes_root/run*_delphes_events.root"); cs = 254.; savename="offZ_mumu_qq.root";}//Mysha-space
  else if(input==6) {fNames.push_back("z_ll_012jets_matched/Events/delphes_root/z_ll_012jets_part*_delphes_events.root"); cs = 3546.; savename="onZ_mumu_qq.root";}//Mysha-space
  else if(input==7) {fNames.push_back("z_ll_bb_matched/Events/delphes_root/tag_1_delphes_events.root");     cs = 61.4; savename="onZ_mumu_bb.root";}//Mysha-space
  if(BottomReq==2) savename.Prepend("2b/");
  else if(BottomReq==1) savename.Prepend("1b/");
  else if(!BottomReq) savename.Prepend("0b/");
  TString prefix;
  if(input<4) prefix="/fdata/hepx/store/user/zprime/"; //Peisi's stuff
  else prefix="/fdata/scratch/mexanick/"; //Mysha's stuff
  TChain chain_("Delphes");
  
  for(unsigned fIdx=0; fIdx<fNames.size(); ++fIdx){
    // Create chain of root trees
    fNames[fIdx].Prepend(prefix);
    if(verbose_)cout<<fNames[fIdx]<<endl;
    chain_.Add(fNames[fIdx]);
  }

  if(verbose_) cout<<"enter Analysis"<<endl;
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain_);
  Long64_t numberOfEntries = treeReader->GetEntries();
  //numberOfEntries=10000; //comment in for test purposes
  float weight = 1000.*cs/numberOfEntries;
  std::cout<<"event weight: "<<weight<<std::endl;

  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMPT = treeReader->UseBranch("MissingET");

  // Book histograms
  TH1 *histJetPT = new TH1F("jet_pt", "jet P_{T}", 100, 0.0, 1000.0);
  TH1 *histJetPT_nonB2 = new TH1F("jet2_nonB", "jet2 P_{T}, not b-tagged;p_{T} [GeV];events/fb^{-1}", 100, 0.0, 1000.0);
  TH1 *histBjetPT = new TH1F("Bjet_pt", "Bjet P_{T}", 100, 0.0, 1000.0);
  TH1 *histMass = new TH1F("mass", "OS leading di-muon mass;M_{inv}(#mu_{1}, #mu_{2});events/fb^{-1}", 100, 0.0, 1000.0);
  TH1 *histMPT = new TH1F("histMPT","missing p_{T};MPT;events/fb^{-1}", 100, 0.0, 1000.0);
  TH2 *histSimilarBranchMass = new TH2F("SimilarBranchMass","most similar mu+b branches permutation mass;SBM [GeV];events/fb^{-1}", 100, 0., 300., 100, 0., 300.);
  TH2 *histDirectionalDeviation = new TH2F("DirectionalDeviation","relative transverse distance of MPT to di-muon system;MPT/p_{T}(#mu#mu)*sin(#Delta#phi);MPT/p_{T}(#mu#mu)*cos(#Delta#phi)",100,-10.,10.,100,-10.,10.);
  TH2 *histDeltaHTLT = new TH2F("DeltaHTLT","difference between hadronic and leptonic scalar transverse momenta sum;HT;LT",100.,0.,1000.,100.,0.,1000.);
  TH1 *histDHTLT = new TH1F("DHTLT","HT-LT;HT-LT;events/fb^{-1}",100,-500.,500.);
  TH1 *massAfterCuts = new TH1F("massAfterCuts","mass after DHTLT<0, MET/M_{#mu#mu}<0.2, max(SBM)>170;M_{#mu#mu} [GeV];events/fb^{-1}",200,0.,1000.);
  massAfterCuts->Sumw2();
  TH1 *massAfterCuts_NoB_2j = new TH1F("massAfterCuts_NoB_2j","mass after DHTLT<0, MET/M_{#mu#mu}<0.2, max(SBM)>170;M_{#mu#mu} [GeV];events/fb^{-1}",100,0.,1000.);
  massAfterCuts_NoB_2j->Sumw2();
  TH1 *DSBM = new TH1F("DSBM","#Delta SBM;|#Delta SBM| [GeV];events/fb^{-1}",100,0.,300.);
  TH1 *maxSBM = new TH1F("maxSBM","max(SBM);SBM;events/fb^{-1}",100,0.,300.);
  TH1 *minSBM = new TH1F("minSBM","min(SBM);SBM;events/fb^{-1}",100,0.,300.);
  TH1 *METvLT = new TH1F("METvLT","MET normalized to scalar sum of lepton momenta after max(SBM)>170 & DHTLT<0",100,0.,10.);
  TH1 *Proj = new TH1F("Proj","MPF;MPF;events/fb^{-1}",61,-3.05,3.05);
  TH1 *METvsMmm = new TH1F("METvsMmm","MET/M^{#mu#mu};MET/M^{#mu#mu};events/fb^{-1}",100,0.,1.);
  /*TH1 *Jet1_dR_1b = new TH1F("Jet1_dR_1b","angular distance between non-btagged jet and closest muon;events/fb^{-1}",100,0.,0.5);
  TH1 *Jet2_dR_1b = new TH1F("Jet2_dR_1b","angular distance between non-btagged jet and closest muon;events/fb^{-1}",100,0.,0.5);
  TH1 *Jet1_dR_2b = new TH1F("Jet1_dR_2b","angular distance between non-btagged jet and closest muon;events/fb^{-1}",100,0.,0.5);
  TH1 *Jet2_dR_2b = new TH1F("Jet2_dR_2b","angular distance between non-btagged jet and closest muon;events/fb^{-1}",100,0.,0.5);*/
  TProfile *cutflow = new TProfile("cutflow","2b, 2 OS muons, MET<100, DHTLT<0, max(SBM)>170, DSBM>50, ;cut number; percent passing",6,-0.5,5.5);
  /*TProfile *cutflow_ABC = new TProfile("cutflow_ABC","(2b & 2OS muons) & MET<100, DHTLT<0, max(SBM)>170;cut number;percent passing",4,-0.5,3.5);
  TProfile *cutflow_ACB = new TProfile("cutflow_ACB","(2b & 2OS muons) & MET<100, max(SBM)>170, DHTLT<0;cut number;percent passing",4,-0.5,3.5);
  TProfile *cutflow_BAC = new TProfile("cutflow_BAC","(2b & 2OS muons) & DHTLT<0, MET<100, max(SBM)>170;cut number;percent passing",4,-0.5,3.5);
  TProfile *cutflow_BCA = new TProfile("cutflow_BCA","(2b & 2OS muons) & DHTLT<0, max(SBM)>170, MET<100;cut number;percent passing",4,-0.5,3.5);*/
  TProfile *cutflow_CAB = new TProfile("cutflow_CAB","(2b & 2OS muons) & max(SBM)>170, MET<100, DHTLT<0;cut number;percent passing",4,-0.5,3.5);
  //TProfile *cutflow_CBA = new TProfile("cutflow_CBA","(2b & 2OS muons) & max(SBM)>170, DHTLT<0, MET<100;cut number;percent passing",4,-0.5,3.5);

  TH2F *Mmm_METvsMmm  = new TH2F("Mmm_METvsMmm","M(#mu^{+}#mu^{-}) vs E_{T}^{miss}/M(#mu^{+}#mu^{-});E_{T}*{miss}/M(#mu^{+}#mu^{-});M(#mu^{+}#mu^{-}) [GeV];events/fb^{-1}",100,0.,1.,200,0.,1000.);
  TH2F *Mmm_maxSBM    = new TH2F("Mmm_maxSBM"  ,"M(#mu^{+}#mu^{-}) vs max(SBM);max(SBM) [GeV];M(#mu^{+}#mu^{-}) [GeV];events/fb^{-1}",100,0.,300.,200,0.,1000.);
  TH2F *Mmm_DHTLT     = new TH2F("Mmm_DHTLT"   ,"M(#mu^{+}#mu^{-}) vs H_{T}-L_{T};H_{T}-L_{T} [GeV];M(#mu^{+}#mu^{-}) [GeV];events/fb^{-1}",100,-500.,500.,200,0.,1000.);
  

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    if(entry%50000==0) std::cout<<entry<<"/"<<numberOfEntries<<std::endl;
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    vector<int> bPosition;
    int Jet2_nonB=-1; int Jet3_nonB=-1;
    float Jet2_nonB_PT=0.; float Jet3_nonB_PT=0.;

    // If event contains at least 1 jet
    for(unsigned Nj=0; Nj<branchJet->GetEntries(); ++Nj){
      // Take first jet
      Jet *jet = (Jet*) branchJet->At(Nj);
      if(jet->PT <30.) continue;
      if(jet->BTag){
	bPosition.push_back(Nj);
	histBjetPT->Fill(jet->PT,weight);
      }
      else if(Jet2_nonB_PT<jet->PT){
	  Jet2_nonB=Nj;
	  Jet2_nonB_PT=jet->PT;
      }
      else if(Jet3_nonB_PT<jet->PT){
	  Jet3_nonB=Nj;
	  Jet3_nonB_PT=jet->PT;
      }

      // Plot jet transverse momentum
      histJetPT->Fill(jet->PT,weight);

      // Print jet transverse momentum
      if(verbose_) cout << "Jet pt: "<<jet->PT << endl;
    }
    if(bPosition.size()==1) histJetPT_nonB2->Fill(Jet2_nonB_PT,weight);

    bool DiBottom   = false;
    bool Bottom_Jet = false;
    bool NoB_2Jet   = false;
    if(bPosition.size()>=2) DiBottom=true;
    if(BottomReq==2) cutflow->Fill(0.,DiBottom,weight);
    else if(BottomReq==1){
      if(bPosition.size()==1 && Jet2_nonB_PT>0.) Bottom_Jet =true;
      cutflow->Fill(0.,DiBottom || Bottom_Jet,weight);
    }
    if(!bPosition.size() && Jet3_nonB_PT>0.) bool NoB_2Jet = true;

    Muon *mu1, *mu2;
    bool OSdiMuon=false;

    // If event contains at least 2 Muons
    if(branchMuon->GetEntries() > 1)
    {
      // Take first two Muons
      mu1 = (Muon *) branchMuon->At(0);
      mu2 = (Muon *) branchMuon->At(1);

      if(mu1->Charge!=mu2->Charge)OSdiMuon=true;

      /*cutflow_ABC->Fill(0.,OSdiMuon && (DiBottom || Bottom_Jet),weight);
      cutflow_ACB->Fill(0.,OSdiMuon && (DiBottom || Bottom_Jet),weight);
      cutflow_BAC->Fill(0.,OSdiMuon && (DiBottom || Bottom_Jet),weight);
      cutflow_BCA->Fill(0.,OSdiMuon && (DiBottom || Bottom_Jet),weight);*/
      cutflow_CAB->Fill(0.,OSdiMuon && (DiBottom || Bottom_Jet),weight);
      //cutflow_CBA->Fill(0.,OSdiMuon && (DiBottom || Bottom_Jet),weight);
      cutflow->Fill(1.,OSdiMuon,weight);
      if(!BottomReq){
	if(OSdiMuon){
	  const float mass=((mu1->P4()) + (mu2->P4())).M();
	  massAfterCuts->Fill(mass,weight);
	}
	continue;
      }

      if(!OSdiMuon) continue;
      if(!(DiBottom || Bottom_Jet || NoB_2Jet)) continue;
      
      //general cuts initialization
      bool PassMPT=false;
      bool PassDHTLT=false;
      bool PassSBM=false;

      //Di-Muon mass and MET
      float Mmm=((mu1->P4()) + (mu2->P4())).M();
      MissingET *MPT = (MissingET *) branchMPT->At(0);

      //SBM calculation
      Jet *b1, *b2;
      if(DiBottom){
        b1 = (Jet *) branchJet->At(bPosition[0]);
        b2 = (Jet *) branchJet->At(bPosition[1]);
      }
      else if(Bottom_Jet){
        b1 = (Jet *) branchJet->At(bPosition[0]);
        b2 = (Jet *) branchJet->At(Jet2_nonB);
      }
      else if(NoB_2Jet){
        b1 = (Jet *) branchJet->At(Jet2_nonB);
        b2 = (Jet *) branchJet->At(Jet3_nonB);
      }
             
      float Pm1_1 = ((mu1->P4())+(b1->P4())).M();
      float Pm1_2 = ((mu2->P4())+(b2->P4())).M();
      float Pm2_1 = ((mu2->P4())+(b1->P4())).M();
      float Pm2_2 = ((mu1->P4())+(b2->P4())).M();
      float SBM1=0.;
      float SBM2=0.;
      if( fabs(Pm1_1-Pm1_2) < fabs(Pm2_1-Pm2_2) ){
        SBM1=Pm1_1;
        SBM2=Pm1_2;
      }
      else{
        SBM1=Pm2_1;
        SBM2=Pm2_2;
      }
      if(SBM1>170. || SBM2>170.) PassSBM=true;

      //DHTLT calculation
      float DHTLT=b1->PT+b2->PT-mu1->PT-mu2->PT;
      if(DHTLT<0) PassDHTLT=true;

      //new MPT cut relative to di-muon mass
      if(MPT->MET<Mmm/5.) PassMPT=true;

      // Plot their invariant mass
      if(DiBottom || Bottom_Jet){
	histMass->Fill(Mmm,weight);
	histMPT->Fill(MPT->MET,weight);
	histSimilarBranchMass->Fill(SBM1,SBM2,weight);
	histDeltaHTLT->Fill(b1->PT+b2->PT,mu1->PT+mu2->PT,weight);
	histDHTLT->Fill(DHTLT,weight);
	Mmm_DHTLT->Fill(DHTLT,Mmm,weight);

	float Dphi=((mu1->P4())+(mu2->P4())).Phi();
	Dphi-=MPT->Phi;
	Dphi=TVector2::Phi_mpi_pi(Dphi);
	float scale=MPT->MET/(((mu1->P4())+(mu2->P4())).Pt());
	histDirectionalDeviation->Fill(scale*sin(Dphi),scale*cos(Dphi),weight);

	//testing how often the muon is inside a jet radius
	/*const float dR11 = sqrt(pow(TVector2::Phi_mpi_pi(mu1->Phi-b1->Phi),2)+pow(mu1->Eta+b1->Phi,2));
	const float dR12 = sqrt(pow(TVector2::Phi_mpi_pi(mu2->Phi-b1->Phi),2)+pow(mu2->Eta+b1->Phi,2));
	float dR1 =1.;
	if(dR11<dR12) dR1=dR11; else dR1=dR12;
	const float dR21 = sqrt(pow(TVector2::Phi_mpi_pi(mu1->Phi-b2->Phi),2)+pow(mu1->Eta+b2->Phi,2));
	const float dR22 = sqrt(pow(TVector2::Phi_mpi_pi(mu2->Phi-b2->Phi),2)+pow(mu2->Eta+b2->Phi,2));
	float dR2 =1.;
	if(dR21<dR22) dR2=dR21; else dR2=dR22;
	if(Bottom_Jet){
	  Jet1_dR_1b->Fill(dR1,weight);
	  Jet2_dR_1b->Fill(dR2,weight);
	}
	else{
	  Jet1_dR_2b->Fill(dR1,weight);
	  Jet2_dR_2b->Fill(dR2,weight);
	}*/

	METvsMmm->Fill(MPT->MET/Mmm,weight);
	Mmm_METvsMmm->Fill(MPT->MET/Mmm,Mmm,weight);
	cutflow->Fill(2.,PassMPT,weight);
	cutflow->Fill(3.,PassDHTLT,weight);
	cutflow->Fill(4.,PassSBM,weight);
	if(PassSBM && PassDHTLT){
	  METvLT->Fill(scale,weight);
	  Proj->Fill(1-scale*cos(Dphi),weight);
	}

	//A:MPT, B:DHTLT, C: max(SBM) fill all permutations
	/*cutflow_ABC->Fill(1.,PassMPT,weight);
	cutflow_ACB->Fill(1.,PassMPT,weight);
	if(PassMPT){
	  cutflow_ABC->Fill(2.,PassDHTLT,weight);
	  if(PassDHTLT) cutflow_ABC->Fill(3.,PassSBM,weight);
	  cutflow_ACB->Fill(2.,PassSBM,weight);
	  if(PassSBM) cutflow_ACB->Fill(3.,PassDHTLT,weight);
	}
	cutflow_BAC->Fill(1.,PassDHTLT,weight);
	cutflow_BCA->Fill(1.,PassDHTLT,weight);
	if(PassDHTLT){
	  cutflow_BAC->Fill(2.,PassMPT,weight);
	  if(PassMPT) cutflow_BAC->Fill(3.,PassSBM,weight);
	  cutflow_BCA->Fill(2.,PassSBM,weight);
	  if(PassSBM) cutflow_BCA->Fill(3.,PassMPT,weight);
	}
	cutflow_CAB->Fill(1.,PassSBM,weight);
	cutflow_CBA->Fill(1.,PassSBM,weight);*/
	if(PassSBM){
	  cutflow_CAB->Fill(2.,PassMPT,weight);
	  if(PassMPT) cutflow_CAB->Fill(3.,PassDHTLT,weight);
	  /*cutflow_CBA->Fill(2.,PassDHTLT,weight);
	  if(PassDHTLT) cutflow_CBA->Fill(3.,PassMPT,weight);*/
	}
	if(SBM1<SBM2){
	  minSBM->Fill(SBM1,weight);
	  maxSBM->Fill(SBM2,weight);
	  Mmm_maxSBM->Fill(SBM2,Mmm,weight);
	}
	else{
	  minSBM->Fill(SBM2,weight);
	  maxSBM->Fill(SBM1,weight);
	  Mmm_maxSBM->Fill(SBM1,Mmm,weight);
	}
	if(PassMPT && PassDHTLT && PassSBM){
	  float DeltaSBM=fabs(SBM1-SBM2);
	  bool PassDSBM=false;
	  if(DeltaSBM>50.) PassDSBM=true;
	  /*if(PassDSBM)*/ massAfterCuts->Fill(Mmm,weight);
	  cutflow->Fill(5.,PassDSBM,weight);
	  DSBM->Fill(DeltaSBM,weight);
	}//PassAll
      }//OSDiMuon + DiBottom
      else if(NoB_2Jet){ //inversion control region
	if(PassMPT && PassDHTLT && PassSBM) massAfterCuts_NoB_2j->Fill(Mmm,weight);
      }
    }//mu>1
    else{
      cutflow->Fill(1.,0.,weight);
      /*cutflow_ABC->Fill(0.,0.,weight);
      cutflow_ACB->Fill(0.,0.,weight);
      cutflow_BAC->Fill(0.,0.,weight);
      cutflow_BCA->Fill(0.,0.,weight);*/
      cutflow_CAB->Fill(0.,0.,weight);
      //cutflow_CBA->Fill(0.,0.,weight);
    }
  }//event loop

  // Show resulting histograms
  TFile *savefile = new TFile(savename,"RECREATE");
  TCanvas *results=new TCanvas("results","",1800,900);
  results->SetLogy(1);
  results->Divide(6,3);
  results->cd(1);
  histJetPT->Draw();
  histJetPT->Write();
  results->cd(2);
  histBjetPT->Draw();
  histBjetPT->Write();
  results->cd(3);
  histMass->Draw();
  histMass->Write();
  results->cd(4);
  histMPT->Draw();
  histMPT->Write();
  results->cd(5);
  histSimilarBranchMass->Draw("colz");
  histSimilarBranchMass->Write();
  results->cd(6);
  histDirectionalDeviation->Draw("colz");
  histDirectionalDeviation->Write();
  results->cd(7);
  histDeltaHTLT->Draw("colz");
  histDeltaHTLT->Write();
  results->cd(8);
  histDHTLT->Draw();
  histDHTLT->Write();
  results->cd(9);
  DSBM->Draw();
  DSBM->Write();
  results->cd(10);
  maxSBM->Draw();
  maxSBM->Write();
  results->cd(11);
  minSBM->Draw();
  minSBM->Write();
  results->cd(12);
  /*cutflow_ABC->Write();
  cutflow_ACB->Write();
  cutflow_BAC->Write();
  cutflow_BCA->Write();*/
  cutflow_CAB->Write();
  /*cutflow_CBA->Write();
  cutflow_ABC->SetLineColor(1);
  cutflow_ACB->SetLineColor(2);
  cutflow_BAC->SetLineColor(3);
  cutflow_BCA->SetLineColor(4);*/
  cutflow_CAB->SetLineColor(5);
  //cutflow_CBA->SetLineColor(6);
  /*(cutflow_ABC->Draw();
  cutflow_ACB->Draw("same");
  cutflow_BAC->Draw("same");
  cutflow_BCA->Draw("same");
  cutflow_CAB->Draw("same");
  cutflow_CBA->Draw("same");*/
  cutflow_CAB->Draw();
  results->cd(13);
  histJetPT_nonB2->Write();
  histJetPT_nonB2->Draw();
  results->cd(14);
  METvLT->Write();
  METvLT->Draw();
  results->cd(15);
  Proj->Write();
  Proj->Draw();
  results->cd(16);
  METvsMmm->Write();
  METvsMmm->Draw();
  results->cd(17);
  massAfterCuts->Draw();
  massAfterCuts->Write();
  results->cd(18);
  cutflow->Draw();
  cutflow->Write();
  results->Write();
  massAfterCuts_NoB_2j->Write();
  /*Jet1_dR_1b->Write();
  Jet2_dR_1b->Write();
  Jet1_dR_2b->Write();
  Jet2_dR_2b->Write();*/
  Mmm_METvsMmm->Write();
  Mmm_maxSBM->Write();
  Mmm_DHTLT->Write();
  savefile->Close();
}

