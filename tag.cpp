#define bTag_AntiKt4EMPFlowJets_cxx
#include "bTag_AntiKt4EMPFlowJets.h"

#include </usr/local/root/include/TH1.h>
#include </usr/local/root/include/TH2.h>

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

//extern int thfff;
int thfff = -1;

void bTag_AntiKt4EMPFlowJets::Loop()
{
  if (fChain == 0) return;
   // book histograms
  TH1* h_dz = new TH1F("dz","",50,-0.5,0.5);
  TH1* h_njet = new TH1F("njet","",20,0.,20.);
  const int ntag = 2;
    const char* chtag[ntag] = {
      "mv2c10", "ip3d"
    };
    double xmin[ntag] = {-1.,-20.};
    double xmax[ntag] = { 1., 50.};

    const int nf = 3; // l,c,b
    TH1F* h_w[ntag*nf];
    for (int jf = 0; jf<nf; ++jf) {
      for (int itag = 0; itag<ntag; ++itag) {
  	TString sweight = "w_"; sweight+=chtag[itag];
  	sweight+="_"; sweight+=jf;
  	h_w[itag*nf+jf] = new TH1F(sweight,"weight",2000,xmin[itag],xmax[itag]);
      }
    }
     TH1* h_all_dl1 = new TH1F("all_dl1","",nf,0.,nf);
     TH1* h_tag_dl1 = (TH1*)h_all_dl1->Clone("tag_dl1");
     TH1* h_all_mv2c10 = (TH1*)h_all_dl1->Clone("all_mv2c10");
     TH1* h_tag_mv2c10 = (TH1*)h_all_dl1->Clone("tag_mv2c10");
     TH1* h_ratio = (TH1*)h_all_dl1->Clone("ratio");

     vector<double> wtag;

     Long64_t nentries = fChain->GetEntriesFast();

     for (Long64_t jentry=0; jentry<nentries; jentry++) {
       Long64_t ientry = LoadTree(jentry);
       if (ientry < 0) break;
       if (jentry%1000==0) cout << jentry << '\r'; cout.flush();
       fChain->GetEntry(jentry);
       //std::cout<<"I am here"<<std::endl;
       // skip events with badly reconstructed PV
       double dz = PVz-truth_PVz;
       h_dz->Fill(dz);
       //if (fabs(dz)>0.1) continue;

           int njet = 0; if (jet_pt) njet = jet_pt->size();
           h_njet->Fill(njet);
           if (njet==0) continue;

           // fill tag weights
           wtag.clear();
           for (int ijet = 0; ijet<njet; ++ijet) {
             wtag.push_back((*jet_mv2c10)[ijet]);
             wtag.push_back((*jet_ip3d_llr)[ijet]);
             }

	     	 for (int ijet = 0; ijet<njet; ++ijet) {
	     	   // select good jets
	     	   if ((*jet_isPU)[ijet]==1) continue;
	     	   if ((*jet_truthMatch)[ijet]==0) continue;
	     	   if ((*jet_pt)[ijet]<120e3 && (*jet_JVT)[ijet]<0.59) continue;
	     	   if ((*jet_aliveAfterOR)[ijet]==0) continue;
	     	   if ((*jet_aliveAfterORmu)[ijet]==0) continue;
	     	   if ((*jet_pt)[ijet]<20e3) continue;
	     	   if (fabs((*jet_eta)[ijet])>2.5) continue;
	     	   // figure out their flavor
	     	   int jf = (*jet_LabDr_HadF)[ijet];
		   // don't forget to skip taus!!!
      if (jf==4) jf = 1; else if (jf==5) jf = 2; else if (jf!=0) continue;
      // fill performance histograms
      for (int itag = 0; itag<ntag; ++itag) {
        h_w[itag*nf+jf]->Fill(wtag[ijet*ntag+itag]);
      }
      h_all_mv2c10->Fill(jf);
      if (wtag[ijet*ntag]>0.83) h_tag_mv2c10->Fill(jf);
      h_ratio->Divide(h_tag_mv2c10,h_all_mv2c10);
    }
  }

  std::cout << "mv2c10 b-tagging efficiency: " << h_tag_mv2c10->GetBinContent(3)/h_all_mv2c10->GetBinContent(3) << std::endl;
  std::cout << "mv2c10 c rejection: " << h_all_mv2c10->GetBinContent(2)/h_tag_mv2c10->GetBinContent(2) << std::endl;
  std::cout << "mv2c10 light rejection: " << h_all_mv2c10->GetBinContent(1)/h_tag_mv2c10->GetBinContent(1) << std::endl;
}

void tag()
{
  const char* chtree = "bTag_AntiKt4EMPFlowJets";
  TChain* tt = new TChain(chtree);
  tt->Add("flav_Akt4EMPf.root");
  TTree* mytree = (TTree*)gROOT->FindObject(chtree);
  Long64_t nentries = mytree->GetEntries();
  bTag_AntiKt4EMPFlowJets* t = new bTag_AntiKt4EMPFlowJets(mytree); // all branches are set
  TFile* fout = new TFile("a1.root","RECREATE");
  t->Loop();
  fout->Write();
  delete fout;
}
