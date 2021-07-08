void w()
{
  // load the tree

  const char* chtree = "bTag_AntiKt4EMPFlowJets";
  TChain* tt = new TChain(chtree);
  tt->Add("flav_Akt4EMPf.root");
  TTree* tree = (TTree*)gROOT->FindObject(chtree);

  // setup useful branches

  /**
  TString stag = "mv2c10"; double xmin = -1., xmax = 1.;
  vector<double>* jet_tagweight = 0; //double for mv2c10
  tree->SetBranchAddress("jet_mv2c10",&jet_tagweight);
  */
  
  TString stag = "ip3d_llr"; double xmin = -15., xmax = 35.;
  vector<float>* jet_tagweight = 0; // float for jet_ip3d_llr 
  tree->SetBranchAddress("jet_ip3d_llr",&jet_tagweight);
  
  vector<int>* jet_LabDr_HadF = 0;
  tree->SetBranchAddress("jet_LabDr_HadF",&jet_LabDr_HadF);
  vector<int>* jet_isPU = 0;
  tree->SetBranchAddress("jet_isPU",&jet_isPU);
  vector<int>* jet_truthMatch =0;
  vector<float>* jet_pt=0;
  tree->SetBranchAddress("jet_pt", &jet_pt);
  vector<int>*  jet_aliveAfterOR=0;
  tree->SetBranchAddress("jet_aliveAfterOR", &jet_aliveAfterOR);
  vector<int>*  jet_aliveAfterORmu=0;
  tree->SetBranchAddress("jet_aliveAfterORmu", &jet_aliveAfterORmu);
  vector<float>* jet_eta=0;
  tree->SetBranchAddress("jet_eta", &jet_eta);
  vector<float>* jet_JVT=0;
  tree->SetBranchAddress("jet_JVT", &jet_JVT);
  
  Float_t         truth_PVz;
  Float_t         PVz;
  
  tree->SetBranchAddress(TString("jet_")+stag,&jet_tagweight);

  // book histograms

  TFile* fout = new TFile(TString("a1_"+stag+"_reg.root"),"RECREATE");
  const int nfl = 3; // light, c, b
  TH1* h_w[nfl];
  h_w[0] = new TH1F(stag+"_0","",100,xmin,xmax);
  h_w[1] = (TH1*)h_w[0]->Clone(stag+"_1");
  h_w[2] = (TH1*)h_w[0]->Clone(stag+"_2");
  TH1* h_dz = new TH1F("dz","",50,-0.5,0.5);
  TH1* h_njet = new TH1F("njet","",20,0.,20.);
  // loop over the tree entries

  Long64_t nentries = tree->GetEntriesFast();
  for (Long64_t jentry=0; jentry<nentries; jentry++) {
    if (tree->LoadTree(jentry)<0) break;
    tree->GetEntry(jentry);
    if (jet_tagweight==0) continue;
    double dz = PVz-truth_PVz;
    h_dz->Fill(dz);
    
    int njet = 0; if (jet_pt) njet = jet_pt->size();
           h_njet->Fill(njet);
           if (njet==0) continue;

    for (int ijet = 0; ijet<njet; ++ijet) {
	     	   // select good jets
	     	   if ((*jet_isPU)[ijet]==1) continue;
	     	   //if ((*jet_truthMatch)[ijet]==0) continue;
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
      h_w[jf]->Fill((*jet_tagweight)[ijet]);
      
    }
  }
  fout->Write();
  delete fout;
}
 
