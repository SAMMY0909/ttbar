// running event loop in a single macro without Makefile

void r()
{
  // open input file and get a handle to the tree
  TFile* ff = new TFile("flav_Akt4EMPf75.root");
  TTree* tree = (TTree*)ff->Get("bTag_AntiKt4EMPFlowJets");

  // set branches
  float PVz, truth_PVz;
  vector<float>* jet_pt = 0;
  vector<float>* jet_eta = 0;
  vector<int>* jet_isPU = 0;
  vector<int>* jet_truthMatch = 0;
  vector<int>* jet_aliveAfterOR = 0;
  vector<int>* jet_aliveAfterORmu = 0;
  vector<int>* jet_LabDr_HadF = 0;
  vector<float>* jet_JVT = 0;
  vector<double>* jet_mv2c10 = 0;
  vector<float>* jet_ip3d_llr = 0;
  vector<float>* jet_ip2d_llr = 0;
  vector<float>* sv1_llr = 0;
  vector<double>* jet_dl1_pb = 0;
  vector<double>* jet_dl1_pu = 0;
  vector<double>* jet_dl1_pc = 0;
  vector<vector<float> >* jet_bH_Lxy = 0;

  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");

  tree->SetBranchAddress("PVz", &PVz);
  tree->SetBranchAddress("truth_PVz", &truth_PVz);
  tree->SetBranchAddress("jet_pt", &jet_pt);
  tree->SetBranchAddress("jet_eta", &jet_eta);
  tree->SetBranchAddress("jet_isPU", &jet_isPU);
  tree->SetBranchAddress("jet_truthMatch", &jet_truthMatch);
  tree->SetBranchAddress("jet_aliveAfterOR", &jet_aliveAfterOR);
  tree->SetBranchAddress("jet_aliveAfterORmu", &jet_aliveAfterORmu);
  tree->SetBranchAddress("jet_LabDr_HadF", &jet_LabDr_HadF);
  tree->SetBranchAddress("jet_JVT", &jet_JVT);
  tree->SetBranchAddress("jet_mv2c10", &jet_mv2c10);
  tree->SetBranchAddress("jet_ip3d_llr", &jet_ip3d_llr);
  tree->SetBranchAddress("jet_ip2d_llr", &jet_ip2d_llr);
  tree->SetBranchAddress("sv1_llr", &sv1_llr);
  tree->SetBranchAddress("jet_dl1_pb", &jet_dl1_pb);
  tree->SetBranchAddress("jet_dl1_pc", &jet_dl1_pc);
  tree->SetBranchAddress("jet_dl1_pu", &jet_dl1_pu);
  tree->SetBranchAddress("jet_bH_Lxy", &jet_bH_Lxy);

  // open output file
  TFile* ff_out = new TFile("a1.root","RECREATE");

  // book histograms
  TH1* h_dz = new TH1F("dz","",50,-0.5,0.5);
  TH1* h_njet = new TH1F("njet","",20,0.,20.);

  const int ntag = 5;
  const char* chtag[ntag] = {
    "mv2c10", "ip3d", "ip2d", "sv1", "dl1"
  };
  double xmin[ntag] = {-1.,-20.,-20.,-5.,-25.};
  double xmax[ntag] = { 1., 50., 50.,15., 25.};

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

  TH1* h_all_mv2c10_lxy[nf];
  TH1* h_tag_mv2c10_lxy[nf];
  TH1* h_all_dl1_lxy[nf];
  TH1* h_tag_dl1_lxy[nf];
  for (int jf = 0; jf<nf; ++jf) {
    h_all_mv2c10_lxy[jf] = new TH1F(TString("all_mv2c10_")+jf,"",50,0.,50.);
    h_tag_mv2c10_lxy[jf] = new TH1F(TString("tag_mv2c10_")+jf,"",50,0.,50.);
    h_all_dl1_lxy[jf] = new TH1F(TString("all_dl1_")+jf,"",50,0.,50.);
    h_tag_dl1_lxy[jf] = new TH1F(TString("tag_dl1_")+jf,"",50,0.,50.);
  }

  vector<double> wtag;

  Long64_t nentries = tree->GetEntries();

  for (Long64_t jentry=0; jentry<nentries; jentry++) {
    tree->GetEntry(jentry);
    if (jentry%1000==0) cout << jentry << '\r'; cout.flush();

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
      wtag.push_back((*jet_ip2d_llr)[ijet]);
      wtag.push_back((*sv1_llr)[ijet]);
      double frac_c = 0.080;
      double a = (*jet_dl1_pb)[ijet];
      double b = frac_c*(*jet_dl1_pc)[ijet] + (1-frac_c)*(*jet_dl1_pu)[ijet];
      double w_dl1 = -999; if (a>0 && b>0) w_dl1 = log(a/b);
      wtag.push_back(w_dl1);
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
      h_all_dl1->Fill(jf);
      if (wtag[ijet*ntag+ntag-1]>2.01) h_tag_dl1->Fill(jf);
      
// find B-hadron decay length
      double lxy = 0; if (!(*jet_bH_Lxy)[ijet].empty()) lxy = (*jet_bH_Lxy)[ijet][0];
      h_all_mv2c10_lxy[jf]->Fill(lxy);
      if (wtag[ijet*ntag]>0.83) h_tag_mv2c10_lxy[jf]->Fill(lxy);
      h_all_dl1_lxy[jf]->Fill(lxy);
      if (wtag[ijet*ntag+ntag-1]>2.01) h_tag_dl1_lxy[jf]->Fill(lxy);
      // fill performance vs Lxy
    }
  }

  cout << "mv2c10 b-tagging efficiency: " << h_tag_mv2c10->GetBinContent(3)/h_all_mv2c10->GetBinContent(3) << endl;
  cout << "mv2c10 c rejection: " << h_all_mv2c10->GetBinContent(2)/h_tag_mv2c10->GetBinContent(2) << endl;
  cout << "mv2c10 light rejection: " << h_all_mv2c10->GetBinContent(1)/h_tag_mv2c10->GetBinContent(1) << endl;

  cout << "dl1 b-tagging efficiency: " << h_tag_dl1->GetBinContent(3)/h_all_dl1->GetBinContent(3) << endl;
  cout << "dl1 c rejection: " << h_all_dl1->GetBinContent(2)/h_tag_dl1->GetBinContent(2) << endl;
  cout << "dl1 light rejection: " << h_all_dl1->GetBinContent(1)/h_tag_dl1->GetBinContent(1) << endl;

  ff_out->Write();
  delete ff_out;
}
