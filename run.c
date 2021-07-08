using namespace std;
void run()
{ 
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");  
  TCanvas *c1 = new TCanvas("c1","c1",1200,800) ;
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  gStyle->SetTitleX(0.1f);
  gStyle->SetTitleW(0.8f);
  //c1->cd();
   TLegend* l = new TLegend(0.60,.90,0.85,0.95);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(42);
    //l->SetMargin(0.4);
  
  float frac_c = 0.080;
  
  TH1F* h_Lxy_all = new TH1F("Lxy_all","B-Tag Eff vs Lxy-5",25,0,50);
  TH1F* h_Lxy_tag = (TH1F*)h_Lxy_all->Clone("Lxy_tag");
  TH1F* h_ratio = (TH1F*)h_Lxy_all->Clone("ratio");
  TH1F* h_dz = new TH1F("dz","",50,-0.5,0.5);
  TH1F* h_njet = new TH1F("njet","",20,0.,20.);
  
  TH1F* h_Lxy_all1 = new TH1F("Lxy_all1","B-Tag Eff vs Lxy-7.5",25,0,50);
  TH1F* h_Lxy_tag1 = (TH1F*)h_Lxy_all1->Clone("Lxy_tag1");
  TH1F* h_ratio1 = (TH1F*)h_Lxy_tag1->Clone("ratio1");
  TH1F* h_dz1 = new TH1F("dz1","",50,-0.5,0.5);
  TH1F* h_njet1 = new TH1F("njet1","",20,0.,20.);
 
  /**
  TH1F* h_Lxy_all2 = new TH1F("Lxy_all2","B-Tag Eff vs Lxy-7.5",25,0,50);
  TH1F* h_Lxy_tag2 = (TH1F*)h_Lxy_all2->Clone("Lxy_tag2");
  TH1F* h_ratio2 = (TH1F*)h_Lxy_tag2->Clone("ratio2");
  TH1F* h_dz2 = new TH1F("dz1","",50,-0.5,0.5);
  TH1F* h_njet2 = new TH1F("njet1","",20,0.,20.);
  */
  int ctr=0;
  ctr=0;
  if (ctr==0){ //unmodified IP file
  // load the tree
  const char* chtree = "bTag_AntiKt4EMPFlowJets";
  TChain* tt = new TChain(chtree);
  tt->Add("flav_Akt4EMPf.root");
  TTree* tree = (TTree*)gROOT->FindObject(chtree);

  // setup useful branches
  vector<double>* jet_tagweight = 0; //for mv2c10 or mv2c10rnn
  //vector<float>* jet_tagweight = 0; //for ip3dllr
  tree->SetBranchAddress("jet_mv2c10",&jet_tagweight);
  //tree->SetBranchAddress("jet_ip3d_llr",&jet_tagweight);
  //tree->SetBranchAddress("jet_mv2c10rnn", &jet_tagweight);
  vector<int>* jet_isPU = 0;
  tree->SetBranchAddress("jet_isPU",&jet_isPU);
  vector<vector<float>>*  jet_bH_Lxy = 0;
  tree->SetBranchAddress("jet_bH_Lxy",&jet_bH_Lxy);
  vector<float>* jet_pt = 0;
  tree->SetBranchAddress("jet_pt",&jet_pt);
  vector<int>* jet_LabDr_HadF = 0;
  tree->SetBranchAddress("jet_LabDr_HadF",&jet_LabDr_HadF);
  vector<int>* jet_truthMatch =0;
  tree->SetBranchAddress("jet_truthMatch", &jet_truthMatch);
  vector<int>*  jet_aliveAfterOR=0;
  tree->SetBranchAddress("jet_aliveAfterOR", &jet_aliveAfterOR);
  vector<int>*  jet_aliveAfterORmu=0;
  tree->SetBranchAddress("jet_aliveAfterORmu", &jet_aliveAfterORmu);
  vector<float>* jet_eta=0;
  tree->SetBranchAddress("jet_eta", &jet_eta);
  vector<float>* jet_JVT=0;
  tree->SetBranchAddress("jet_JVT", &jet_JVT);
  vector<double>* jet_dl1_pb = 0;
  vector<double>* jet_dl1_pu = 0;
  vector<double>* jet_dl1_pc = 0;
  tree->SetBranchAddress("jet_dl1_pb", &jet_dl1_pb);
  tree->SetBranchAddress("jet_dl1_pc", &jet_dl1_pc);
  tree->SetBranchAddress("jet_dl1_pu", &jet_dl1_pu);
  Float_t         truth_PVz;
  Float_t         PVz;

  //TFile* fout1 = new TFile("Lxy_threehist.root","RECREATE");

  // loop over the tree entries
  Long64_t nentries = tree->GetEntriesFast();
  for (Long64_t jentry=0; jentry<nentries; jentry++) {
    if (tree->LoadTree(jentry)<0) break;
    tree->GetEntry(jentry);
    if (jet_tagweight==0) continue;
    double dz = PVz-truth_PVz;
    h_dz->Fill(dz);
    // check if there are jets
    int njet = 0; if (jet_pt) njet = jet_pt->size();
           h_njet->Fill(njet);
           if (njet==0) continue;

    for (int ijet = 0; ijet<njet; ++ijet) {
	     	   // select good jets
	     	   if ((*jet_isPU)[ijet]==1) continue;
	     	   if ((*jet_truthMatch)[ijet]==0) continue;
	     	   if ((*jet_pt)[ijet]<120e3 && (*jet_JVT)[ijet]<0.59) continue;
	     	   if ((*jet_aliveAfterOR)[ijet]==0) continue;
	     	   if ((*jet_aliveAfterORmu)[ijet]==0) continue;
	     	   if ((*jet_pt)[ijet]==0) continue;
	     	   if (fabs((*jet_eta)[ijet])>2.5) continue;
	     	   // figure out their flavor
	     	   int jf = (*jet_LabDr_HadF)[ijet];
		   // don't forget to skip taus!!!
      if (jf==5) { jf = 2;float Lxy=0;
                   if (!(*jet_bH_Lxy)[ijet].empty()){Lxy = fabs((*jet_bH_Lxy)[ijet][0]);} 
                   h_Lxy_all->Fill(Lxy);
		  if ((*jet_tagweight)[ijet]>0.83) {//0.83 for mv2c10, 0.87 for mv2c10rnn, 
		      h_Lxy_tag->Fill(Lxy);
            }
        /**
          float a = (*jet_dl1_pb)[ijet];
          float b = frac_c*(*jet_dl1_pc)[ijet] + (1-frac_c)*(*jet_dl1_pu)[ijet];
          float w_dl1 = -999; if (a>0 && b>0) w_dl1 = log(a/b);
          if (w_dl1>0) { h_Lxy_tag->Fill(Lxy); }//2.02 for dl1pb
          */
       }
    }
  }
  //fout1->Write();
  //delete fout1;
  //delete tree;
  //delete chtree;
  delete tt;
}
 
  ctr=1;             
  if (ctr==1){ //modified IP file
  // load the tree
  const char* chtree1 = "bTag_AntiKt4EMPFlowJets";
  TChain* tt1 = new TChain(chtree1);
  tt1->Add("flav_Akt4EMPf75.root");
  TTree* tree1 = (TTree*)gROOT->FindObject(chtree1);

  // setup useful branches
  vector<double>* jet_tagweight = 0; //for mv2c10
  //vector<float>* jet_tagweight = 0; //for ip3dllr
  tree1->SetBranchAddress("jet_mv2c10",&jet_tagweight);
  //tree1->SetBranchAddress("jet_ip3d_llr",&jet_tagweight);
  //tree1->SetBranchAddress("jet_mv2c10rnn", &jet_tagweight);
  vector<int>* jet_isPU = 0;
  tree1->SetBranchAddress("jet_isPU",&jet_isPU);
  vector<vector<float>>*  jet_bH_Lxy = 0;
  tree1->SetBranchAddress("jet_bH_Lxy",&jet_bH_Lxy);
  vector<float>* jet_pt = 0;
  tree1->SetBranchAddress("jet_pt",&jet_pt);
  vector<int>* jet_LabDr_HadF = 0;
  tree1->SetBranchAddress("jet_LabDr_HadF",&jet_LabDr_HadF);
  vector<int>* jet_truthMatch =0;
  tree1->SetBranchAddress("jet_truthMatch", &jet_truthMatch);
  vector<int>*  jet_aliveAfterOR=0;
  tree1->SetBranchAddress("jet_aliveAfterOR", &jet_aliveAfterOR);
  vector<int>*  jet_aliveAfterORmu=0;
  tree1->SetBranchAddress("jet_aliveAfterORmu", &jet_aliveAfterORmu);
  vector<float>* jet_eta=0;
  tree1->SetBranchAddress("jet_eta", &jet_eta);
  vector<float>* jet_JVT=0;
  tree1->SetBranchAddress("jet_JVT", &jet_JVT);
  vector<double>* jet_dl1_pb = 0;
  vector<double>* jet_dl1_pu = 0;
  vector<double>* jet_dl1_pc = 0;
  tree1->SetBranchAddress("jet_dl1_pb", &jet_dl1_pb);
  tree1->SetBranchAddress("jet_dl1_pc", &jet_dl1_pc);
  tree1->SetBranchAddress("jet_dl1_pu", &jet_dl1_pu);
  
  Float_t         truth_PVz;
  Float_t         PVz;

  //TFile* fout1 = new TFile("Lxy_threehist1.root","RECREATE");

  // loop over the tree entries
  Long64_t nentries = tree1->GetEntriesFast();
  for (Long64_t jentry=0; jentry<nentries; jentry++) {
    if (tree1->LoadTree(jentry)<0) break;
    tree1->GetEntry(jentry);
    if (jet_tagweight==0) continue;
    double dz1 = PVz-truth_PVz;
    h_dz1->Fill(dz1);
    // check if there are jets
    int njet1 = 0; if (jet_pt) njet1 = jet_pt->size();
           h_njet1->Fill(njet1);
           if (njet1==0) continue;

    for (int ijet = 0; ijet<njet1; ++ijet) {
	     	   // select good jets
	     	   if ((*jet_isPU)[ijet]==1) continue;
	     	   if ((*jet_truthMatch)[ijet]==0) continue;
	     	   if ((*jet_pt)[ijet]<120e3 && (*jet_JVT)[ijet]<0.59) continue;
	     	   if ((*jet_aliveAfterOR)[ijet]==0) continue;
	     	   if ((*jet_aliveAfterORmu)[ijet]==0) continue;
	     	   if ((*jet_pt)[ijet]==0) continue;
	     	   if (fabs((*jet_eta)[ijet])>2.5) continue;
	     	   // figure out their flavor
	     	   int jf = (*jet_LabDr_HadF)[ijet];
		   // don't forget to skip taus!!!
      if (jf==5) { jf = 2;float Lxy=0;
                   if (!(*jet_bH_Lxy)[ijet].empty()) {Lxy = fabs((*jet_bH_Lxy)[ijet][0]);} 
                   h_Lxy_all1->Fill(Lxy);
		   /**
                   if ((*jet_tagweight)[ijet]>0.83) {//0.83 for mv2c10, 0.87 for mv2c10rnn, 
		      h_Lxy_tag1->Fill(Lxy);
            }
              */
          float a = (*jet_dl1_pb)[ijet];
          float b = frac_c*(*jet_dl1_pc)[ijet] + (1-frac_c)*(*jet_dl1_pu)[ijet];
          float w_dl1 = -999; if (a>0 && b>0) w_dl1 = log(a/b);
          if (w_dl1>0) {h_Lxy_tag1->Fill(Lxy);}//2.02 for dl1
	   		  
      }
    }
  }
  //fout1->Write();
  //delete fout1;
  //delete tree1;
  //delete chtree1;
   delete tt1;
  
  }
 
/** 
  h_Lxy_all->SetLineColor(2);
  h_Lxy_all->SetLineWidth(2);
  h_Lxy_all->SetFillColor(48);
  h_Lxy_all->SetFillStyle(3001);
  h_Lxy_all->SetMarkerStyle(20);
  h_Lxy_all->Draw("E HIST");
  
  h_Lxy_tag->SetLineColor(3);
  h_Lxy_tag->SetLineWidth(2);
  h_Lxy_tag->SetFillColor(48);
  h_Lxy_tag->SetFillStyle(3001);
  h_Lxy_tag->SetMarkerStyle(20);
  h_Lxy_tag->Draw("E HIST SAME");
 
  h_Lxy_all1->SetLineColor(2);
  h_Lxy_all1->SetLineWidth(2);
  h_Lxy_all1->SetFillColor(48);
  h_Lxy_all1->SetFillStyle(3001);
  h_Lxy_all1->SetMarkerStyle(20);
  h_Lxy_all1->Draw("SAME");
  
  h_Lxy_tag1->SetLineColor(3);
  h_Lxy_tag1->SetLineWidth(2);
  h_Lxy_tag1->SetFillColor(48);
  h_Lxy_tag1->SetFillStyle(3001);
  h_Lxy_tag1->SetMarkerStyle(20);
  h_Lxy_tag1->Draw("SAME");
*/ 
 
  h_ratio->Divide(h_Lxy_tag,h_Lxy_all,1,1,"B");
  // h_ratio->Divide(h_Lxy_tag,h_Lxy_all);
  h_ratio->SetLineColor(1);
  h_ratio->SetLineWidth(3);
  //h_ratio->SetFillColor(48);
  //h_ratio->SetFillStyle(3001);
  h_ratio->SetMarkerStyle(20);
  
  //h_ratio->SetTitle("B-Tag Eff vs Lxy(MV2c10:0.83);Lxy(mm);Efficiency");
  //h_ratio->SetTitle("B-Tag Eff vs Lxy(IP3D);Lxy(mm);Efficiency");
  //h_ratio->SetTitle("B-Tag Eff vs Lxy(MV2c10rnn:0.87);Lxy(mm);Efficiency");
  h_ratio->SetTitle("B-Tag Eff vs Lxy(DL1:2.02);Lxy(mm);Efficiency");
  
  //l->AddEntry(h_ratio,"d0max=5.0mm, z0max=5.0mm","l");
  
  h_ratio1->Divide(h_Lxy_tag1,h_Lxy_all1,1,1,"B");
  double ymax = 1.1*TMath::Max(h_ratio->GetMaximum(),h_ratio1->GetMaximum());
  
  h_ratio->SetMaximum(ymax);
  h_ratio->Draw("HIST E0");
  
  h_ratio1->SetLineColor(2);
  h_ratio1->SetLineWidth(2);
  //h_ratio1->SetFillColor(54);
  //h_ratio1->SetFillStyle(3001);
  h_ratio1->SetMarkerStyle(20);
  //h_ratio1->SetTitle("B-Tag Eff vs Lxy(IP3D/MV2c10);Lxy(mm);Efficiency");
  h_ratio1->Draw("HIST E0 SAME");
  //h_ratio1->Draw("SAME");
    
  l->AddEntry(h_ratio,"d0max=5.0mm, z0max=5.0mm","l");
  l->AddEntry(h_ratio1,"d0max=7.5mm, z0max=7.5mm","l");
  l->Draw();

  c1->Modified();
  c1->Update();
  //c1->SaveAs("TaggedvsAll_PU_MV2c10_1k.png");
  //c1->SaveAs("TaggedvsAll_PU_IP3D.png");
  //c1->SaveAs("TaggedvsAll_PU_Mv2c10rnn.png");
  c1->SaveAs("TaggedvsAll_PU_DL1_1k.png");
}
