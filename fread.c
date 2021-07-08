using namespace std;
void fread(){
    gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
    const char* chtree1 = "bTag_AntiKt4EMPFlowJets";
    TChain* tt1 = new TChain(chtree1);
    tt1->Add("flav_Akt4EMPf75.root");
    TTree* tree1 = (TTree*)gROOT->FindObject(chtree1);
    vector<vector<float>>*  jet_bH_Lxy = 0;
    tree1->SetBranchAddress("jet_bH_Lxy",&jet_bH_Lxy);
    vector<float>* jet_pt = 0;
    tree1->SetBranchAddress("jet_pt",&jet_pt);
    
    Long64_t nentries = tree1->GetEntriesFast();
    for (Int_t i=0; i<nentries; i++) {
        int njet1 = 0; if (jet_pt!=0) njet1 = jet_pt->size();
           if (njet1==0) continue;
           for (int ijet = 0; ijet<njet1; ++ijet) {
                float Lxy = (*jet_bH_Lxy)[i][0];
                cout<<Lxy<<" ";
        
    }
    cout<<"\n";
  }
}
