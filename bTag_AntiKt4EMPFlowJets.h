//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Mar 21 17:18:20 2021 by ROOT version 6.23/01
// from TChain bTag_AntiKt4EMPFlowJets/
//////////////////////////////////////////////////////////

#ifndef bTag_AntiKt4EMPFlowJets_h
#define bTag_AntiKt4EMPFlowJets_h

#include </usr/local/root/include/TROOT.h>
#include </usr/local/root/include/TChain.h>
#include </usr/local/root/include/TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
using std::vector;

class bTag_AntiKt4EMPFlowJets {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           runnb;
   Int_t           eventnb;
   Int_t           mcchan;
   Float_t         mcwg;
   Float_t         avgmu;
   Int_t           actmu;
   Float_t         PVx;
   Float_t         PVy;
   Float_t         PVz;
   Float_t         truth_PVx;
   Float_t         truth_PVy;
   Float_t         truth_PVz;
   Int_t           njets;
   vector<float>   *jet_pt;
   vector<float>   *jet_eta;
   vector<float>   *jet_phi;
   vector<float>   *jet_E;
   vector<float>   *jet_pt_orig;
   vector<float>   *jet_eta_orig;
   vector<float>   *jet_phi_orig;
   vector<float>   *jet_E_orig;
   vector<int>     *jet_LabDr_HadF;
   vector<int>     *jet_DoubleHadLabel;
   vector<float>   *jet_JVT;
   vector<float>   *jet_m;
   vector<float>   *jet_nConst;
   vector<float>   *jet_dRiso;
   vector<int>     *jet_truthMatch;
   vector<int>     *jet_isPU;
   vector<int>     *jet_aliveAfterOR;
   vector<int>     *jet_aliveAfterORmu;
   vector<int>     *jet_isBadMedium;
   vector<float>   *jet_truthPt;
   vector<float>   *jet_dRminToB;
   vector<float>   *jet_dRminToC;
   vector<float>   *jet_dRminToT;
   vector<int>     *jet_truthMatch_sm75GeV;
   vector<double>  *jet_dl1_pb;
   vector<double>  *jet_dl1_pc;
   vector<double>  *jet_dl1_pu;
   vector<double>  *jet_dl1rmu_pb;
   vector<double>  *jet_dl1rmu_pc;
   vector<double>  *jet_dl1rmu_pu;
   vector<double>  *jet_dl1r_pb;
   vector<double>  *jet_dl1r_pc;
   vector<double>  *jet_dl1r_pu;
   vector<double>  *jet_mv2c10;
   vector<double>  *jet_mv2c10mu;
   vector<double>  *jet_mv2c10rnn;
   vector<double>  *jet_mv2c100;
   vector<double>  *jet_mv2cl100;
   vector<double>  *jet_mv2c10Flip;
   vector<double>  *jet_MV2rnnFlip;
   vector<double>  *jet_MV2r;
   vector<double>  *jet_DL1Flip_pb;
   vector<double>  *jet_DL1Flip_pc;
   vector<double>  *jet_DL1Flip_pu;
   vector<double>  *jet_mv2c100Flip;
   vector<double>  *jet_mv2cl100Flip;
   vector<double>  *jet_dl1rFlip_pb;
   vector<double>  *jet_dl1rFlip_pc;
   vector<double>  *jet_dl1rFlip_pu;
   vector<double>  *jet_MV2rFlip;
   vector<float>   *jet_ip2d_pb;
   vector<float>   *jet_ip2d_pc;
   vector<float>   *jet_ip2d_pu;
   vector<float>   *jet_ip2d_llr;
   vector<float>   *jet_ip3d_pb;
   vector<float>   *jet_ip3d_pc;
   vector<float>   *jet_ip3d_pu;
   vector<float>   *jet_ip3d_llr;
   vector<float>   *jet_ip2;
   vector<float>   *jet_ip2_c;
   vector<float>   *jet_ip2_cu;
   vector<float>   *jet_ip2_nan;
   vector<float>   *jet_ip2_c_nan;
   vector<float>   *jet_ip2_cu_nan;
   vector<float>   *jet_ip3;
   vector<float>   *jet_ip3_c;
   vector<float>   *jet_ip3_cu;
   vector<float>   *jet_ip3_nan;
   vector<float>   *jet_ip3_c_nan;
   vector<float>   *jet_ip3_cu_nan;
   vector<float>   *jet_rnnip_pb;
   vector<float>   *jet_rnnip_pc;
   vector<float>   *jet_rnnip_pu;
   vector<float>   *jet_rnnip_ptau;
   vector<float>   *jet_ip3dneg_pb;
   vector<float>   *jet_ip3dneg_pc;
   vector<float>   *jet_ip3dneg_pu;
   vector<float>   *jet_ip3dneg_llr;
   vector<float>   *jet_ip3neg;
   vector<float>   *jet_ip3neg_c;
   vector<float>   *jet_ip3neg_cu;
   vector<float>   *jet_ip3neg_nan;
   vector<float>   *jet_ip3neg_c_nan;
   vector<float>   *jet_ip3neg_cu_nan;
   vector<float>   *jet_ip3dflip_pb;
   vector<float>   *jet_ip3dflip_pc;
   vector<float>   *jet_ip3dflip_pu;
   vector<float>   *jet_ip3dflip_llr;
   vector<float>   *jet_ip3flip;
   vector<float>   *jet_ip3flip_c;
   vector<float>   *jet_ip3flip_cu;
   vector<float>   *jet_ip3flip_nan;
   vector<float>   *jet_ip3flip_c_nan;
   vector<float>   *jet_ip3flip_cu_nan;
   vector<int>     *jet_sv1_Nvtx;
   vector<float>   *jet_sv1_ntrkv;
   vector<float>   *jet_sv1_n2t;
   vector<float>   *jet_sv1_m;
   vector<float>   *jet_sv1_efc;
   vector<float>   *jet_sv1_sig3d;
   vector<float>   *sv1_llr;
   vector<float>   *jet_sv1_normdist;
   vector<float>   *jet_sv1_deltaR;
   vector<float>   *jet_sv1_Lxy;
   vector<float>   *jet_sv1_L3d;
   vector<vector<float> > *jet_sv1_vtx_x;
   vector<vector<float> > *jet_sv1_vtx_y;
   vector<vector<float> > *jet_sv1_vtx_z;
   vector<int>     *jet_nBHadr;
   vector<int>     *jet_nCHadr;
   vector<vector<int> > *jet_bH_pdgId;
   vector<vector<int> > *jet_bH_parent_pdgId;
   vector<vector<float> > *jet_bH_pt;
   vector<vector<float> > *jet_bH_eta;
   vector<vector<float> > *jet_bH_phi;
   vector<vector<float> > *jet_bH_E;
   vector<vector<float> > *jet_bH_charge;
   vector<vector<float> > *jet_bH_Lxy;
   vector<vector<float> > *jet_bH_x;
   vector<vector<float> > *jet_bH_y;
   vector<vector<float> > *jet_bH_z;
   vector<vector<float> > *jet_bH_dRjet;
   vector<vector<int> > *jet_cH_pdgId;
   vector<vector<int> > *jet_cH_parent_pdgId;
   vector<vector<float> > *jet_cH_pt;
   vector<vector<float> > *jet_cH_eta;
   vector<vector<float> > *jet_cH_phi;
   vector<vector<float> > *jet_cH_E;
   vector<vector<float> > *jet_cH_charge;
   vector<vector<float> > *jet_cH_Lxy;
   vector<vector<float> > *jet_cH_x;
   vector<vector<float> > *jet_cH_y;
   vector<vector<float> > *jet_cH_z;
   vector<vector<float> > *jet_cH_dRjet;
   
   // List of branches
   TBranch        *b_runnb;   //!
   TBranch        *b_eventnb;   //!
   TBranch        *b_mcchan;   //!
   TBranch        *b_mcwg;   //!
   TBranch        *b_avgmu;   //!
   TBranch        *b_actmu;   //!
   TBranch        *b_PVx;   //!
   TBranch        *b_PVy;   //!
   TBranch        *b_PVz;   //!
   TBranch        *b_truth_PVx;   //!
   TBranch        *b_truth_PVy;   //!
   TBranch        *b_truth_PVz;   //!
   TBranch        *b_njets;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_E;   //!
   TBranch        *b_jet_pt_orig;   //!
   TBranch        *b_jet_eta_orig;   //!
   TBranch        *b_jet_phi_orig;   //!
   TBranch        *b_jet_E_orig;   //!
   TBranch        *b_jet_LabDr_HadF;   //!
   TBranch        *b_jet_DoubleHadLabel;   //!
   TBranch        *b_jet_JVT;   //!
   TBranch        *b_jet_m;   //!
   TBranch        *b_jet_nConst;   //!
   TBranch        *b_jet_dRiso;   //!
   TBranch        *b_jet_truthMatch;   //!
   TBranch        *b_jet_isPU;   //!
   TBranch        *b_jet_aliveAfterOR;   //!
   TBranch        *b_jet_aliveAfterORmu;   //!
   TBranch        *b_jet_isBadMedium;   //!
   TBranch        *b_jet_truthPt;   //!
   TBranch        *b_jet_dRminToB;   //!
   TBranch        *b_jet_dRminToC;   //!
   TBranch        *b_jet_dRminToT;   //!
   TBranch        *b_jet_truthMatch_sm75GeV;   //!
   TBranch        *b_jet_dl1_pb;   //!
   TBranch        *b_jet_dl1_pc;   //!
   TBranch        *b_jet_dl1_pu;   //!
   TBranch        *b_jet_dl1rmu_pb;   //!
   TBranch        *b_jet_dl1rmu_pc;   //!
   TBranch        *b_jet_dl1rmu_pu;   //!
   TBranch        *b_jet_dl1r_pb;   //!
   TBranch        *b_jet_dl1r_pc;   //!
   TBranch        *b_jet_dl1r_pu;   //!
   TBranch        *b_jet_mv2c10;   //!
   TBranch        *b_jet_mv2c10mu;   //!
   TBranch        *b_jet_mv2c10rnn;   //!
   TBranch        *b_jet_mv2c100;   //!
   TBranch        *b_jet_mv2cl100;   //!
   TBranch        *b_jet_mv2c10Flip;   //!
   TBranch        *b_jet_MV2rnnFlip;   //!
   TBranch        *b_jet_MV2r;   //!
   TBranch        *b_jet_DL1Flip_pb;   //!
   TBranch        *b_jet_DL1Flip_pc;   //!
   TBranch        *b_jet_DL1Flip_pu;   //!
   TBranch        *b_jet_mv2c100Flip;   //!
   TBranch        *b_jet_mv2cl100Flip;   //!
   TBranch        *b_jet_dl1rFlip_pb;   //!
   TBranch        *b_jet_dl1rFlip_pc;   //!
   TBranch        *b_jet_dl1rFlip_pu;   //!
   TBranch        *b_jet_MV2rFlip;   //!
   TBranch        *b_jet_ip2d_pb;   //!
   TBranch        *b_jet_ip2d_pc;   //!
   TBranch        *b_jet_ip2d_pu;   //!
   TBranch        *b_jet_ip2d_llr;   //!
   TBranch        *b_jet_ip3d_pb;   //!
   TBranch        *b_jet_ip3d_pc;   //!
   TBranch        *b_jet_ip3d_pu;   //!
   TBranch        *b_jet_ip3d_llr;   //!
   TBranch        *b_jet_ip2;   //!
   TBranch        *b_jet_ip2_c;   //!
   TBranch        *b_jet_ip2_cu;   //!
   TBranch        *b_jet_ip2_nan;   //!
   TBranch        *b_jet_ip2_c_nan;   //!
   TBranch        *b_jet_ip2_cu_nan;   //!
   TBranch        *b_jet_ip3;   //!
   TBranch        *b_jet_ip3_c;   //!
   TBranch        *b_jet_ip3_cu;   //!
   TBranch        *b_jet_ip3_nan;   //!
   TBranch        *b_jet_ip3_c_nan;   //!
   TBranch        *b_jet_ip3_cu_nan;   //!
   TBranch        *b_jet_rnnip_pb;   //!
   TBranch        *b_jet_rnnip_pc;   //!
   TBranch        *b_jet_rnnip_pu;   //!
   TBranch        *b_jet_rnnip_ptau;   //!
   TBranch        *b_jet_ip3dneg_pb;   //!
   TBranch        *b_jet_ip3dneg_pc;   //!
   TBranch        *b_jet_ip3dneg_pu;   //!
   TBranch        *b_jet_ip3dneg_llr;   //!
   TBranch        *b_jet_ip3neg;   //!
   TBranch        *b_jet_ip3neg_c;   //!
   TBranch        *b_jet_ip3neg_cu;   //!
   TBranch        *b_jet_ip3neg_nan;   //!
   TBranch        *b_jet_ip3neg_c_nan;   //!
   TBranch        *b_jet_ip3neg_cu_nan;   //!
   TBranch        *b_jet_ip3dflip_pb;   //!
   TBranch        *b_jet_ip3dflip_pc;   //!
   TBranch        *b_jet_ip3dflip_pu;   //!
   TBranch        *b_jet_ip3dflip_llr;   //!
   TBranch        *b_jet_ip3flip;   //!
   TBranch        *b_jet_ip3flip_c;   //!
   TBranch        *b_jet_ip3flip_cu;   //!
   TBranch        *b_jet_ip3flip_nan;   //!
   TBranch        *b_jet_ip3flip_c_nan;   //!
   TBranch        *b_jet_ip3flip_cu_nan;   //!
   TBranch        *b_jet_sv1_Nvtx;   //!
   TBranch        *b_jet_sv1_ntrkv;   //!
   TBranch        *b_jet_sv1_n2t;   //!
   TBranch        *b_jet_sv1_m;   //!
   TBranch        *b_jet_sv1_efc;   //!
   TBranch        *b_jet_sv1_sig3d;   //!
   TBranch        *b_sv1_llr;   //!
   TBranch        *b_jet_sv1_normdist;   //!
   TBranch        *b_jet_sv1_deltaR;   //!
   TBranch        *b_jet_sv1_Lxy;   //!
   TBranch        *b_jet_sv1_L3d;   //!
   TBranch        *b_jet_sv1_vtx_x;   //!
   TBranch        *b_jet_sv1_vtx_y;   //!
   TBranch        *b_jet_sv1_vtx_z;   //!
   TBranch        *b_jet_nBHadr;   //!
   TBranch        *b_jet_nCHadr;   //!
   TBranch        *b_jet_bH_pdgId;   //!
   TBranch        *b_jet_bH_parent_pdgId;   //!
   TBranch        *b_jet_bH_pt;   //!
   TBranch        *b_jet_bH_eta;   //!
   TBranch        *b_jet_bH_phi;   //!
   TBranch        *b_jet_bH_E;   //!
   TBranch        *b_jet_bH_charge;   //!
   TBranch        *b_jet_bH_Lxy;   //!
   TBranch        *b_jet_bH_x;   //!
   TBranch        *b_jet_bH_y;   //!
   TBranch        *b_jet_bH_z;   //!
   TBranch        *b_jet_bH_dRjet;   //!
   TBranch        *b_jet_cH_pdgId;   //!
   TBranch        *b_jet_cH_parent_pdgId;   //!
   TBranch        *b_jet_cH_pt;   //!
   TBranch        *b_jet_cH_eta;   //!
   TBranch        *b_jet_cH_phi;   //!
   TBranch        *b_jet_cH_E;   //!
   TBranch        *b_jet_cH_charge;   //!
   TBranch        *b_jet_cH_Lxy;   //!
   TBranch        *b_jet_cH_x;   //!
   TBranch        *b_jet_cH_y;   //!
   TBranch        *b_jet_cH_z;   //!
   TBranch        *b_jet_cH_dRjet;   //!
   
   bTag_AntiKt4EMPFlowJets(TTree *tree=0);
   virtual ~bTag_AntiKt4EMPFlowJets();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef bTag_AntiKt4EMPFlowJets_cxx
bTag_AntiKt4EMPFlowJets::bTag_AntiKt4EMPFlowJets(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("flav_Akt4EMPf.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("flav_Akt4EMPf.root");
      }
      f->GetObject("bTag_AntiKt4EMPFlowJets",tree);

  }
   Init(tree);
}

bTag_AntiKt4EMPFlowJets::~bTag_AntiKt4EMPFlowJets()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t bTag_AntiKt4EMPFlowJets::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t bTag_AntiKt4EMPFlowJets::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void bTag_AntiKt4EMPFlowJets::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   jet_pt = 0;
   jet_eta = 0;
   jet_phi = 0;
   jet_E = 0;
   jet_pt_orig = 0;
   jet_eta_orig = 0;
   jet_phi_orig = 0;
   jet_E_orig = 0;
   jet_LabDr_HadF = 0;
   jet_DoubleHadLabel = 0;
   jet_JVT = 0;
   jet_m = 0;
   jet_nConst = 0;
   jet_dRiso = 0;
   jet_truthMatch = 0;
   jet_isPU = 0;
   jet_aliveAfterOR = 0;
   jet_aliveAfterORmu = 0;
   jet_isBadMedium = 0;
   jet_truthPt = 0;
   jet_dRminToB = 0;
   jet_dRminToC = 0;
   jet_dRminToT = 0;
   jet_truthMatch_sm75GeV = 0;
   jet_dl1_pb = 0;
   jet_dl1_pc = 0;
   jet_dl1_pu = 0;
   jet_dl1rmu_pb = 0;
   jet_dl1rmu_pc = 0;
   jet_dl1rmu_pu = 0;
   jet_dl1r_pb = 0;
   jet_dl1r_pc = 0;
   jet_dl1r_pu = 0;
   jet_mv2c10 = 0;
   jet_mv2c10mu = 0;
   jet_mv2c10rnn = 0;
   jet_mv2c100 = 0;
   jet_mv2cl100 = 0;
   jet_mv2c10Flip = 0;
   jet_MV2rnnFlip = 0;
   jet_MV2r = 0;
   jet_DL1Flip_pb = 0;
   jet_DL1Flip_pc = 0;
   jet_DL1Flip_pu = 0;
   jet_mv2c100Flip = 0;
   jet_mv2cl100Flip = 0;
   jet_dl1rFlip_pb = 0;
   jet_dl1rFlip_pc = 0;
   jet_dl1rFlip_pu = 0;
   jet_MV2rFlip = 0;
   jet_ip2d_pb = 0;
   jet_ip2d_pc = 0;
   jet_ip2d_pu = 0;
   jet_ip2d_llr = 0;
   jet_ip3d_pb = 0;
   jet_ip3d_pc = 0;
   jet_ip3d_pu = 0;
   jet_ip3d_llr = 0;
   jet_ip2 = 0;
   jet_ip2_c = 0;
   jet_ip2_cu = 0;
   jet_ip2_nan = 0;
   jet_ip2_c_nan = 0;
   jet_ip2_cu_nan = 0;
   jet_ip3 = 0;
   jet_ip3_c = 0;
   jet_ip3_cu = 0;
   jet_ip3_nan = 0;
   jet_ip3_c_nan = 0;
   jet_ip3_cu_nan = 0;
   jet_rnnip_pb = 0;
   jet_rnnip_pc = 0;
   jet_rnnip_pu = 0;
   jet_rnnip_ptau = 0;
   jet_ip3dneg_pb = 0;
   jet_ip3dneg_pc = 0;
   jet_ip3dneg_pu = 0;
   jet_ip3dneg_llr = 0;
   jet_ip3neg = 0;
   jet_ip3neg_c = 0;
   jet_ip3neg_cu = 0;
   jet_ip3neg_nan = 0;
   jet_ip3neg_c_nan = 0;
   jet_ip3neg_cu_nan = 0;
   jet_ip3dflip_pb = 0;
   jet_ip3dflip_pc = 0;
   jet_ip3dflip_pu = 0;
   jet_ip3dflip_llr = 0;
   jet_ip3flip = 0;
   jet_ip3flip_c = 0;
   jet_ip3flip_cu = 0;
   jet_ip3flip_nan = 0;
   jet_ip3flip_c_nan = 0;
   jet_ip3flip_cu_nan = 0;
   jet_sv1_Nvtx = 0;
   jet_sv1_ntrkv = 0;
   jet_sv1_n2t = 0;
   jet_sv1_m = 0;
   jet_sv1_efc = 0;
   jet_sv1_sig3d = 0;
   sv1_llr = 0;
   jet_sv1_normdist = 0;
   jet_sv1_deltaR = 0;
   jet_sv1_Lxy = 0;
   jet_sv1_L3d = 0;
   jet_sv1_vtx_x = 0;
   jet_sv1_vtx_y = 0;
   jet_sv1_vtx_z = 0;
   jet_nBHadr = 0;
   jet_nCHadr = 0;
   jet_bH_pdgId = 0;
   jet_bH_parent_pdgId = 0;
   jet_bH_pt = 0;
   jet_bH_eta = 0;
   jet_bH_phi = 0;
   jet_bH_E = 0;
   jet_bH_charge = 0;
   jet_bH_Lxy = 0;
   jet_bH_x = 0;
   jet_bH_y = 0;
   jet_bH_z = 0;
   jet_bH_dRjet = 0;
   jet_cH_pdgId = 0;
   jet_cH_parent_pdgId = 0;
   jet_cH_pt = 0;
   jet_cH_eta = 0;
   jet_cH_phi = 0;
   jet_cH_E = 0;
   jet_cH_charge = 0;
   jet_cH_Lxy = 0;
   jet_cH_x = 0;
   jet_cH_y = 0;
   jet_cH_z = 0;
   jet_cH_dRjet = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runnb", &runnb, &b_runnb);
   fChain->SetBranchAddress("eventnb", &eventnb, &b_eventnb);
   fChain->SetBranchAddress("mcchan", &mcchan, &b_mcchan);
   fChain->SetBranchAddress("mcwg", &mcwg, &b_mcwg);
   fChain->SetBranchAddress("avgmu", &avgmu, &b_avgmu);
   fChain->SetBranchAddress("actmu", &actmu, &b_actmu);
   fChain->SetBranchAddress("PVx", &PVx, &b_PVx);
   fChain->SetBranchAddress("PVy", &PVy, &b_PVy);
   fChain->SetBranchAddress("PVz", &PVz, &b_PVz);
   fChain->SetBranchAddress("truth_PVx", &truth_PVx, &b_truth_PVx);
   fChain->SetBranchAddress("truth_PVy", &truth_PVy, &b_truth_PVy);
   fChain->SetBranchAddress("truth_PVz", &truth_PVz, &b_truth_PVz);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_E", &jet_E, &b_jet_E);
   fChain->SetBranchAddress("jet_pt_orig", &jet_pt_orig, &b_jet_pt_orig);
   fChain->SetBranchAddress("jet_eta_orig", &jet_eta_orig, &b_jet_eta_orig);
   fChain->SetBranchAddress("jet_phi_orig", &jet_phi_orig, &b_jet_phi_orig);
   fChain->SetBranchAddress("jet_E_orig", &jet_E_orig, &b_jet_E_orig);
   fChain->SetBranchAddress("jet_LabDr_HadF", &jet_LabDr_HadF, &b_jet_LabDr_HadF);
   fChain->SetBranchAddress("jet_DoubleHadLabel", &jet_DoubleHadLabel, &b_jet_DoubleHadLabel);
   fChain->SetBranchAddress("jet_JVT", &jet_JVT, &b_jet_JVT);
   fChain->SetBranchAddress("jet_m", &jet_m, &b_jet_m);
   fChain->SetBranchAddress("jet_nConst", &jet_nConst, &b_jet_nConst);
   fChain->SetBranchAddress("jet_dRiso", &jet_dRiso, &b_jet_dRiso);
   fChain->SetBranchAddress("jet_truthMatch", &jet_truthMatch, &b_jet_truthMatch);
   fChain->SetBranchAddress("jet_isPU", &jet_isPU, &b_jet_isPU);
   fChain->SetBranchAddress("jet_aliveAfterOR", &jet_aliveAfterOR, &b_jet_aliveAfterOR);
   fChain->SetBranchAddress("jet_aliveAfterORmu", &jet_aliveAfterORmu, &b_jet_aliveAfterORmu);
   fChain->SetBranchAddress("jet_isBadMedium", &jet_isBadMedium, &b_jet_isBadMedium);
   fChain->SetBranchAddress("jet_truthPt", &jet_truthPt, &b_jet_truthPt);
   fChain->SetBranchAddress("jet_dRminToB", &jet_dRminToB, &b_jet_dRminToB);
   fChain->SetBranchAddress("jet_dRminToC", &jet_dRminToC, &b_jet_dRminToC);
   fChain->SetBranchAddress("jet_dRminToT", &jet_dRminToT, &b_jet_dRminToT);
   fChain->SetBranchAddress("jet_truthMatch_sm75GeV", &jet_truthMatch_sm75GeV, &b_jet_truthMatch_sm75GeV);
   fChain->SetBranchAddress("jet_dl1_pb", &jet_dl1_pb, &b_jet_dl1_pb);
   fChain->SetBranchAddress("jet_dl1_pc", &jet_dl1_pc, &b_jet_dl1_pc);
   fChain->SetBranchAddress("jet_dl1_pu", &jet_dl1_pu, &b_jet_dl1_pu);
   fChain->SetBranchAddress("jet_dl1rmu_pb", &jet_dl1rmu_pb, &b_jet_dl1rmu_pb);
   fChain->SetBranchAddress("jet_dl1rmu_pc", &jet_dl1rmu_pc, &b_jet_dl1rmu_pc);
   fChain->SetBranchAddress("jet_dl1rmu_pu", &jet_dl1rmu_pu, &b_jet_dl1rmu_pu);
   fChain->SetBranchAddress("jet_dl1r_pb", &jet_dl1r_pb, &b_jet_dl1r_pb);
   fChain->SetBranchAddress("jet_dl1r_pc", &jet_dl1r_pc, &b_jet_dl1r_pc);
   fChain->SetBranchAddress("jet_dl1r_pu", &jet_dl1r_pu, &b_jet_dl1r_pu);
   fChain->SetBranchAddress("jet_mv2c10", &jet_mv2c10, &b_jet_mv2c10);
   fChain->SetBranchAddress("jet_mv2c10mu", &jet_mv2c10mu, &b_jet_mv2c10mu);
   fChain->SetBranchAddress("jet_mv2c10rnn", &jet_mv2c10rnn, &b_jet_mv2c10rnn);
   fChain->SetBranchAddress("jet_mv2c100", &jet_mv2c100, &b_jet_mv2c100);
   fChain->SetBranchAddress("jet_mv2cl100", &jet_mv2cl100, &b_jet_mv2cl100);
   fChain->SetBranchAddress("jet_mv2c10Flip", &jet_mv2c10Flip, &b_jet_mv2c10Flip);
   fChain->SetBranchAddress("jet_MV2rnnFlip", &jet_MV2rnnFlip, &b_jet_MV2rnnFlip);
   fChain->SetBranchAddress("jet_MV2r", &jet_MV2r, &b_jet_MV2r);
   fChain->SetBranchAddress("jet_DL1Flip_pb", &jet_DL1Flip_pb, &b_jet_DL1Flip_pb);
   fChain->SetBranchAddress("jet_DL1Flip_pc", &jet_DL1Flip_pc, &b_jet_DL1Flip_pc);
   fChain->SetBranchAddress("jet_DL1Flip_pu", &jet_DL1Flip_pu, &b_jet_DL1Flip_pu);
   fChain->SetBranchAddress("jet_mv2c100Flip", &jet_mv2c100Flip, &b_jet_mv2c100Flip);
   fChain->SetBranchAddress("jet_mv2cl100Flip", &jet_mv2cl100Flip, &b_jet_mv2cl100Flip);
   fChain->SetBranchAddress("jet_dl1rFlip_pb", &jet_dl1rFlip_pb, &b_jet_dl1rFlip_pb);
   fChain->SetBranchAddress("jet_dl1rFlip_pc", &jet_dl1rFlip_pc, &b_jet_dl1rFlip_pc);
   fChain->SetBranchAddress("jet_dl1rFlip_pu", &jet_dl1rFlip_pu, &b_jet_dl1rFlip_pu);
   fChain->SetBranchAddress("jet_MV2rFlip", &jet_MV2rFlip, &b_jet_MV2rFlip);
   fChain->SetBranchAddress("jet_ip2d_pb", &jet_ip2d_pb, &b_jet_ip2d_pb);
   fChain->SetBranchAddress("jet_ip2d_pc", &jet_ip2d_pc, &b_jet_ip2d_pc);
   fChain->SetBranchAddress("jet_ip2d_pu", &jet_ip2d_pu, &b_jet_ip2d_pu);
   fChain->SetBranchAddress("jet_ip2d_llr", &jet_ip2d_llr, &b_jet_ip2d_llr);
   fChain->SetBranchAddress("jet_ip3d_pb", &jet_ip3d_pb, &b_jet_ip3d_pb);
   fChain->SetBranchAddress("jet_ip3d_pc", &jet_ip3d_pc, &b_jet_ip3d_pc);
   fChain->SetBranchAddress("jet_ip3d_pu", &jet_ip3d_pu, &b_jet_ip3d_pu);
   fChain->SetBranchAddress("jet_ip3d_llr", &jet_ip3d_llr, &b_jet_ip3d_llr);
   fChain->SetBranchAddress("jet_ip2", &jet_ip2, &b_jet_ip2);
   fChain->SetBranchAddress("jet_ip2_c", &jet_ip2_c, &b_jet_ip2_c);
   fChain->SetBranchAddress("jet_ip2_cu", &jet_ip2_cu, &b_jet_ip2_cu);
   fChain->SetBranchAddress("jet_ip2_nan", &jet_ip2_nan, &b_jet_ip2_nan);
   fChain->SetBranchAddress("jet_ip2_c_nan", &jet_ip2_c_nan, &b_jet_ip2_c_nan);
   fChain->SetBranchAddress("jet_ip2_cu_nan", &jet_ip2_cu_nan, &b_jet_ip2_cu_nan);
   fChain->SetBranchAddress("jet_ip3", &jet_ip3, &b_jet_ip3);
   fChain->SetBranchAddress("jet_ip3_c", &jet_ip3_c, &b_jet_ip3_c);
   fChain->SetBranchAddress("jet_ip3_cu", &jet_ip3_cu, &b_jet_ip3_cu);
   fChain->SetBranchAddress("jet_ip3_nan", &jet_ip3_nan, &b_jet_ip3_nan);
   fChain->SetBranchAddress("jet_ip3_c_nan", &jet_ip3_c_nan, &b_jet_ip3_c_nan);
   fChain->SetBranchAddress("jet_ip3_cu_nan", &jet_ip3_cu_nan, &b_jet_ip3_cu_nan);
   fChain->SetBranchAddress("jet_rnnip_pb", &jet_rnnip_pb, &b_jet_rnnip_pb);
   fChain->SetBranchAddress("jet_rnnip_pc", &jet_rnnip_pc, &b_jet_rnnip_pc);
   fChain->SetBranchAddress("jet_rnnip_pu", &jet_rnnip_pu, &b_jet_rnnip_pu);
   fChain->SetBranchAddress("jet_rnnip_ptau", &jet_rnnip_ptau, &b_jet_rnnip_ptau);
   fChain->SetBranchAddress("jet_ip3dneg_pb", &jet_ip3dneg_pb, &b_jet_ip3dneg_pb);
   fChain->SetBranchAddress("jet_ip3dneg_pc", &jet_ip3dneg_pc, &b_jet_ip3dneg_pc);
   fChain->SetBranchAddress("jet_ip3dneg_pu", &jet_ip3dneg_pu, &b_jet_ip3dneg_pu);
   fChain->SetBranchAddress("jet_ip3dneg_llr", &jet_ip3dneg_llr, &b_jet_ip3dneg_llr);
   fChain->SetBranchAddress("jet_ip3neg", &jet_ip3neg, &b_jet_ip3neg);
   fChain->SetBranchAddress("jet_ip3neg_c", &jet_ip3neg_c, &b_jet_ip3neg_c);
   fChain->SetBranchAddress("jet_ip3neg_cu", &jet_ip3neg_cu, &b_jet_ip3neg_cu);
   fChain->SetBranchAddress("jet_ip3neg_nan", &jet_ip3neg_nan, &b_jet_ip3neg_nan);
   fChain->SetBranchAddress("jet_ip3neg_c_nan", &jet_ip3neg_c_nan, &b_jet_ip3neg_c_nan);
   fChain->SetBranchAddress("jet_ip3neg_cu_nan", &jet_ip3neg_cu_nan, &b_jet_ip3neg_cu_nan);
   fChain->SetBranchAddress("jet_ip3dflip_pb", &jet_ip3dflip_pb, &b_jet_ip3dflip_pb);
   fChain->SetBranchAddress("jet_ip3dflip_pc", &jet_ip3dflip_pc, &b_jet_ip3dflip_pc);
   fChain->SetBranchAddress("jet_ip3dflip_pu", &jet_ip3dflip_pu, &b_jet_ip3dflip_pu);
   fChain->SetBranchAddress("jet_ip3dflip_llr", &jet_ip3dflip_llr, &b_jet_ip3dflip_llr);
   fChain->SetBranchAddress("jet_ip3flip", &jet_ip3flip, &b_jet_ip3flip);
   fChain->SetBranchAddress("jet_ip3flip_c", &jet_ip3flip_c, &b_jet_ip3flip_c);
   fChain->SetBranchAddress("jet_ip3flip_cu", &jet_ip3flip_cu, &b_jet_ip3flip_cu);
   fChain->SetBranchAddress("jet_ip3flip_nan", &jet_ip3flip_nan, &b_jet_ip3flip_nan);
   fChain->SetBranchAddress("jet_ip3flip_c_nan", &jet_ip3flip_c_nan, &b_jet_ip3flip_c_nan);
   fChain->SetBranchAddress("jet_ip3flip_cu_nan", &jet_ip3flip_cu_nan, &b_jet_ip3flip_cu_nan);
   fChain->SetBranchAddress("jet_sv1_Nvtx", &jet_sv1_Nvtx, &b_jet_sv1_Nvtx);
   fChain->SetBranchAddress("jet_sv1_ntrkv", &jet_sv1_ntrkv, &b_jet_sv1_ntrkv);
   fChain->SetBranchAddress("jet_sv1_n2t", &jet_sv1_n2t, &b_jet_sv1_n2t);
   fChain->SetBranchAddress("jet_sv1_m", &jet_sv1_m, &b_jet_sv1_m);
   fChain->SetBranchAddress("jet_sv1_efc", &jet_sv1_efc, &b_jet_sv1_efc);
   fChain->SetBranchAddress("jet_sv1_sig3d", &jet_sv1_sig3d, &b_jet_sv1_sig3d);
   fChain->SetBranchAddress("sv1_llr", &sv1_llr, &b_sv1_llr);
   fChain->SetBranchAddress("jet_sv1_normdist", &jet_sv1_normdist, &b_jet_sv1_normdist);
   fChain->SetBranchAddress("jet_sv1_deltaR", &jet_sv1_deltaR, &b_jet_sv1_deltaR);
   fChain->SetBranchAddress("jet_sv1_Lxy", &jet_sv1_Lxy, &b_jet_sv1_Lxy);
   fChain->SetBranchAddress("jet_sv1_L3d", &jet_sv1_L3d, &b_jet_sv1_L3d);
   fChain->SetBranchAddress("jet_sv1_vtx_x", &jet_sv1_vtx_x, &b_jet_sv1_vtx_x);
   fChain->SetBranchAddress("jet_sv1_vtx_y", &jet_sv1_vtx_y, &b_jet_sv1_vtx_y);
   fChain->SetBranchAddress("jet_sv1_vtx_z", &jet_sv1_vtx_z, &b_jet_sv1_vtx_z);
   fChain->SetBranchAddress("jet_nBHadr", &jet_nBHadr, &b_jet_nBHadr);
   fChain->SetBranchAddress("jet_nCHadr", &jet_nCHadr, &b_jet_nCHadr);
   fChain->SetBranchAddress("jet_bH_pdgId", &jet_bH_pdgId, &b_jet_bH_pdgId);
   fChain->SetBranchAddress("jet_bH_parent_pdgId", &jet_bH_parent_pdgId, &b_jet_bH_parent_pdgId);
   fChain->SetBranchAddress("jet_bH_pt", &jet_bH_pt, &b_jet_bH_pt);
   fChain->SetBranchAddress("jet_bH_eta", &jet_bH_eta, &b_jet_bH_eta);
   fChain->SetBranchAddress("jet_bH_phi", &jet_bH_phi, &b_jet_bH_phi);
   fChain->SetBranchAddress("jet_bH_E", &jet_bH_E, &b_jet_bH_E);
   fChain->SetBranchAddress("jet_bH_charge", &jet_bH_charge, &b_jet_bH_charge);
   fChain->SetBranchAddress("jet_bH_Lxy", &jet_bH_Lxy, &b_jet_bH_Lxy);
   fChain->SetBranchAddress("jet_bH_x", &jet_bH_x, &b_jet_bH_x);
   fChain->SetBranchAddress("jet_bH_y", &jet_bH_y, &b_jet_bH_y);
   fChain->SetBranchAddress("jet_bH_z", &jet_bH_z, &b_jet_bH_z);
   fChain->SetBranchAddress("jet_bH_dRjet", &jet_bH_dRjet, &b_jet_bH_dRjet);
   fChain->SetBranchAddress("jet_cH_pdgId", &jet_cH_pdgId, &b_jet_cH_pdgId);
   fChain->SetBranchAddress("jet_cH_parent_pdgId", &jet_cH_parent_pdgId, &b_jet_cH_parent_pdgId);
   fChain->SetBranchAddress("jet_cH_pt", &jet_cH_pt, &b_jet_cH_pt);
   fChain->SetBranchAddress("jet_cH_eta", &jet_cH_eta, &b_jet_cH_eta);
   fChain->SetBranchAddress("jet_cH_phi", &jet_cH_phi, &b_jet_cH_phi);
   fChain->SetBranchAddress("jet_cH_E", &jet_cH_E, &b_jet_cH_E);
   fChain->SetBranchAddress("jet_cH_charge", &jet_cH_charge, &b_jet_cH_charge);
   fChain->SetBranchAddress("jet_cH_Lxy", &jet_cH_Lxy, &b_jet_cH_Lxy);
   fChain->SetBranchAddress("jet_cH_x", &jet_cH_x, &b_jet_cH_x);
   fChain->SetBranchAddress("jet_cH_y", &jet_cH_y, &b_jet_cH_y);
   fChain->SetBranchAddress("jet_cH_z", &jet_cH_z, &b_jet_cH_z);
   fChain->SetBranchAddress("jet_cH_dRjet", &jet_cH_dRjet, &b_jet_cH_dRjet);
   Notify();
}

Bool_t bTag_AntiKt4EMPFlowJets::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void bTag_AntiKt4EMPFlowJets::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t bTag_AntiKt4EMPFlowJets::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef bTag_AntiKt4EMPFlowJets_cxx
