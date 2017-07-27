#include <iostream>
#include <string>
#include <sstream>
using namespace std;
void reproduceAN2(){
	gStyle->SetOptStat(0);
	TChain myTree("myTree");
  myTree.Add("OniaTree_JpsiMM_5p02TeV_TuneCUETP8M1_ptJpsi912_noCUT.root");
  cout<<"Entries = "<<myTree.GetEntries()<<endl;
    //const Float_t xbins[] = {3,6.5,9,12};
  const Float_t xbins[] = {3,6.5,9,12};
  std::string acceptanceCut_num = "((Reco_QQ_mumi_4mom.Pt()>3.5&&TMath::Abs(Reco_QQ_mumi_4mom.Eta())<1.2)||(Reco_QQ_mumi_4mom.Pt()>(5.77-1.8*TMath::Abs(Reco_QQ_mumi_4mom.Eta()))&&TMath::Abs(Reco_QQ_mumi_4mom.Eta())>=1.2&&TMath::Abs(Reco_QQ_mumi_4mom.Eta())<2.1)||(Reco_QQ_mumi_4mom.Pt()>1.8&&TMath::Abs(Reco_QQ_mumi_4mom.Eta())>=2.1&&TMath::Abs(Reco_QQ_mumi_4mom.Eta())<2.4))&&((Reco_QQ_mupl_4mom.Pt()>3.5&&TMath::Abs(Reco_QQ_mupl_4mom.Eta())<1.2)||(Reco_QQ_mupl_4mom.Pt()>(5.77-1.8*TMath::Abs(Reco_QQ_mupl_4mom.Eta()))&&TMath::Abs(Reco_QQ_mupl_4mom.Eta())>=1.2&&TMath::Abs(Reco_QQ_mupl_4mom.Eta())<2.1)||(Reco_QQ_mupl_4mom.Pt()>1.8&&TMath::Abs(Reco_QQ_mupl_4mom.Eta())>=2.1&&TMath::Abs(Reco_QQ_mupl_4mom.Eta())<2.4))"; 
  string acceptanceCut_denom = "((Gen_QQ_mumi_4mom.Pt()>3.5&&TMath::Abs(Gen_QQ_mumi_4mom.Eta())<1.2)||(Gen_QQ_mumi_4mom.Pt()>(5.77-1.8*TMath::Abs(Gen_QQ_mumi_4mom.Eta()))&&TMath::Abs(Gen_QQ_mumi_4mom.Eta())>=1.2&&TMath::Abs(Gen_QQ_mumi_4mom.Eta())<2.1)||(Gen_QQ_mumi_4mom.Pt()>1.8&&TMath::Abs(Gen_QQ_mumi_4mom.Eta())>=2.1&&TMath::Abs(Gen_QQ_mumi_4mom.Eta())<2.4))&&((Gen_QQ_mupl_4mom.Pt()>3.5&&TMath::Abs(Gen_QQ_mupl_4mom.Eta())<1.2)||(Gen_QQ_mupl_4mom.Pt()>(5.77-1.8*TMath::Abs(Gen_QQ_mupl_4mom.Eta()))&&TMath::Abs(Gen_QQ_mupl_4mom.Eta())>=1.2&&TMath::Abs(Gen_QQ_mupl_4mom.Eta())<2.1)||(Gen_QQ_mupl_4mom.Pt()>1.8&&TMath::Abs(Gen_QQ_mupl_4mom.Eta())>=2.1&&TMath::Abs(Gen_QQ_mupl_4mom.Eta())<2.4))"; 
  string analysisCut = "(Reco_QQ_trig&1)==1&&(HLTriggers&1)==1&&Reco_QQ_mupl_TMOneStaTight>=1&&Reco_QQ_mupl_nPixWMea>0&&Reco_QQ_mupl_nTrkWMea>5&&Reco_QQ_mupl_dxy<.3&&Reco_QQ_mupl_dz<20&&Reco_QQ_mupl_dz>-20&&Reco_QQ_mupl_isGoodMuon&&Reco_QQ_mumi_TMOneStaTight>=1&&Reco_QQ_mumi_nPixWMea>0&&Reco_QQ_mumi_nTrkWMea>5&&Reco_QQ_mumi_dxy<.3&&Reco_QQ_mumi_dz<20&&Reco_QQ_mumi_dz>-20&&Reco_QQ_mumi_isGoodMuon&&Reco_QQ_VtxProb>.01";
  string ppDecayCut = "Reco_QQ_ctau3D<.012+.23/Reco_QQ_4mom.Pt()";
  string massChargeCut = "Reco_QQ_4mom.M()>2.946&&Reco_QQ_4mom.M()<3.346";
  string matchCut = "TMath::Sqrt(TMath::Power(Reco_QQ_mupl_4mom.Eta()-Gen_QQ_mupl_4mom.Eta(),2)+TMath::Power(Reco_QQ_mupl_4mom.Phi()-Gen_QQ_mupl_4mom.Phi(),2))<.03&&TMath::Sqrt(TMath::Power(Reco_QQ_mumi_4mom.Eta()-Gen_QQ_mumi_4mom.Eta(),2)+TMath::Power(Reco_QQ_mumi_4mom.Phi()-Gen_QQ_mumi_4mom.Phi(),2))<.03";
  string otherCut_num = "TMath::Abs(Reco_QQ_4mom.Rapidity())<1.6&&TMath::Abs(Reco_QQ_4mom.Rapidity())<2.4";
  string otherCut_denom = "TMath::Abs(Gen_QQ_4mom.Rapidity())<1.6&&TMath::Abs(Gen_QQ_4mom.Rapidity())<2.4";

  TH1D *histoPt_denom = new TH1D("histoPt_denom","Dimuon P_{T}",100,0,12);
  histoPt_denom->GetXaxis()->SetTitle("P_{T} (GeV)");
  histoPt_denom->SetLineColor(kBlack);
  histoPt_denom->SetFillColor(kBlack);
  histoPt_denom->SetFillStyle(4050);
  histoPt_denom->SetLineWidth(3);

  TH1D *histoPt_num = new TH1D("histoPt_num","Dimuon P_{T} (After Cuts)",100,0,12);
  histoPt_num->GetXaxis()->SetTitle("P_{T} (GeV)");
  histoPt_num->GetXaxis()->SetTitle("P_{T} (GeV)");
  histoPt_num->SetLineColor(kRed);
  histoPt_num->SetFillColorAlpha(kRed,.35);
  histoPt_num->SetFillStyle(4050);
  histoPt_num->SetLineWidth(3);

  TH1D *histoMass_num = new TH1D("histoMass_num","Dimuon Mass",100,2.5,3.8);
  histoMass_num->GetXaxis()->SetTitle("M (GeV)");

  TH1D *histoPt_eff = new TH1D("histoPt_eff","Dimuon P_{T} Eff.",100,0,12);
  histoPt_eff->GetXaxis()->SetTitle("P_{T} (GeV)");
  histoPt_eff->GetYaxis()->SetRange(0,1);

  const char *histo_DL_cut = (analysisCut+"&&"+massChargeCut+"&&"+otherCut_num).c_str();
  const char *denom_cut = (acceptanceCut_denom+"&&"+otherCut_denom).c_str();
  const char *num_cut = (acceptanceCut_num+"&&"+ppDecayCut+"&&"+massChargeCut+"&&"+matchCut+"&&"+otherCut_num+"&&"+analysisCut).c_str();

  TCanvas *c0 = new TCanvas();
  c0->SetTickx();
  c0->SetTicky();
  myTree.Draw("Gen_QQ_4mom.Pt()>>histoPt_denom",denom_cut,"COLZ");
  Double_t scale = 1/histoPt_denom->Integral();
  histoPt_denom->Scale(scale*histoPt_denom->GetEntries());

   	//TCanvas *c1 = new TCanvas();
  myTree.Draw("Reco_QQ_4mom.Pt()>>histoPt_num",num_cut,"COLZSAME");
  scale = 1/histoPt_num->Integral();
  histoPt_num->Scale(scale*histoPt_num->GetEntries());

  cout<<"Entries in num: "<<histoPt_num->GetEntries()<<endl;
  cout<<"Entries in denom: "<<histoPt_denom->GetEntries()<<endl;

  TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->SetHeader("Gen and Reco P_{T}"); // option "C" allows to center the header
  legend->AddEntry(histoPt_num,"Reco P_{T} (numerator) ","f");
  legend->AddEntry(histoPt_denom,"Gen P_{T} (denominator) ","f");
  legend->Draw("SAME");

  TCanvas *c2 = new TCanvas();
  myTree.Draw("Reco_QQ_4mom.M()>>histoMass_num",num_cut,"COLZ");
	//TF1 *f = new TF1("f","gaus(0)+[3]+[4]*x",2.5,3.8);
	//f->SetParameters();

  histoPt_eff->Sumw2();
  histoPt_eff->Divide(histoPt_num,histoPt_denom,1,1,"B");

  TCanvas *c3 = new TCanvas();
  histoPt_eff->SetMaximum(1.5);
  histoPt_eff->SetMinimum(0);
  histoPt_eff->Draw("COLZ");

  TCanvas *c4 = new TCanvas();
  TH1D* histo_DL = new TH1D("histo_DL","Decay Length",100,-.15,.2);
  myTree.Draw("Reco_QQ_ctau3D>>histo_DL",histo_DL_cut);
  Double_t scale0 = 1/histo_DL->Integral();
  histo_DL->Scale(scale0*histo_DL->GetEntries());
  histo_DL->Draw();

  cout<<"Efficiency: "<<histoPt_num->GetEntries()/histoPt_denom->GetEntries()<<endl;

  //cout<<"Bin #    	Eff.    	Err. 		0th bin is underflow"<<endl;
 	/*for(int i = 0; i<=histoPt_eff->GetNbinsX();i++){

 		cout<<i<<"\t"<<histoPt_eff->GetBinContent(i)<<"\t"<<histoPt_eff->GetBinError(i)<<endl;
 	}*/
}