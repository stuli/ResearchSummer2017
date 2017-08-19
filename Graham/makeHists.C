/******************************
Makes histograms from data and MC. Reweights MC based on single-muon data/MC ratios in Pt and Rapidity.
Histograms are saved to a file.
******************************/

void makeHists()
{
	//Constant parameters
	const int nPtBins = 6; //For macros, Root wants array sizes to be constants
	const int mMin = 6;
	const int mMax = 20;
	
	//Maximum number of entries to read for MC
	int maxMCEntries = 10000000;
	
	//Data and MC files to read in
	TString inDataFileName = "RD2013_pa_1st_run_merged.root";
	TString inMCFileName = "upsilonFlatDimuonMass2BodyNtuple.root";
	
	TString outFileName = "histograms.root";
	
	//Some useful cuts for data
	TCut muplIdCut = "Reco_QQ_mupl_nTrkWMea > 5 && Reco_QQ_mupl_nPixWMea > 1 && Reco_QQ_mupl_normChi2_inner < 1.8 && TMath::Abs(Reco_QQ_mupl_dxy) < 3 && TMath::Abs(Reco_QQ_mupl_dz) < 30 && Reco_QQ_mupl_TrkMuArb && Reco_QQ_mupl_TMOneStaTight && Reco_QQ_VtxProb > 0.01";
	TCut mumiIdCut = "Reco_QQ_mumi_nTrkWMea > 5 && Reco_QQ_mumi_nPixWMea > 1 && Reco_QQ_mumi_normChi2_inner < 1.8 && TMath::Abs(Reco_QQ_mumi_dxy) < 3 && TMath::Abs(Reco_QQ_mumi_dz) < 30 && Reco_QQ_mumi_TrkMuArb && Reco_QQ_mumi_TMOneStaTight";
	TCut muIdCut = muplIdCut+mumiIdCut;
	TCut trigCut = "(HLTriggers&1)==1 && (Reco_QQ_trig&1)==1";
	TCut kinCut = "Reco_QQ_mupl_4mom.Pt() > 4 && Reco_QQ_mumi_4mom.Pt() > 4 && TMath::Abs(Reco_QQ_mupl_4mom.Eta()) < 1.93 && TMath::Abs(Reco_QQ_mumi_4mom.Eta()) < 1.93";
	TCut sameSignCut = "Reco_QQ_sign != 0";
	TCut oppSignCut = "Reco_QQ_sign == 0";
	
	//Convenient compilation of cuts
	TCut activeCut = muIdCut+trigCut+kinCut;
	
	//////////Make data hists//////////
	
	TFile* inDataFile = TFile::Open(inDataFileName,"READ");
	TTree* dataTree = (TTree*)inDataFile->Get("myTree");
	
	//Make same-sign mass histogram from data
	TH1D* mH[nPtBins];
	mH[0] = new TH1D("mH_0","p_{T} < 5",100,mMin,mMax); //Index 0 is integrated pt 0-5 bin
	dataTree->Project("mH_0","Reco_QQ_4mom.M()",activeCut + sameSignCut + "Reco_QQ_4mom.Pt() < 5");
	for (int h = 1; h < nPtBins; h++)
	{
		mH[h] = new TH1D(Form("mH_%d",h),Form("%d < p_{T} < %d",h-1,h),100,mMin,mMax);
		TCut binCut = Form("Reco_QQ_4mom.Pt() > %d && Reco_QQ_4mom.Pt() < %d",h-1,h);
		dataTree->Project(Form("mH_%d",h),"Reco_QQ_4mom.M()",activeCut + sameSignCut + binCut);
	}
	
	//Make opposite-sign mass histogram from data
	TH1D* mH_oppS[nPtBins];
	mH_oppS[0] = new TH1D("mH_oppS_0","p_{T} < 5",100,mMin,mMax); //Index 0 is integrated pt 0-5 bin
	dataTree->Project("mH_oppS_0","Reco_QQ_4mom.M()",activeCut + oppSignCut + "Reco_QQ_4mom.Pt() < 5");
	for (int h = 1; h < nPtBins; h++)
	{
		mH_oppS[h] = new TH1D(Form("mH_oppS_%d",h),Form("%d < p_{T} < %d",h-1,h),100,mMin,mMax);
		TCut binCut = Form("Reco_QQ_4mom.Pt() > %d && Reco_QQ_4mom.Pt() < %d",h-1,h);
		dataTree->Project(Form("mH_oppS_%d",h),"Reco_QQ_4mom.M()",activeCut + oppSignCut + binCut);
	}
	
	//PtEta histogram parameters
	double ptRap_rMin = -1.93;
	double ptRap_rMax = 1.93;
	int ptRap_rnBins = 100;
	double ptRap_rBinSize = (ptRap_rMax - ptRap_rMin)/((double)ptRap_rnBins);
	double ptRap_ptMin = 4;
	double ptRap_ptMax = 14;
	int ptRap_ptnBins = 100;
	double ptRap_ptBinSize = (ptRap_ptMax - ptRap_ptMin)/((double)ptRap_ptnBins);
	
	//Make data PtEta hist without kincut
	TH2D* muplH = new TH2D("muplH","Data MuPlus Eta and p_{T}",ptRap_rnBins,ptRap_rMin,ptRap_rMax,ptRap_ptnBins,ptRap_ptMin,ptRap_ptMax);
	dataTree->Project("muplH","Reco_QQ_mupl_4mom.Pt():Reco_QQ_mupl_4mom.Eta()",muIdCut+trigCut+"Reco_QQ_4mom.Pt()<5");
	
	//////////Make MC hists//////////
	
	//Need to loop through MC manually because entries will be reweighted
	TFile* MCFile = TFile::Open(inMCFileName,"READ");
	TNtuple* MCtuple = (TNtuple*)MCFile->Get("upsilonNtuple");
	Float_t UpsM, UpsRap, UpsPt, MuPPt, MuMPt, MuPy, MuMy;
	MCtuple->SetBranchAddress("UpsM",&UpsM);
	MCtuple->SetBranchAddress("UpsRap",&UpsRap);
	MCtuple->SetBranchAddress("UpsPt",&UpsPt);
	MCtuple->SetBranchAddress("MuPPt",&MuPPt);
	MCtuple->SetBranchAddress("MuMPt",&MuMPt);
	MCtuple->SetBranchAddress("MuPy",&MuPy);
	MCtuple->SetBranchAddress("MuMy",&MuMy);
	
	//Weighting functions for MC
	TF1* mfunc = new TF1("mfunc","[0]*exp(-x/[1])",6,20);
	mfunc->SetParameters(100,4.2);
	mfunc->SetParNames("A","#lambda");
	TF1* rfunc = new TF1("upsRapidityFunc","gaus(0)+gaus(3)",-5,5);
	rfunc->SetParameters(70,0.0,1.0,550,0.0,0.23);
	rfunc->SetParNames("A1","mean1","sigma1","A2","mean2","sigma2");
	TF1* ptfunc = new TF1("ptfunc","[0]*x/(exp(x/[1])+1)",0,5);
	ptfunc->SetParameters(100,2.71);
	ptfunc->SetParNames("A","T");
	
	//Make MC PtRap histogram
	TH2D* MCmuplH = new TH2D("MCmuplH","MC MuPlus Rapidity and p_{T}",ptRap_rnBins,ptRap_rMin,ptRap_rMax,ptRap_ptnBins,ptRap_ptMin,ptRap_ptMax);
	
	Int_t nentries1 = MCtuple->GetEntries();
	if (nentries1 > maxMCEntries)
		nentries1 = maxMCEntries;
	for (Int_t i = 0; i < nentries1; i++)
	{
		MCtuple->GetEntry(i);	
		
		//Apply weighting
		double mWeight = mfunc->Eval(UpsM);
		double rWeight = rfunc->Eval(UpsRap);
		double ptWeight = ptfunc->Eval(UpsPt);
		
		double weight = mWeight * rWeight * ptWeight;
		
		//Apply kinematic cuts
		if (MuPPt>4 && MuMPt>4 && abs(MuPy)<1.93 && abs(MuMy)<1.93)
		{
			MCmuplH->Fill(MuPy,MuPPt,weight);
		}
	}
	
	//Take ratio of data to MC
	TH2D* muplRatioH = (TH2D*)muplH->Clone("muplRatioH");
	muplRatioH->SetTitle("Ratio of data/MC");
	muplRatioH->Divide(MCmuplH);
	
	//Make reweighted MC histograms
	TH2D* MCmuplWeightedH = new TH2D("MCmuplWeightedH","Weighted MC MuPlus Rapidity and p_{T}",ptRap_rnBins,ptRap_rMin,ptRap_rMax,ptRap_ptnBins,ptRap_ptMin,ptRap_ptMax);
	TH1D* MCmH[nPtBins]; //Same-sign dimuon mass histograms from MC
	MCmH[0] = new TH1D("MCmH_0","pT < 5",100,mMin,mMax); //Index 0 is integrated pt 0-5 bin
	for (int h = 1; h < nPtBins; h++)
	{
		MCmH[h] = new TH1D(Form("MCmH_%d",h),Form("%d < pT < %d",h-1,h),100,mMin,mMax);
	}
	
	Int_t nentries = MCtuple->GetEntries();
	if (nentries > maxMCEntries)
		nentries = maxMCEntries;
	for (Int_t i = 0; i < nentries; i++)
	{
		MCtuple->GetEntry(i);	
		
		//Apply weighting
		int rbin = TMath::Ceil((MuPy - ptRap_rMin)/ptRap_rBinSize);
		int ptbin = TMath::Ceil((MuPPt - ptRap_ptMin)/ptRap_ptBinSize);
		double muplWeight = muplRatioH->GetBinContent(rbin,ptbin);
		
		double mWeight = mfunc->Eval(UpsM);
		double rWeight = rfunc->Eval(UpsRap);
		double ptWeight = ptfunc->Eval(UpsPt);
		
		double weight = muplWeight * mWeight * rWeight * ptWeight;
		
		//Apply kinematic cuts
		if (MuPPt>4 && MuMPt>4 && abs(MuPy)<1.93 && abs(MuMy)<1.93)
		{
			MCmuplWeightedH->Fill(MuPy,MuPPt,weight);
			
			if (UpsPt < 5)
			{
				MCmH[0]->Fill(UpsM,weight);
				for (int j = 1; j<nPtBins; j++)
				{
					if (UpsPt < j)
					{
						MCmH[j]->Fill(UpsM,weight);
						break;
					}
				}
			}
		}
	}
	
	//Rescale MC histograms to have same integral as data.
	for (int h = 0; h < nPtBins; h++)
		MCmH[h]->Scale((mH[h]->Integral())/(MCmH[h]->Integral()));
	
	//Write histograms to file
	TFile* histOutFile = new TFile(outFileName,"RECREATE");
	for (int h = 0; h < nPtBins; h++)
	{
		mH[h]->Write();
		mH_oppS[h]->Write();
		MCmH[h]->Write();
	}
	muplH->Write();
	MCmuplH->Write();
	muplRatioH->Write();
	MCmuplWeightedH->Write();
}