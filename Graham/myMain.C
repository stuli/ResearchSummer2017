/*
Generate and fit mass/rapidity/pt plots, compare between data and MC. Uses reduced "mrp" tuples for data and MC tuple.
*/

#include "myFitFunctions.C"

void myMain()
{
	TFile* inFile = TFile::Open("mrpTuple_sameS_trigcut.root","READ");
	TNtuple* mrpTuple = (TNtuple*)inFile->Get("mrpTuple");
	
	
	///////////Make hists from trig-cut data
	double minM = 6;
	double maxM = 20;
	const int nbinsM = 14;
	double binsizeM = (maxM - minM)/nbinsM;
	TH1D* rH[nbinsM];
	TH1D* ptH[nbinsM];
	for (int h = 0; h < nbinsM; h++)
	{
		double binminM = minM + h * binsizeM;
		double binmaxM = minM + (h+1)*binsizeM;
		
		rH[h] = new TH1D(Form("rH_%d",h),Form("%d < M < %d",h+6,h+7),100,-2.4,2.4);
		mrpTuple->Project(Form("rH_%d",h),"Rapidity",Form("Pt < 5 && M > %d && M < %d",6+h,7+h));
		
		ptH[h] = new TH1D(Form("ptH_%d",h),Form("%d < M < %d",h+6,h+7),100,0,5);
		mrpTuple->Project(Form("ptH_%d",h),"Pt",Form("Pt < 5 && M > %d && M < %d",6+h,7+h));
	}
	
	////////////Fit rapidity hists
	TF1* rF1[nbinsM];
	for (int h = 0; h < nbinsM; h++)
	{
		rF1[h] = new TF1(Form("rFit1_%d", h),fRapidity,-2.4,2.4,7);
		rF1[h]->SetParNames("A","pgmu","pgsig","ratio","sgmu","sgsig","sgerf");
		if (h <= 5)
		{
			rF1[h]->SetParameters(1000,0.2,0.2,0.8,1,1,2);
			rF1[h]->SetParLimits(0,500,10000);
			rF1[h]->SetParLimits(1,0,0.35);
			rF1[h]->SetParLimits(2,0.025,0.25);
			rF1[h]->SetParLimits(3,0.2,3);
			rF1[h]->SetParLimits(4,0,1.5);
			rF1[h]->SetParLimits(5,0.25,1.2);
			//rF1[h]->SetParLimits(6,0,4);
			rF1[h]->SetParLimits(6,0,0);
			rF1[h]->FixParameter(6,0);
		}
		else
		{
			rF1[h]->SetParameters(1000,0,0.2,0.2,0,1,0);
			rF1[h]->SetParLimits(0,500,10000);
			rF1[h]->SetParLimits(1,0,0);
			rF1[h]->SetParLimits(2,0.025,0.25);
			rF1[h]->SetParLimits(3,0.2,0.5);
			rF1[h]->SetParLimits(4,0,0);
			rF1[h]->SetParLimits(5,0.25,1.2);
			rF1[h]->SetParLimits(6,0,0);
			rF1[h]->FixParameter(1,0);
			rF1[h]->FixParameter(6,0);
		}
		rH[h]->Fit(rF1[h]);
	}
	//Fit pT hists
	TF1* ptF1[nbinsM];
	for (int h = 0; h < nbinsM; h++)
	{
		ptF1[h] = new TF1(Form("ptFit1_%d",h),"[0]*x/(exp(x/[1])+1)+gaus(2)",0,40);
		if (h == 12)
			ptF1[h]->SetParameters(700,2.2,200,4,5); // GeV/c^2
		else
			ptF1[h]->SetParameters(100,2.71,200,3,1); // GeV/c^2
		ptF1[h]->SetParLimits(2,100,1000);
		ptF1[h]->SetParNames("A","T","gausA","gausMu","gausSigma");
		ptH[h]->Fit(ptF1[h]);
	}
	double ptParam[5][nbinsM];
	for (int h = 0; h < nbinsM; h++)
	{
		TF1* hf = ptH[h]->GetFunction(Form("ptFit1_%d", h));
		for (int p = 0; p < 5; p++)
		{
			ptParam[p][h] = hf->GetParameter(p);
		}
	}
	
	////////////Draw data hists
	int hFirst = 0;
	int hLast = nbinsM-1;
	TCanvas* c1 = new TCanvas("c1","Rapidity Plots",800,600);
	int nPlots = hLast-hFirst+1;
	int nColumns = TMath::Floor(TMath::Sqrt(nPlots));
	int nRows = TMath::Ceil(((double)nPlots)/nColumns);
	c1->Divide(nColumns, nRows);
	//Draw rapidity hists
	for (int h = hFirst, i = 1; h <= hLast; h++, i++)
	{
		c1->cd(i);
		
		rH[h]->SetMinimum(0);
		rH[h]->Draw();
		rH[h]->SetXTitle("Rapidity"); 
		rH[h]->SetYTitle("count");
		rH[h]->GetXaxis()->SetTitleSize(0.05);
		rH[h]->GetXaxis()->SetLabelSize(0.06);
		rH[h]->GetYaxis()->SetTitleSize(0.05);
		rH[h]->GetYaxis()->SetLabelSize(0.06);
		rH[h]->SetStats(0);
		//gStyle->SetOptStat(1000000001);
		//gStyle->SetOptFit(1);
		gStyle->SetTitleFontSize(.1);
	}
	//Draw pt hists
	TCanvas* cpt = new TCanvas("cpt","pT Plots",800,600);
	cpt->Divide(nColumns, nRows);
	for (int h = hFirst, i = 1; h <= hLast; h++, i++)
	{
		cpt->cd(i);
		
		ptH[h]->SetMinimum(0);
		ptH[h]->Draw();
		ptH[h]->SetXTitle("pT"); 
		ptH[h]->SetYTitle("count");
		ptH[h]->GetXaxis()->SetTitleSize(0.05);
		ptH[h]->GetXaxis()->SetLabelSize(0.06);
		ptH[h]->GetYaxis()->SetTitleSize(0.05);
		ptH[h]->GetYaxis()->SetLabelSize(0.06);
		ptH[h]->SetStats(0);
		//gStyle->SetOptStat(1000000001);
		//gStyle->SetOptFit(1);
	}
	
	////////////Plot parameters of rapidity distribution
	TF1* hf0 = rH[0]->GetFunction(Form("rFit1_%d", 0));
	int nPars = 7;
	double param[7][nbinsM];
	TGraph* paramGraph[7];
	for (int h = 0; h < nbinsM; h++)
	{
		TF1* hf = rH[h]->GetFunction(Form("rFit1_%d", h));
		for (int p = 0; p < nPars; p++)
		{
			param[p][h] = hf->GetParameter(p);
		}
	}
	double graphxvals[nbinsM];
	for (int h = 0; h < nbinsM; h++)
		graphxvals[h] = h+minM;
	for (int p = 0; p < nPars; p++)
	{
		paramGraph[p] = new TGraph(nbinsM,graphxvals,param[p]);
		//paramGraph[p]->SetTitle(hf0->GetParName());
	}
	paramGraph[0]->SetTitle("A");
	paramGraph[1]->SetTitle("pgmu");
	paramGraph[2]->SetTitle("pgsig");
	paramGraph[3]->SetTitle("ratio");
	paramGraph[4]->SetTitle("sgmu");
	paramGraph[5]->SetTitle("sgsig");
	paramGraph[6]->SetTitle("sgerf");
	
	TCanvas* cpar = new TCanvas("cpar","Parameter Plots",800,600);
	nPlots = nPars;
	nColumns = TMath::Floor(TMath::Sqrt(nPlots));
	nRows = TMath::Ceil(((double)nPlots)/nColumns);
	cpar->Divide(nColumns, nRows);
	for (int p = 0; p < nPars; p++)
	{
		cpar->cd(p+1);
		
		paramGraph[p]->SetMinimum(0);
		paramGraph[p]->SetMarkerStyle(21);
		paramGraph[p]->Draw("AP");
	}
	
	cout << "\nPARAMETERS\n";
	for (int p = 0; p < nPars; p++)
	{
		cout << "{" << param[p][0];
		for (int h = 1; h < nbinsM; h++)
			cout << "," << param[p][h];
		cout << "}\n";
	}
	
	////////////Make hists from all-cut data
	TFile * inFileAllcut = TFile::Open("mrpTuple_sameS_allcut.root","READ");
	TNtuple* mrpTupleAllcut = (TNtuple*)inFileAllcut->Get("mrpTuple");
	
	const int nPtBins = 6;
	TH1D* mH[6];
	mH[0] = new TH1D("mH_0","pT < 5",100,minM,maxM);
	mrpTupleAllcut->Project("mH_0","M","Pt<5");
	for (int h = 1; h < nPtBins; h++)
	{
		mH[h] = new TH1D(Form("mH_%d",h),Form("%d < pT < %d",h-1,h),100,minM,maxM);
		mrpTupleAllcut->Project(Form("mH_%d",h),"M",Form("Pt>%d && Pt<%d",h-1,h));
	}
	TH1D* rAllcutH[nbinsM];
	for (int h = 0; h < nbinsM; h++)
	{
		rAllcutH[h] = new TH1D(Form("rAllcutH_%d",h),Form("%d < M < %d",h+6,h+7),100,-2.4,2.4);
		mrpTupleAllcut->Project(Form("rAllcutH_%d",h),"Rapidity",Form("Pt < 5 && M > %d && M < %d",6+h,7+h));
	}
	
	//Draw mass hists
	TCanvas* cM = new TCanvas("cM","Mass Plots",800,600);
	cM->Divide(2, 3);
	
	for (int h = 0, i = 1; h < nPtBins; h++, i++)
	{
		cM->cd(i);
		
		mH[h]->SetMinimum(0);
		mH[h]->Draw();
		mH[h]->GetXaxis()->SetTitleSize(0.05);
		mH[h]->GetXaxis()->SetLabelSize(0.06);
		mH[h]->GetYaxis()->SetTitleSize(0.05);
		mH[h]->GetYaxis()->SetLabelSize(0.06);
		mH[h]->SetXTitle("Mass (GeV)"); 
		mH[h]->SetYTitle("count");
		mH[h]->SetStats(0);
		//gStyle->SetOptStat(1000000001);
		//gStyle->SetOptFit(1);
	}
	//Draw allcut rapidity hists
	TCanvas* crac = new TCanvas("crac","Rapidity allcut Plots",800,600);
	crac->Divide(3, 5);
	for (int h = hFirst, i = 1; h <= hLast; h++, i++)
	{
		crac->cd(i);
		
		rAllcutH[h]->SetMinimum(0);
		rAllcutH[h]->Draw();
		rAllcutH[h]->GetXaxis()->SetTitleSize(0.05);
		rAllcutH[h]->GetXaxis()->SetLabelSize(0.06);
		rAllcutH[h]->GetYaxis()->SetTitleSize(0.05);
		rAllcutH[h]->GetYaxis()->SetLabelSize(0.06);
		rAllcutH[h]->SetXTitle("Rapidity"); 
		rAllcutH[h]->SetYTitle("count");
		rAllcutH[h]->SetStats(0);
		//gStyle->SetOptStat(1000000001);
		//gStyle->SetOptFit(1);
		gStyle->SetTitleFontSize(.1);
	}
	
	////////////Make MC hists
	TFile* MCFile = TFile::Open("upsilonFlatDimuonMass2BodyNtuple.root","READ");
	TNtuple* MCtuple = (TNtuple*)MCFile->Get("upsilonNtuple");
	Float_t UpsM, UpsRap, UpsPt, MuPPt, MuMPt, MuPy, MuMy;
	MCtuple->SetBranchAddress("UpsM",&UpsM);
	MCtuple->SetBranchAddress("UpsRap",&UpsRap);
	MCtuple->SetBranchAddress("UpsPt",&UpsPt);
	MCtuple->SetBranchAddress("MuPPt",&MuPPt);
	MCtuple->SetBranchAddress("MuMPt",&MuMPt);
	MCtuple->SetBranchAddress("MuPy",&MuPy);
	MCtuple->SetBranchAddress("MuMy",&MuMy);
	
	//Rapidity hists
	TH1D* MCrH[nbinsM];
	TH1D* MCrAllcutH[nbinsM];
	for (int h = 0; h < nbinsM; h++)
	{
		double binminM = minM + h * binsizeM;
		double binmaxM = minM + (h+1)*binsizeM;
		
		MCrH[h] = new TH1D(Form("MCrH_%d",h),Form("%d < M < %d",h+6,h+7),100,-2.4,2.4);
		MCrAllcutH[h] = new TH1D(Form("MCrAllcutH_%d",h),Form("%d < M < %d",h+6,h+7),100,-2.4,2.4);
	}
	//PT hists
	TH1D* MCptH[nbinsM];
	for (int h = 0; h < nbinsM; h++)
	{
		MCptH[h] = new TH1D(Form("MCptH_%d",h),Form("%d < M < %d",h+6,h+7),100,0,5);
	}
	//Mass hists
	TH1D* MCmH[nPtBins];
	MCmH[0] = new TH1D("mH_0","pT < 5",100,minM,maxM);
	for (int h = 1; h < nPtBins; h++)
	{
		MCmH[h] = new TH1D(Form("mH_%d",h),Form("%d < pT < %d",h-1,h),100,minM,maxM);
	}
	
	//Loop through MC entries
	Int_t nentries = MCtuple->GetEntries();
	for (Int_t i = 0; i < nentries; i++)
	{
		MCtuple->GetEntry(i);
		for (int h = 0; h < nbinsM; h++)
		{
			//bool muplcut = (MuPPt>3.5 && abs(MuPy)<1.2) || (MuPPt>(5.77-1.8*abs(MuPy)) && abs(MuPy)>=1.2 && abs(MuPy)<2.1) || (MuPPt>1.8 && abs(MuPy)>=2.1 && abs(MuPy)<2.4);
			//bool mumicut = (MuMPt>3.5 && abs(MuMy)<1.2) || (MuMPt>(5.77-1.8*abs(MuMy)) && abs(MuMy)>=1.2 && abs(MuMy)<2.1) || (MuMPt>1.8 && abs(MuMy)>=2.1 && abs(MuMy)<2.4);
			//bool mucut = muplcut && mumicut;
			
			//Fill Rapidity hists
			if (UpsM > h+6 && UpsM < h+7 && UpsPt < 5)
			{	
				//Apply weighting to mass and rapidity distributions
				rfunc = new TF1("rfunc",fRapidity,-5,5,7);
				rfunc->SetParNames("A","pgmu","pgsig","ratio","sgmu","sgsig","sgerf");
				rfunc->SetParameters(param[0][h],param[1][h],param[2][h],param[3][h],param[4][h],param[5][h],param[6][h]);
				/*TF1* rfunc = new TF1("upsRapidityFunc","gaus(0)+gaus(3)",-5,5);
				rfunc->SetParameters(70,0.0,1.0,550,0.0,0.23);
				rfunc->SetParNames("A1","mean1","sigma1","A2","mean2","sigma2");*/
				mfunc = new TF1("mfunc","[0]*exp(-x/[1])",6,20);
				mfunc->SetParameters(100,3.4);//3.4 before
				mfunc->SetParNames("A","#lambda");
				//Double_t weight = rfunc->Eval(UpsRap)*mfunc->Eval(UpsM);
				ptfunc = new TF1("ptfunc","[0]*x/(exp(x/[1])+1)+gaus(2)",0,5);
				ptfunc->SetParNames("A","T","gausA","gausMu","gausSigma");
				ptfunc->SetParameters(ptParam[0][h],ptParam[1][h],ptParam[2][h],ptParam[3][h],ptParam[4][h]);
				Double_t weight = rfunc->Eval(UpsRap)*mfunc->Eval(UpsM)*ptfunc->Eval(UpsPt);
				
				MCrH[h]->Fill(UpsRap,weight);
				MCptH[h]->Fill(UpsPt,weight);
				
				//Apply all cuts to mass hists
				if (MuPPt>4 && MuMPt>4 && abs(MuPy)<1.93 && abs(MuMy)<1.93 && UpsPt<5.0 && UpsPt>0.0)
				{
					MCmH[0]->Fill(UpsM,weight);
					
					MCrAllcutH[h]->Fill(UpsRap,weight);
					
					for (int j = 1; j<nPtBins; j++)
					{
						if (UpsPt < j)
						{
							MCmH[j]->Fill(UpsM,weight);
							break;
						}
					}
				}
				
				delete rfunc;
				delete mfunc;
				delete ptfunc;
			}
			
		}
	}
	for (int h = 0; h < nbinsM; h++)
		MCrH[h]->Scale((rH[h]->Integral())/(MCrH[h]->Integral()));
	for (int h = 0; h < nbinsM; h++)
		MCrAllcutH[h]->Scale((rAllcutH[h]->Integral())/(MCrAllcutH[h]->Integral()));
	for (int h = 0; h < nbinsM; h++)
		MCptH[h]->Scale((ptH[h]->Integral())/(MCptH[h]->Integral()));
	for (int h = 0; h < nPtBins; h++)
		MCmH[h]->Scale((mH[h]->Integral())/(MCmH[h]->Integral()));
	
	///////////Fit MC mass histograms
	/*TF1* mF1[nPtBins];
	double mParam[7][nPtBins];
	for (int h = 1; h < nPtBins; h++)
	{
		mF1[h] = new TF1(Form("mFit1_%d",h),fSumErfExp,minM,maxM,7);
		mF1[h]->SetParNames("Norm","Ratio","erf #mu","erf #sigma","exp #mu","exp #lambda","decay #lambda");
		mF1[h]->SetParameters(20,0.8,8,0.5,8,1,3.4);
		mF1[h]->SetParLimits(2,6,10);
		mF1[h]->SetParLimits(3,0,20);
		mF1[h]->SetParLimits(4,6,10);
		mF1[h]->SetParLimits(5,0,20);
		mF1[h]->SetParLimits(6,0,100);
		MCmH[h]->Fit(mF1[h]);
		
		TF1* hf = MCmH[h]->GetFunction(Form("mFit1_%d", h));
		for (int p = 0; p < 7; p++)
		{
			param[p][h] = hf->GetParameter(p);
		}
	}
	mF1[0] = new TF1("mFit1_0",fSumErfExpTotal2,minM,maxM,40);
	for (int h = 1; h < nPtBins; h++)
	{
		for (int p = 0; p < 7; p++)
		{
			mF1[0]->FixParameter(5+(h-1)*7+p,mParam[p][h]);
		}
	}
	mF1[0]->SetParLimits(0,0,1000);
	mF1[0]->SetParLimits(1,0,1000);
	mF1[0]->SetParLimits(2,0,1000);
	mF1[0]->SetParLimits(3,0,1000);
	mF1[0]->SetParLimits(4,0,1000);
	MCmH[0]->Fit(mF1[0]);
	*/
	
	
	///////////Draw MC hists
	//Draw MC rapidity hists
	for (int h = hFirst, i = 1; h <= hLast; h++, i++)
	{
		c1->cd(i);
		
		MCrH[h]->SetLineColor(kGreen+1);
		MCrH[h]->SetMinimum(0);
		MCrH[h]->Draw("SAME");
		/*MCrH[h]->GetXaxis()->SetTitleSize(0.1);
		MCrH[h]->GetXaxis()->SetLabelSize(0.09);
		MCrH[h]->GetYaxis()->SetTitleSize(0.1);
		MCrH[h]->GetYaxis()->SetLabelSize(0.09);*/
		MCrH[h]->SetXTitle("Rapidity"); 
		MCrH[h]->SetYTitle("count");
		MCrH[h]->SetStats(0);
		//gStyle->SetOptStat(1000000001);
		//gStyle->SetOptFit(1);
		gStyle->SetTitleFontSize(.1);
		c1->Update();
	}
	for (int h = hFirst, i = 1; h <= hLast; h++, i++)
	{
		crac->cd(i);
		
		MCrAllcutH[h]->SetLineColor(kGreen+1);
		MCrAllcutH[h]->SetMinimum(0);
		MCrAllcutH[h]->Draw("SAME");
		MCrAllcutH[h]->SetXTitle("Rapidity"); 
		MCrAllcutH[h]->SetYTitle("count");
		/*MCrAllcutH[h]->GetXaxis()->SetTitleSize(0.1);
		MCrAllcutH[h]->GetXaxis()->SetLabelSize(0.09);
		MCrAllcutH[h]->GetYaxis()->SetTitleSize(0.1);
		MCrAllcutH[h]->GetYaxis()->SetLabelSize(0.09);*/
		MCrAllcutH[h]->SetStats(0);
		gStyle->SetTitleFontSize(.1);
		crac->Update();
	}
	//Draw MC pt hists
	for (int h = hFirst, i = 1; h <= hLast; h++, i++)
	{
		cpt->cd(i);
		
		MCptH[h]->SetLineColor(kGreen+1);
		MCptH[h]->SetMinimum(0);
		MCptH[h]->Draw("SAME");
		/*MCptH[h]->GetXaxis()->SetTitleSize(0.1);
		MCptH[h]->GetXaxis()->SetLabelSize(0.09);
		MCptH[h]->GetYaxis()->SetTitleSize(0.1);
		MCptH[h]->GetYaxis()->SetLabelSize(0.09);*/
		MCptH[h]->SetXTitle("Rapidity"); 
		MCptH[h]->SetYTitle("count");
		MCptH[h]->SetStats(0);
		cpt->Update();
	}
	//Draw MC mass hists
	for (int h = 0, i = 1; h <= 5; h++, i++)
	{
		cM->cd(i);
		
		MCmH[h]->SetLineColor(kGreen+1);
		MCmH[h]->SetMinimum(0);
		MCmH[h]->Draw("SAME");
		/*MCmH[h]->GetXaxis()->SetTitleSize(0.1);
		MCmH[h]->GetXaxis()->SetLabelSize(0.09);
		MCmH[h]->GetYaxis()->SetTitleSize(0.1);
		MCmH[h]->GetYaxis()->SetLabelSize(0.09);*/
		MCmH[h]->SetXTitle("Mass (GeV)"); 
		MCmH[h]->SetYTitle("count");
		//MCmH[h]->SetStats(0);
		gStyle->SetOptStat(1000000001);
		gStyle->SetOptFit(1);
		cM->Update();
	}
	
	
	
}