{
	TFile* f = TFile::Open("pt_ratio_high_new.root");

	TCanvas* c1 = new TCanvas("c1","c1",1);

	TH1F* h1 = (TH1F*)(f->Get("h1"));
	TH1F* h2 = (TH1F*)(f->Get("h2"));
	h1->Sumw2();
	h2->Sumw2();

	gStyle->SetOptStat(0);
	h1->Rebin(5);
	h2->Rebin(5);
	h1->SetLineColor(kYellow);
       	h1->SetTitle("p_{T} of jets in p_{T} > 70 GeV range");
	h2->SetTitle("p_{T} of jets in p_{T} > 70 GeV range");

	h1->GetYaxis()->SetRangeUser(0.5,2*h1->GetBinContent(h1->GetMaximumBin()));
	h2->GetYaxis()->SetRangeUser(0.5,2*h1->GetBinContent(h1->GetMaximumBin()));

	auto rp = new TRatioPlot(h2,h1,"divsym");
	rp->Draw();
	rp->GetUpperPad()->SetLogy();
	//rp->GetLowerPad()->SetLogy();
	rp->SetH1DrawOpt("hist plc");
	rp->SetH2DrawOpt("hist plc");
	rp->GetUpperRefYaxis()->SetTitle("Entries / 6.5 GeV");
	rp->GetLowerRefYaxis()->SetTitle("ratio");
	rp->GetLowerRefYaxis()->SetTitleOffset(1.57);
	rp->GetUpperRefXaxis()->SetTitle("p_{T} [GeV]");
	rp->SetGridlines(0,0);

	rp->GetCalculationOutputGraph()->Fit("pol0","","",70,200);
	TFitResultPtr fitptr = rp->GetCalculationOutputGraph()->Fit("pol0","S","",70,200);
	Double_t chi2 = fitptr->Chi2();                  // to retrieve the fit chi2
	Double_t par0 = fitptr->Parameter(0);            // retrieve the value for the parameter 0
	Double_t err0 = fitptr->ParError(0);             // retrieve the error for the parameter 0
	Int_t ndf = fitptr->Ndf();
	TLatex text;
	string string1 = "p0 = " + to_string(par0) + "#pm" + to_string(err0);
	string string2 = "chi2 = " + to_string(chi2);
	string string3 = "NDf = " + to_string(ndf);

	rp->GetLowerPad()->cd();
	text.SetTextSize(0.11);
	text.DrawLatex(142.75,0.000625,string1.c_str());
	text.DrawLatex(142.75,0.000475,string2.c_str());
	text.DrawLatex(177,0.000475,string3.c_str());
	c1->Update();
	c1->Print("generator_jet_ratio_high_pt120_200K.eps");

	// c1->Print("generator_jet_ratio_low_pt40_200K_4.pdf");
	// ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
	// TF1* f1 = new TF1("f1","[0]+[1]*x",30,50);
	// f1->SetParameters(-0.1,0);
	// f1->SetParLimits(0,-0.1,0);
	// f1->SetParLimits(1,0,0.001);
	// rp->GetCalculationOutputGraph()->Fit("f1","","",30,50);
	// c1->Update();
	// c1->Print("generator_jet_ratio_low_pt40_200K_3.pdf(");

	// // Pull Plot
	// Double_t* x;
	// Double_t* y;
	// Double_t* err;
	// x = rp->GetCalculationOutputGraph()->GetX();
	// y = rp->GetCalculationOutputGraph()->GetY();
	// err = rp->GetCalculationOutputGraph()->GetEY();
	// for (Int_t i = 1; i <= 25; ++i) {
	// 	Float_t res = y[i] - f1->Eval(x[i]);
	// 	cout << res << '\t' << err[i] << endl;
	// 	hpull->Fill(res/err[i]);
	// }
	// c1->Clear();
	// gStyle->SetOptFit(1011);
	// hpull->Draw();
	// hpull->Fit("gaus");
	// c1->Print("generator_jet_ratio_low_pt40_200K_3.pdf)");
}
