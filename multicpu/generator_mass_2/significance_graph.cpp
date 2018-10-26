{
	gStyle->SetOptTitle(kFALSE);
	gStyle->SetPalette(kSolar);
	
	TCanvas *c = new TCanvas("c","c",1);
	Double_t x[340],y[340],z[340];
	Int_t n = 340;
	Int_t estimation = 0;
	for (Int_t i = 0; i < n; ++i) {
		x[i] = i*0.01+0.1;
		y[i] = 9.27437*sqrt(x[i]/3.0);
		z[i] = 5.22222*sqrt(x[i]/3.0);
		if(abs(z[i]-5) < 0.01) estimation = i; 
	}

	TGraph* gr1 = new TGraph(n,x,y);
	TGraph* gr2 = new TGraph(n,x,z);
	
	gr1->SetTitle("S_{1}");
	gr1->GetYaxis()->SetTitle("Significance [#sigma]");
	gr1->GetXaxis()->SetTitle("#int L [#alpha b^{-1}]");
	gr1->GetXaxis()->SetTitleOffset(1.2);
	gr1->GetXaxis()->SetRangeUser(0,3.5);
	gr1->SetLineWidth(3);

	gr2->SetTitle("S_{cL}");
	gr2->SetLineWidth(3);

	TLine gr3;
	TLine gr4;
	gr3.SetLineColor(kYellow);
	gr3.SetLineWidth(2);
	gr3.SetLineStyle(2);
	gr4.SetLineColor(kRed);
	gr4.SetLineWidth(2);
	gr4.SetLineStyle(2);

	gr1->Draw("AC plc pfc");
	gr2->Draw("C plc pfc same");

	gr3.DrawLine(0,3,3.5,3);
	gr4.DrawLine(0,5,3.5,5);

	if (estimation != 0) {
		TLine gr5;
		gr5.SetLineColor(kGreen);
		gr5.SetLineWidth(2);
		gr5.SetLineStyle(2);
		gr5.DrawLine(x[estimation],0.8,x[estimation],z[estimation]);
		TLatex text;
		text.SetTextSize(0.035);
		text.SetTextColor(kGreen+2);
		text.DrawLatex(0.01*estimation,0.45,"2.76");
	}

	cout << x[estimation] << endl;
	auto leg = new TLegend(0.5,8,1.5,10,"","NB"); 
	leg->AddEntry(gr1,"S_{1}","f");
	leg->AddEntry(gr2,"S_{cL}","f");
	leg->SetBorderSize(0);
	leg->Draw();
	c->Print("significance_graph.eps");
}
