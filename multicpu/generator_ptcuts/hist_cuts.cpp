{
	TFile* f = TFile::Open("cuts.root");
	TH1D* hcounts1 = (TH1D*)(f->Get("hcounts1"));
	TH1D* hcountb1 = (TH1D*)(f->Get("hcountb10"));
	TH1D* hcounts2 = (TH1D*)(f->Get("hcounts2"));
	TH1D* hcountb2 = (TH1D*)(f->Get("hcountb20"));
	TH1D* hcounts3 = (TH1D*)(f->Get("hcounts3"));
	TH1D* hcountb3 = (TH1D*)(f->Get("hcountb30"));
	TH1D* hcounts4 = (TH1D*)(f->Get("hcounts4"));
	TH1D* hcountb4 = (TH1D*)(f->Get("hcountb40"));

	hcounts1->GetXaxis()->SetTitle("p_{T} cut");
	hcounts2->GetXaxis()->SetTitle("p_{T} cut");
	hcounts3->GetXaxis()->SetTitle("p_{T} cut");
	hcounts4->GetXaxis()->SetTitle("p_{T} cut");
	hcountb1->GetXaxis()->SetTitle("p_{T} cut");
	hcountb2->GetXaxis()->SetTitle("p_{T} cut");
	hcountb3->GetXaxis()->SetTitle("p_{T} cut");
	hcountb4->GetXaxis()->SetTitle("p_{T} cut");

	TH1D* hsb1 = new TH1D("hsb1","S/B;p_{T} cut",100,20,60);
	TH1D* hsbsq1 = new TH1D("hsbsq1","S/#sqrt{B};p_{T} cut",100,20,60);
	TH1D* hsb2 = new TH1D("hsb2","S/B;p_{T} cut",100,0,30);
	TH1D* hsbsq2 = new TH1D("hsbsq2","S/#sqrt{B};p_{T} cut",100,0,30);
	TH1D* hsb3 = new TH1D("hsb3","S/B;p_{T} cut",100,0,20);
	TH1D* hsbsq3 = new TH1D("hsbsq3","S/#sqrt{B};p_{T} cut",100,0,20);
	TH1D* hsb4 = new TH1D("hsb4","S/B;p_{T} cut",100,0,50);
	TH1D* hsbsq4 = new TH1D("hsbsq4","S/#sqrt{B};p_{T} cut",100,0,50);

	TH1D* hestim1 = new TH1D("hestim1","#sqrt{2Nlog(1+S/B)-2S};p_{T} cut",100,20,60);
	TH1D* hestim2 = new TH1D("hestim2","#sqrt{2Nlog(1+S/B)-2S};p_{T} cut",100,0,30);
	TH1D* hestim3 = new TH1D("hestim3","#sqrt{2Nlog(1+S/B)-2S};p_{T} cut",100,0,20);
	TH1D* hestim4 = new TH1D("hestim4","#sqrt{2Nlog(1+S/B)-2S};p_{T} cut",100,0,50);

	Double_t countb1[101];
	Double_t counts1[101];
	Double_t countb2[101];
	Double_t counts2[101];
	Double_t countb3[101];
	Double_t counts3[101];
	Double_t countb4[101];
	Double_t counts4[101];

	for (Int_t i = 1; i<101; ++i) {
		counts1[i] = hcounts1->GetBinContent(i);
		countb1[i] = hcountb1->GetBinContent(i);
		counts2[i] = hcounts2->GetBinContent(i);
		countb2[i] = hcountb2->GetBinContent(i);
		counts3[i] = hcounts3->GetBinContent(i);
		countb3[i] = hcountb3->GetBinContent(i);
		counts4[i] = hcounts4->GetBinContent(i);
		countb4[i] = hcountb4->GetBinContent(i);
	}

	for (Int_t i = 1; i<101; ++i) {
		if (countb1[i] >= 0.5) {
			hsb1->SetBinContent(i, counts1[i]/countb1[i]);
			hsbsq1->SetBinContent(i, counts1[i]/sqrt(countb1[i]));
			hestim1->SetBinContent(i, sqrt(2*(counts1[i]+countb1[i])*log(1+counts1[i]/countb1[i])-2*counts1[i]));
		}
		else {
			hsb1->SetBinContent(i, 0);
			hsbsq1->SetBinContent(i, 0);
			hestim1->SetBinContent(i, 0);
		}
		if (countb2[i] >= 0.5) {
			hsb2->SetBinContent(i, counts2[i]/countb2[i]);
			hsbsq2->SetBinContent(i, counts2[i]/sqrt(countb2[i]));
			hestim2->SetBinContent(i, sqrt(2*(counts2[i]+countb2[i])*log(1+counts2[i]/countb2[i])-2*counts2[i]));
		}
		else {
			hsb2->SetBinContent(i, 0);
			hsbsq2->SetBinContent(i, 0);
			hsbsq2->SetBinContent(i, 0);
		}
		if (countb3[i] >= 0.5) {
			hsb3->SetBinContent(i, counts3[i]/countb3[i]);
			hsbsq3->SetBinContent(i, counts3[i]/sqrt(countb3[i]));
			hestim3->SetBinContent(i, sqrt(2*(counts3[i]+countb3[i])*log(1+counts3[i]/countb3[i])-2*counts3[i]));
		}
		else {
			hsb3->SetBinContent(i, 0);
			hsbsq3->SetBinContent(i, 0);
			hsbsq3->SetBinContent(i, 0);
		}
		if (countb4[i] >= 0.5) {
			hsb4->SetBinContent(i, counts4[i]/countb4[i]);
			hsbsq4->SetBinContent(i, counts4[i]/sqrt(countb4[i]));
			hestim4->SetBinContent(i, sqrt(2*(counts4[i]+countb4[i])*log(1+counts4[i]/countb4[i])-2*counts4[i]));
		}
		else {
			hsb4->SetBinContent(i, 0);
			hsbsq4->SetBinContent(i, 0);
			hsbsq4->SetBinContent(i, 0);
		}
	}

	TLine line1;
	line1.SetLineColor(kYellow);
        line1.SetLineWidth(2);
        line1.SetLineStyle(2);
	TLine line2;
	line2.SetLineColor(kRed);
        line2.SetLineWidth(2);

	TLatex text1, text2;
	text1.SetTextColor(kYellow+2); text2.SetTextColor(kRed);
	text1.SetTextSize(0.037); text2.SetTextSize(0.04);
	gStyle->SetOptStat(0);

	TCanvas* c = new TCanvas("c","c",1);
	c->Divide(2,2);
	c->cd(1); hcounts1->Draw();
	c->cd(2); hcountb1->Draw();
	c->cd(3); hsb1->Draw();
	c->cd(4); hsbsq1->Draw();
	line2.DrawLine(20,5,60,5); line1.DrawLine(57.5,0,57.5,9.7);
	text2.DrawLatex(19.05,4.7,"5"); text1.DrawLatex(56.7,-0.65,"57.5");
	c->Print("generator_cuts_1.eps");
	c->cd(1); hcounts2->Draw();
	c->cd(2); hcountb2->Draw();
	c->cd(3); hsb2->Draw();
	c->cd(4); hsbsq2->Draw();
	line1.DrawLine(15,9.2,15,13.5);
	text1.DrawLatex(14.3,-0.65,"15");
	c->Print("generator_cuts_2.eps");
	c->cd(1); hcounts2->Draw();
	c->cd(2); hcountb3->Draw();
	c->cd(3); hsb3->Draw();
	c->cd(4); hsbsq3->Draw();
	line1.DrawLine(57.5,0,57.5,9.7);
	text1.DrawLatex(56.7,-0.65,"57.5");
	c->Print("generator_cuts_3.eps");
	c->cd(1); hcounts4->Draw();
	c->cd(2); hcountb4->Draw();
	c->cd(3); hsb4->Draw();
	c->cd(4); hsbsq4->Draw();
	line1.DrawLine(28,8.4,28,13.1);
	text1.DrawLatex(27.3,-0.65,"28");
	c->Print("generator_cuts_4.eps");
	c->cd(1); hestim1->Draw();
	c->cd(2); hestim2->Draw();
	c->cd(3); hestim3->Draw();
	c->cd(4); hestim4->Draw();
	c->Print("generator_cuts_5.eps");
}
