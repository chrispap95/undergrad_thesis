{
	TFile* f = TFile::Open("ptratio_vs_ptjet_low.root");
	TH2F* h1 = (TH2F*)(f->Get("h2d"));
	TFile* f2 = TFile::Open("ptratio_vs_ptjet_mid.root");
	TH2F* h2 = (TH2F*)(f2->Get("h2d"));
	TFile* f3 = TFile::Open("ptratio_vs_ptjet_high.root");
	TH2F* h3 = (TH2F*)(f3->Get("h2d"));
	TCanvas* c = new TCanvas("c","c",1);

	h3->Add(h1);
	h3->Add(h2);
	h3->RebinX(5);
	h3->RebinY(5);
	gStyle->SetOptStat(0);
	h3->Draw("colz");
	h3->SetTitle("p_{T}^{isol. #pi^{0}}/p_{T}^{jet} vs p_{T}^{jet}");
 	c->Print("ptratio_vs_ptjet.eps");
	TLatex l;

	c->Clear();
	c->Divide(2,2);

	c->cd(1);
	TCutG* cutg1 = new TCutG("cutg1",4);
	cutg1->SetVarX("x");
	cutg1->SetVarY("y");
	cutg1->SetPoint(0,30,0);
	cutg1->SetPoint(1,30,1);
	cutg1->SetPoint(2,50,0);
	cutg1->SetPoint(3,50,1);
	TH1D* hpy1 = h3->ProjectionY("hpy1",0,-1,"[cutg1]");
	hpy1->SetTitle("p_{T}^{#gamma}/p_{T}^{Jet} projection");
	hpy1->GetXaxis()->SetTitleOffset(1);
	hpy1->Rebin();
	hpy1->Draw(); 
	l.DrawLatex(0.5,20,"p_{T}^{Jet} in [30,50] GeV");

	c->cd(2);
	TCutG* cutg2 = new TCutG("cutg2",4);
	cutg2->SetVarX("x");
	cutg2->SetVarY("y");
	cutg2->SetPoint(0,50,0);
	cutg2->SetPoint(1,50,1);
	cutg2->SetPoint(2,70,0);
	cutg2->SetPoint(3,70,1);
	TH1D* hpy2 = h3->ProjectionY("hpy2",0,-1,"[cutg2]");
	hpy2->SetTitle("p_{T}^{#gamma}/p_{T}^{Jet} projection");
	hpy2->GetXaxis()->SetTitleOffset(1);
	hpy2->Rebin();
	hpy2->Draw();
	l.DrawLatex(0.5,10,"p_{T}^{Jet} in [50,70] GeV");

	c->cd(3);
	TCutG* cutg3 = new TCutG("cutg3",4);
	cutg3->SetVarX("x");
	cutg3->SetVarY("y");
	cutg3->SetPoint(0,70,0);
	cutg3->SetPoint(1,70,1);
	cutg3->SetPoint(2,90,0);
	cutg3->SetPoint(3,90,1);
	TH1D* hpy3 = h3->ProjectionY("hpy3",0,-1,"[cutg3]");
	hpy3->SetTitle("p_{T}^{#gamma}/p_{T}^{Jet} projection");
	hpy3->GetXaxis()->SetTitleOffset(1);
	hpy3->Rebin();
	hpy3->Draw();
	l.DrawLatex(0.55,5,"p_{T}^{Jet} in [70,90] GeV");

	c->cd(4);
	TCutG* cutg4 = new TCutG("cutg4",4);
	cutg4->SetVarX("x");
	cutg4->SetVarY("y");
	cutg4->SetPoint(0,90,0);
	cutg4->SetPoint(1,90,1);
	cutg4->SetPoint(2,110,0);
	cutg4->SetPoint(3,110,1);
	TH1D* hpy4 = h3->ProjectionY("hpy4",0,-1,"[cutg4]");
	hpy4->SetTitle("p_{T}^{#gamma}/p_{T}^{Jet} projection");
	hpy4->GetXaxis()->SetTitleOffset(1);
	hpy4->Rebin();
	hpy4->Draw();
	l.DrawLatex(0.4,1,"p_{T}^{Jet} in [90,110] GeV");
	c->Print("ptratio_projections.eps");

	TFile* fout = new TFile("ptratio_distributions.root","RECREATE");
	hpy1->Write();
	hpy2->Write();
	hpy3->Write();
	hpy4->Write();
	fout->Close();
}
