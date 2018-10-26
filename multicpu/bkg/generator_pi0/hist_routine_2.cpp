{
	TFile* f = TFile::Open("hisol02.root");
	TH1D* h04s = (TH1D*)(f->Get("h04s"));
	TH1D* h04b = (TH1D*)(f->Get("h04b"));

	h04s->SetTitle("Relative Isolation; rel. I");
	h04b->SetTitle("Relative Isolation; rel. I");

	h04s->GetXaxis()->SetRangeUser(0,0.3);
	h04b->GetXaxis()->SetRangeUser(0,0.3);

	gStyle->SetOptStat(0000);
	gStyle->SetPalette(kBird);
	TCanvas* c = new TCanvas("c","c",1);
	h04s->Draw("hist plc");
	h04s->GetYaxis()->SetTitle("Signal");
	c->Update();
	Float_t rightmax = 1.1*h04b->GetMaximum();
	Float_t scale = gPad->GetUymax()/rightmax;
	h04b->Scale(scale);
	h04b->Draw("hist same plc");
	TGaxis* axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymax(),0,rightmax,510,"+L");
	axis->SetTextColor(kYellow+2);
	axis->SetLineColor(kYellow);
	axis->SetLabelColor(kYellow+1);
	axis->SetLabelSize(0.033);
	axis->SetTextSize(0.033);
	axis->SetTitle("Background");
	axis->Draw();
	TLine* line = new TLine();
	line->SetLineStyle(5);
	line->SetLineColor(kRed);
	line->DrawLine(0.08,gPad->GetUymin(),0.08,gPad->GetUymax());
	TLatex text;
	text.SetTextColor(kRed);
	text.SetTextSize(0.035);
	text.DrawLatex(0.085,3500,"Rel. Isolation cut");
	c->Print("isolation_plots02_2.eps");
}
