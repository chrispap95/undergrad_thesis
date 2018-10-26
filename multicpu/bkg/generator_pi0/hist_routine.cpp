{
	TFile* f = TFile::Open("hisol02.root");
	TH1D* h03s = (TH1D*)(f->Get("h03s"));
	TH1D* h04s = (TH1D*)(f->Get("h04s"));
	TH1D* h03b = (TH1D*)(f->Get("h03b"));
	TH1D* h04b = (TH1D*)(f->Get("h04b"));

	h03b->Rebin(15);
	h04b->Rebin(15);
	h03s->GetXaxis()->SetRangeUser(0,3);
	h04s->GetXaxis()->SetRangeUser(0,0.2);
	h03b->SetTitle("Isolation of #gamma for background");
	h03s->SetTitle("Isolation of #gamma for signal");
	h04b->SetTitle("Relative Isolation of #gamma for background");
	h04s->SetTitle("Relative Isolation of #gamma for signal");

	TCanvas* c = new TCanvas("c","c",1);
	c->Divide(2,2);
	c->cd(1); h03b->Draw();
	c->cd(2); h04b->Draw();
	c->cd(3); h03s->Draw();
	c->cd(4); h04s->Draw();
	c->Print("isolation_plots02.eps");
}
