{
	TFile* f = TFile::Open("photon_cut.root");
	TH1D* hcounts = (TH1D*)(f->Get("hcounts"));
	TH1D* hcountb = (TH1D*)(f->Get("hcountb"));

	TH1D* hsb = new TH1D("hsb","S/B",100,20,60);
	TH1D* hsbsq = new TH1D("hsbsq","S/#sqrt{B}",100,20,60);

	Int_t countb_norm[101],countb[101];
	Int_t counts_norm[101],counts[101];

	for (Int_t i = 1; i<101; ++i) {
		counts[i] = hcounts->GetBinContent(i);
		countb[i] = hcountb->GetBinContent(i);
	}

	for (Int_t i = 1; i<101; ++i) {
		countb_norm[i] = countb[i];
		counts_norm[i] = counts[i];
	}
	float normb = 169800*0.0003;
	float norms = 0.027306;
	for (Int_t i = 1; i<101; ++i) {
		countb_norm[i] *= normb;
		counts_norm[i] *= norms;
	}

	for (Int_t i = 1; i<101; ++i) {
		if (countb_norm[i] >= 1) {
			hsb->SetBinContent(i, (double)counts_norm[i]/(double)countb_norm[i]);
			hsbsq->SetBinContent(i, counts_norm[i]/sqrt(countb_norm[i]));
		}
		else {
			hsb->SetBinContent(i, 0);
			hsbsq->SetBinContent(i, 0);
		}
	}

	TCanvas* c = new TCanvas("c","c",1);
	c->Divide(2,2);
	c->cd(1); hcounts->Draw();
	c->cd(2); hcountb->Draw();
	c->cd(3); hsb->Draw();
	c->cd(4); hsbsq->Draw();
}
