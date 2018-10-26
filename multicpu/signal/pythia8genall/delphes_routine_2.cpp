void delphes_routine_1() {
  TFile* f1 = TFile::Open("/data/cpapage/higgsTojpsiGamma-sim8x.root");
  TTree* t1 = (TTree*)(f1->Get("Delphes"));

  TCanvas* c1 = new TCanvas("c1","Delphes Results",1);
  TH1F* h12 = new TH1F("h12", "p_{T} of #mu min; p_{T} (GeV/c)", 100, 0, 200);
  TH1F* h22 = new TH1F("h22", "#phi of #mu min; #phi", 100, -4, 4);
  TH1F* h32 = new TH1F("h32", "#eta of #mu min; #eta", 100, -3, 3);
  TH1F* h42 = new TH1F("h42", "p_{T} of #mu max; p_{T} (GeV/c)", 100, 0, 200);
  TH1F* h52 = new TH1F("h52", "#phi of #mu max; #phi", 100, -4, 4);
  TH1F* h62 = new TH1F("h62", "#eta of #mu max; #eta", 100, -3, 3);
  TH1F* h72 = new TH1F("h12", "p_{T} of #mu min; p_{T} (GeV/c)", 100, 0, 200);
  TH1F* h82 = new TH1F("h22", "#phi of #mu min; #phi", 100, -4, 4);
  TH1F* h92 = new TH1F("h32", "#eta of #mu min; #eta", 100, -3, 3);
  TH1F* h13 = new TH1F("h13", "p_{T} of #gamma; p_{T} (GeV/c)", 100, 0, 400);
  TH1F* h23 = new TH1F("h23", "#phi of #gamma; #phi", 100, -4, 4);
  TH1F* h33 = new TH1F("h33", "#eta of #gamma; #eta", 100, -3, 3);
  TH1F* h43 = new TH1F("h43", "E of #gamma; E (GeV)", 100, 0, 1000);

  THStack* hs12 = new THStack("hs12", "p_{T} of #mu min; p_{T} (GeV/c)");
  THStack* hs22 = new THStack("h22", "#phi of #mu min; #phi");
  THStack* hs32 = new THStack("h32", "#eta of #mu min; #eta");

  c1->Divide(2,2);
  c1->cd(1); t1->Draw("Muon.PT>>h12");
  c1->cd(2); t1->Draw("Muon.Phi>>h22");
  c1->cd(3); t1->Draw("Muon.Eta>>h32");
  c1->Print("simulation_signal2.pdf(");
  c1->cd(1); t1->Draw("Photon.PT>>h13");
  c1->cd(2); t1->Draw("Photon.Phi>>h23");
  c1->cd(3); t1->Draw("Photon.Eta>>h33");
  c1->cd(4); t1->Draw("Photon.E>>h43");
  c1->Print("simulation_signal2.pdf)");

	TFile* f2 = new TFile("signal_sim_hists.root", "RECREATE");
	h12->Write();
	h22->Write();
	h32->Write();
	h13->Write();
	h23->Write();
	h33->Write();
	h43->Write();
	f2->Close();
}
