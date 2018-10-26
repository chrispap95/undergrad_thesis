void delphes_routine_1() {
  TFile* f1 = TFile::Open("/data/cpapage/higgsTojpsiGamma-sim8x.root");
  TTree* t1 = (TTree*)(f1->Get("Delphes"));

  TCanvas* c1 = new TCanvas("c1","Delphes Results",1);
  TH1F* h12 = new TH1F("h12", "p_{T} of #mu; p_{T} (GeV/c)", 100, 0, 200);
  TH1F* h22 = new TH1F("h22", "#phi of #mu; #phi", 100, -4, 4);
  TH1F* h32 = new TH1F("h32", "#eta of #mu; #eta", 100, -3, 3);
  //TH1F* h42 = new TH1F("h42", "E of #mu; E(GeV)", 100, 0, 600);
  TH1F* h13 = new TH1F("h13", "p_{T} of #gamma; p_{T} (GeV/c)", 100, 0, 400);
  TH1F* h23 = new TH1F("h23", "#phi of #gamma; #phi", 100, -4, 4);
  TH1F* h33 = new TH1F("h33", "#eta of #gamma; #eta", 100, -3, 3);
  TH1F* h43 = new TH1F("h43", "E of #gamma; E (GeV)", 100, 0, 1000);

  c1->Divide(2,2);
  c1->cd(1); t1->Draw("Muon.PT>>h12");
  c1->cd(2); t1->Draw("Muon.Phi>>h22");
  c1->cd(3); t1->Draw("Muon.Eta>>h32");
  //c1->cd(4); t1->Draw("Particle.PT>>h42", "Particle.PID == 13");
  c1->Print("simulation_histograms.pdf(");
  c1->cd(1); t1->Draw("Photon.PT>>h13");
  c1->cd(2); t1->Draw("Photon.Phi>>h23");
  c1->cd(3); t1->Draw("Photon.Eta>>h33");
  c1->cd(4); t1->Draw("Photon.E>>h43");
  c1->Print("simulation_histograms.pdf)");
}
