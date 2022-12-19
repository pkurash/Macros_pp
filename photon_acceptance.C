 void photon_acceptance(TString mcPeriod = "LHC16k")
 {

   std::vector<int> runNumbers; 
   ifstream ff;
   
   ff.open(Form("datasets/%s_good.txt", mcPeriod.Data()));
//   ff.open(Form("datasets/LHC16k_good.txt"));
   Int_t nn;
   while (ff >> nn) runNumbers.emplace_back(nn);
   
   const Int_t nRuns = runNumbers.size();
   auto histAcc =  new TH1F("histAcc",  "fraction of detected #gamma-s E > 0.3", nRuns, 0, nRuns*1.0);
   auto histAcc2 = new TH1F("histAcc2", "fraction of detected #gamma-s E > 1.0", nRuns, 0, nRuns*1.0);
   auto histAcc3 = new TH1F("histAcc3", "fraction of detected #gamma-s E > 2.0", nRuns, 0, nRuns*1.0);
   auto histAcc4 = new TH1F("histAcc4", "fraction of detected #gamma-s E > 4.0", nRuns, 0, nRuns*1.0);
   auto histAcc5 = new TH1F("histAcc5", "fraction of detected #gamma-s E > 7.0", nRuns, 0, nRuns*1.0);
   
   TFile *ftemp;
   Double_t xb[2]  = {0.3, 20};
   Double_t xb2[2] = {1.0, 20};
   Double_t xb3[2] = {2.0, 20};
   Double_t xb4[2] = {4.0, 20};
   Double_t xb5[2] = {7.0, 20};

   Int_t iBin = 0; 

   for (auto run : runNumbers) {
     iBin ++;

     printf("%d %d\n", iBin, run);

     ftemp = TFile::Open(Form("../MergedResults_MC/%s/%d/AnalysisResults.root", mcPeriod.Data(), run));

     if (!ftemp) {
       histAcc->SetBinContent(iBin + 1, 0);
       histAcc->SetBinError(iBin + 1, 0);
       continue;
     }

     auto dd  = (THashList*)ftemp->Get("Data");
     auto dd2 = (THashList*)ftemp->Get("Data2");

     auto htot  = (TH2F*)dd2->FindObject("hGammaMC_true");
     auto htot_mc  = (TH1D*)htot ->ProjectionX("htot1", 71, 170);
     auto ht1  = (TH1D*)htot_mc ->Rebin(1, "ht1", xb);

     auto hpass = (TH2F*)dd2->FindObject("hCaloPhotonPdgvsPt_all");
     auto hpass_mc = (TH1D*)hpass->ProjectionX("hpass1", 4022, 4023);
     auto hp1 = (TH1D*)hpass_mc->Rebin(1, "hp1", xb);
     auto hp2 = (TH1D*)hpass_mc->Rebin(1, "hp2", xb2);
     auto ht2 = (TH1D*)htot_mc ->Rebin(1, "ht2", xb2);
     auto hp3 = (TH1D*)hpass_mc->Rebin(1, "hp3", xb3);
     auto ht3 = (TH1D*)htot_mc ->Rebin(1, "ht3", xb3);
     auto hp4 = (TH1D*)hpass_mc->Rebin(1, "hp4", xb4);
     auto ht4 = (TH1D*)htot_mc ->Rebin(1, "ht4", xb4);
     auto hp5 = (TH1D*)hpass_mc->Rebin(1, "hp5", xb5);
     auto ht5 = (TH1D*)htot_mc ->Rebin(1, "ht5", xb5);

     hp1->Divide(hp1, ht1, 1, 1, "b");
     hp2->Divide(hp2, ht2, 1, 1, "b");
     hp3->Divide(hp3, ht3, 1, 1, "b");
     hp4->Divide(hp4, ht4, 1, 1, "b");
     hp5->Divide(hp5, ht5, 1, 1, "b");
     
     histAcc ->SetBinContent(iBin, hp1->GetBinContent(1));
     histAcc ->SetBinError(iBin, hp1->GetBinError(1));

     histAcc2->SetBinContent(iBin, hp2->GetBinContent(1));
     histAcc2->SetBinError(iBin, hp2->GetBinError(1));
     histAcc3->SetBinContent(iBin, hp3->GetBinContent(1));
     histAcc3->SetBinError(iBin, hp3->GetBinError(1));
     histAcc4->SetBinContent(iBin, hp4->GetBinContent(1));
     histAcc4->SetBinError(iBin, hp4->GetBinError(1));
     histAcc5->SetBinContent(iBin, hp5->GetBinContent(1));
     histAcc5->SetBinError(iBin, hp5->GetBinError(1));

     delete dd;
     delete dd2;
     delete ftemp;
   } 

   auto cc = new TCanvas("cc", "cc", 700, 10*nRuns);
   histAcc->SetAxisRange(0., 0.055, "Y");
   histAcc->SetLineColor(kRed); histAcc->SetMarkerColor(histAcc->GetLineColor());
   histAcc2->SetLineColor(kBlue); histAcc2->SetMarkerColor(histAcc2->GetLineColor());
   histAcc3->SetLineColor(kGreen); histAcc3->SetMarkerColor(histAcc3->GetLineColor());
   histAcc4->SetLineColor(kMagenta); histAcc4->SetMarkerColor(histAcc4->GetLineColor());
   histAcc5->SetLineColor(kBlack); histAcc5->SetMarkerColor(histAcc5->GetLineColor());

   gStyle->SetOptTitle(0);
   gStyle->SetOptStat(0);
   cc->cd();
   histAcc->Draw();
   histAcc2->Draw("same");
   histAcc3->Draw("same");
   histAcc4->Draw("same");
   histAcc5->Draw("same");
   
   auto fout = TFile::Open(Form("ROOT/Efficiencies/Acceptance_%s.root", mcPeriod.Data()), "recreate");
   
   printf("\n saving output to file %s\n", fout->GetName());
   
   histAcc->Write();

   histAcc2->Write();
   histAcc3->Write();
   histAcc4->Write();
   histAcc5->Write();

   cc->Write();

   fout->Close();

 }
