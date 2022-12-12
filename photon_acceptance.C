 void photon_acceptance(TString mcPeriod = "LHC16k")
 {

   std::vector<int> runNumbers; 
   ifstream ff;
   
//   ff.open(Form("datasets_MC/runs_%s.dat", mcPeriod.Data()));
   ff.open(Form("datasets/LHC16k_good.txt"));
   Int_t nn;
   while (ff >> nn) runNumbers.emplace_back(nn);
   
   const Int_t nRuns = runNumbers.size();
   auto histAcc = new TH1F("histAcc", "PHOS acceptance per run", nRuns, 0, nRuns*1.0);
   
   TFile *ftemp;
   Double_t xb[2] = {0.3, 20};
   Int_t iBin;

   for (auto run : runNumbers) {
     iBin ++;

     printf("%d %d\n", iBin, run);

     ftemp = TFile::Open(Form("../MergedResults_MC/%s/%d/AnalysisResults.root", mcPeriod.Data(), run));

     if (!ftemp) continue;

     auto dd  = (THashList*)ftemp->Get("Data");
     auto dd2 = (THashList*)ftemp->Get("Data2");

     auto htot  = (TH2F*)dd2->FindObject("hGammaMC_true");
     auto hpass = (TH2F*)dd2->FindObject("hCaloPhotonPdgvsPt_all");

     auto htot1  = (TH1D*)htot ->ProjectionX("htot1", 71, 170);
     auto hpass1 = (TH1D*)hpass->ProjectionX("hpass1", 4022, 4023);
    
     auto ht1 = (TH1D*)htot1 ->Rebin(1, "ht1", xb);
     auto hp1 = (TH1D*)hpass1->Rebin(1, "hp1", xb);

     hp1->Divide(hp1, ht1, 1, 1, "b");
     
     histAcc->SetBinContent(iBin + 1, hp1->GetBinContent(1));
     histAcc->SetBinError(iBin + 1, hp1->GetBinError(1));

     delete dd;
     delete dd2;
     delete ftemp;
   } 
   
   auto fout = TFile::Open(Form("ROOT/Efficiency/Acceptance_%s.root", mcPeriod.Data()), "recreate");
   
   printf("\n saving output to file %s\n", fout->GetName());
   
   histAcc->Write();

   fout->Close();

 }
