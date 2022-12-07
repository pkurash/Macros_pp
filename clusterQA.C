/* 
 * Macro for filling out QA data for extracting 
 * 1) number of analyzed events,
 * 2) average number of PHOS calophotons per event in a given run, 
 * 3) mean energy of PHOS calophotons in a run
 * 4) mean number of PHOS celss per calophoton per run
 *
 */
void clusterQA(TString period = "LHC16g", 
                TString cut = "all", 
	        Double_t emin = 0.3, 
	        Int_t nCellMin = 3)
{
   std::vector<Int_t> runNumbers;
   Int_t nn;

   ifstream ff;

    printf("%s-%s.txt\n", period.Data());
    ff.open(Form("../../PHOSGamma_PP/datasets/%s.txt", period.Data()));
    while ( ff >> nn) { 
      cout << nn << endl;
      runNumbers.emplace_back(nn);
    }
    ff.close();
    cout << runNumbers.size() << " runs" << endl;

    const Int_t Nruns = runNumbers.size();

    TFile *f1;
    TList *Data;

    TH2F *hClusterEvsN[5];
    TH1F *hEn[5], *hNClust[5], *hNCell[5], *hNClust2[5], *hClustPt[5], *hcounter;
    Double_t xa[2]={emin, 20.};

    TH1F *hEvents = new TH1F("hEvents", "number of analyzed events;;N_{ev}", Nruns, 0, Nruns);

    for (Int_t iHist = 0; iHist < 5; iHist ++ ) {
      hEn[iHist] = new TH1F(Form("hEn%d", iHist), Form("%s;;<E_{cl}>, GeV/c", iHist == 0 ? "all modules" : Form("module %d", iHist)), Nruns, 0, Nruns); 
      hNClust[iHist] = new TH1F(Form("hNClust%d", iHist), Form("%s;;N_{cl}", iHist == 0 ? "all modules" : Form("module %d", iHist)), Nruns, 0, Nruns); 
      hNClust2[iHist] = new TH1F(Form("hNClust2%d", iHist), Form("%s;;N_{cl}", iHist == 0 ? "all modules" : Form("module %d", iHist)), Nruns, 0, Nruns); 
      hNCell[iHist] = new TH1F(Form("hNCell%d", iHist),Form("%s;;N_{cell}", iHist == 0 ? "all modules" : Form("module %d", iHist)), Nruns, 0, Nruns); 
    }

    Int_t iFile = 0;
 
    for (auto run : runNumbers) {
      
      printf("run: %d\n", run);
      
      f1 = TFile::Open(Form("../MergedResults_Run2/%s/%d/AnalysisResults.root", period.Data(), run));    
      if (!f1) { 
        printf("ERROR: no file for run %d found!\n", run);
	cout << Form("../MergedResults_Run2/%s/%d/AnalysisResults.root", period.Data(), run) << endl;
        return;
      }

      Data = (TList*)f1->Get("Data");
 
      hcounter =       (TH1F*)Data->FindObject("hSelEvents");
      Int_t nev =hcounter->GetBinContent(9);
      if (nev == 0) {
         printf("ERROR: no events in the file for run %d!\n", run);
	 return;
      }

      hEvents->SetBinContent(iFile + 1, nev);
      hEvents->GetXaxis()->SetBinLabel(iFile + 1, Form("%d", run));

      for (Int_t iHist = 0; iHist < 5; iHist ++) {
          hClusterEvsN[iHist]   = (TH2F*)Data->FindObject(Form("hClusterEvsN_%s%s", 
	                                       cut.Data(), iHist == 0 ? "" : Form("_M%d", iHist)));
          hClusterEvsN[iHist]->GetXaxis()->SetRangeUser(emin, 0.01*hClusterEvsN[iHist]->GetNbinsX());
          hClusterEvsN[iHist]->GetYaxis()->SetRange(nCellMin, hClusterEvsN[iHist]->GetNbinsY());
          
	  hClustPt[iHist] = (TH1F*)Data->FindObject(Form("hClustPt_%s%s", cut.Data(),
	                                       iHist == 0 ? "" : Form("_M%d", iHist)));
          TH1F *hpx = (TH1F*)hClustPt[iHist]->Rebin(1, "h", xa);
	  hpx->Scale(1./nev);

	  hEn[iHist]->SetBinContent(iFile + 1, hClusterEvsN[iHist]->GetMean(1));
          hEn[iHist]->SetBinError(iFile + 1, 
	                  hClusterEvsN[iHist]->GetMean(1)/TMath::Sqrt(hClusterEvsN[iHist]->Integral()));
          hEn[iHist]->GetXaxis()->SetBinLabel(iFile + 1, Form("%d", run));

          TH1D *hClusterE = (TH1D*)hClusterEvsN[iHist]->ProjectionX("hClusterE", 1, hClusterEvsN[iHist]->GetNbinsY());
	  TH1D *hex = (TH1D*)hClusterE->Rebin(1, "hex", xa);
	  hex->Scale(1./nev);
          hNClust[iHist]->SetBinContent(iFile + 1, hex->GetBinContent(1));
	  hNClust[iHist]->SetBinError(iFile + 1, hex->GetBinError(1));
	  hNClust[iHist]->GetXaxis()->SetBinLabel(iFile + 1, Form("%d", run));

	  hNClust2[iHist]->SetBinContent(iFile + 1, hpx->GetBinContent(1));
	  hNClust2[iHist]->SetBinError(iFile + 1, hpx->GetBinError(1));
	  hNClust2[iHist]->GetXaxis()->SetBinLabel(iFile + 1, Form("%d", run));


          TH1D *hClusterN = (TH1D*)hClusterEvsN[iHist]->ProjectionY("hClusterN", 1, hClusterEvsN[iHist]->GetNbinsX());
	  hClusterN->GetXaxis()->SetRangeUser(3, 40);
	  Double_t ncpcl = hClusterN->GetMean(); 
	  hNCell[iHist]->SetBinContent(iFile + 1, ncpcl);
	  hNCell[iHist]->SetBinError(iFile + 1, ncpcl/TMath::Sqrt(hClusterN->GetEntries()));
          hNCell[iHist]->GetXaxis()->SetBinLabel(iFile + 1, Form("%d", run));
      }

      iFile ++;

      delete f1;
      delete Data;
   } 
   
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0); 
   gStyle->SetPaperSize(10.,10.);

   gSystem->Exec(Form("mkdir -p figures/Cluster_QA_%s/", period.Data()));

   hEn[0] ->SetLineColor(kBlack);
   hEn[0] ->SetMarkerColor(kBlack);
   hEn[0] ->SetMarkerSize(0.8);
   hEn[0] ->SetMarkerStyle(kFullCircle);

   hEn[1]->SetLineColor(kBlue);
   hEn[1]->SetMarkerColor(hEn[1]->GetLineColor());
   hEn[1]->SetMarkerSize(0.8);
   hEn[1]->SetMarkerStyle(kFullCircle);

   hEn[2]->SetLineColor(kGreen);
   hEn[2]->SetMarkerColor(hEn[2]->GetLineColor());
   hEn[2]->SetMarkerSize(0.8);
   hEn[2]->SetMarkerStyle(kFullCircle);

   hEn[3]->SetLineColor(kOrange);
   hEn[3]->SetMarkerColor(hEn[3]->GetLineColor());
   hEn[3]->SetMarkerSize(0.8);
   hEn[3]->SetMarkerStyle(kFullCircle);

   hEn[4]->SetLineColor(kMagenta);
   hEn[4]->SetMarkerColor(hEn[4]->GetLineColor());
   hEn[4]->SetMarkerSize(0.8);
   hEn[4]->SetMarkerStyle(kFullCircle);

   TCanvas *cen = new TCanvas("cen", "cen", 1200, 600);

   cen->SetGridx();
   cen->SetGridy();

   for (Int_t iHist = 0; iHist < 0; iHist ++) {
     hEn[iHist]->Draw(Form("%s", iHist == 0 ? "" : "same"));
   }

   TPaveText *t1 = new TPaveText(0.2, 0.92, 0.8, 1.0, "brNDC"); // left-up   
   t1->SetBorderSize(0);
   t1->SetFillColor(kWhite);
   t1->AddText(Form("Average energy of clusters, %s", period.Data()));
   t1->Draw();
  
   cen->BuildLegend(0.82, 0.2, 0.9, 0.45);
  
   cen->Print(Form("figures/Cluster_QA_%s/cen_%s.eps", period.Data(), cut.Data())); 
   cen->Print(Form("figures/Cluster_QA_%s/cen_%s.tex", period.Data(), cut.Data())); 
   cen->Print(Form("figures/Cluster_QA_%s/cen_%s.pdf", period.Data(), cut.Data())); 

   hNClust[0]->SetLineColor(kBlack);
   hNClust[0]->SetMarkerColor(hNClust[0]->GetLineColor());
   hNClust[0]->SetMarkerSize(0.8);
   hNClust[0]->SetMarkerStyle(kFullCircle);

   hNClust[1]->SetLineColor(kBlue);
   hNClust[1]->SetMarkerColor(hNClust[1]->GetLineColor());
   hNClust[1]->SetMarkerSize(0.8);
   hNClust[1]->SetMarkerStyle(kFullCircle);

   hNClust[2]->SetLineColor(kGreen);
   hNClust[2]->SetMarkerColor(hNClust[2]->GetLineColor());
   hNClust[2]->SetMarkerSize(0.8);
   hNClust[2]->SetMarkerStyle(kFullCircle);

   hNClust[3]->SetLineColor(kOrange);
   hNClust[3]->SetMarkerColor(hNClust[3]->GetLineColor());
   hNClust[3]->SetMarkerSize(0.8);
   hNClust[3]->SetMarkerStyle(kFullCircle);

   hNClust[4]->SetLineColor(kMagenta);
   hNClust[4]->SetMarkerColor(hNClust[4]->GetLineColor());
   hNClust[4]->SetMarkerSize(0.8);
   hNClust[4]->SetMarkerStyle(kFullCircle);

   TCanvas *cN = new TCanvas("cN", "cN", 1200, 600);

   cN->SetGridx();
   cN->SetGridy();

   for (Int_t iHist = 0; iHist < 5; iHist ++) {
     hNClust[iHist]->GetYaxis()->SetRangeUser(0, 0.07);
     hNClust[iHist]->Draw(Form("%s", iHist == 0 ? "" : "same"));
   }

   cN->BuildLegend(0.82, 0.4, 0.9, 0.65);

   TPaveText *t2 = new TPaveText(0.2, 0.92, 0.8, 1.0, "brNDC"); // left-up   
   t2->SetBorderSize(0);
   t2->SetFillColor(kWhite);
   t2->AddText(Form("Average number of clusters per event, %s, E > %0.1f, N_{cells} #geq %d", period.Data(), emin, nCellMin));
   t2->Draw();
 
   cN->Print(Form("figures/Cluster_QA_%s/cN_%s.eps", period.Data(), cut.Data()));   
   cN->Print(Form("figures/Cluster_QA_%s/cN_%s.tex", period.Data(), cut.Data()));    
   cN->Print(Form("figures/Cluster_QA_%s/cN_%s.pdf", period.Data(), cut.Data()));        

   hNClust2[0]->SetLineColor(kBlack);
   hNClust2[0]->SetMarkerColor(hNClust[0]->GetLineColor());
   hNClust2[0]->SetMarkerSize(0.8);
   hNClust2[0]->SetMarkerStyle(kFullCircle);

   hNClust2[1]->SetLineColor(kBlue);
   hNClust2[1]->SetMarkerColor(hNClust[1]->GetLineColor());
   hNClust2[1]->SetMarkerSize(0.8);
   hNClust2[1]->SetMarkerStyle(kFullCircle);

   hNClust2[2]->SetLineColor(kGreen);
   hNClust2[2]->SetMarkerColor(hNClust[2]->GetLineColor());
   hNClust2[2]->SetMarkerSize(0.8);
   hNClust2[2]->SetMarkerStyle(kFullCircle);

   hNClust2[3]->SetLineColor(kOrange);
   hNClust2[3]->SetMarkerColor(hNClust[3]->GetLineColor());
   hNClust2[3]->SetMarkerSize(0.8);
   hNClust2[3]->SetMarkerStyle(kFullCircle);

   hNClust2[4]->SetLineColor(kMagenta);
   hNClust2[4]->SetMarkerColor(hNClust[4]->GetLineColor());
   hNClust2[4]->SetMarkerSize(0.8);
   hNClust2[4]->SetMarkerStyle(kFullCircle);

   TCanvas *cN2 = new TCanvas("cN2", "cN2", 1200, 600);

   cN2->SetGridx();
   cN2->SetGridy();

   for (Int_t iHist = 0; iHist < 5; iHist ++) {
     hNClust2[iHist]->Draw(Form("%s", iHist == 0 ? "" : "same"));
   }

   cN2->BuildLegend(0.82, 0.4, 0.9, 0.65);

   TPaveText *t22 = new TPaveText(0.2, 0.92, 0.8, 1.0, "brNDC"); // left-up   
   t22->SetBorderSize(0);
   t22->SetFillColor(kWhite);
   t22->AddText(Form("Average number of clusters per event, %s, E > %0.1f, N_{cells} #geq %d", period.Data(), emin, nCellMin));
   t22->Draw();
 
   cN2->Print(Form("figures/Cluster_QA_%s/cN2_%s.eps", period.Data(), cut.Data()));   
   cN2->Print(Form("figures/Cluster_QA_%s/cN2_%s.tex", period.Data(), cut.Data()));    
   cN2->Print(Form("figures/Cluster_QA_%s/cN2_%s.pdf", period.Data(), cut.Data()));        
 
 
   hNCell[0] ->SetLineColor(kBlack);
   hNCell[0] ->SetMarkerColor(hNCell[0]->GetLineColor());
   hNCell[0] ->SetMarkerSize(0.8);
   hNCell[0] ->SetMarkerStyle(kFullCircle);

   hNCell[1]->SetLineColor(kBlue);
   hNCell[1]->SetMarkerColor(hNCell[1]->GetLineColor());
   hNCell[1]->SetMarkerSize(0.8);
   hNCell[1]->SetMarkerStyle(kFullCircle);

   hNCell[2]->SetLineColor(kGreen);
   hNCell[2]->SetMarkerColor(hNCell[2]->GetLineColor());
   hNCell[2]->SetMarkerSize(0.8);
   hNCell[2]->SetMarkerStyle(kFullCircle);

   hNCell[3]->SetLineColor(kOrange);
   hNCell[3]->SetMarkerColor(hNCell[3]->GetLineColor());
   hNCell[3]->SetMarkerSize(0.8);
   hNCell[3]->SetMarkerStyle(kFullCircle);

   hNCell[4]->SetLineColor(kMagenta);
   hNCell[4]->SetMarkerColor(hNCell[4]->GetLineColor());
   hNCell[4]->SetMarkerSize(0.8);
   hNCell[4]->SetMarkerStyle(kFullCircle);

   TCanvas *cNcell = new TCanvas("cNcell", "cNcell", 1200, 600);

   cNcell->SetGridx();
   cNcell->SetGridy();

   for (Int_t iHist = 0; iHist < 5; iHist ++ ) {
	      hNCell[iHist]->Draw(Form("%s", iHist == 0 ? "" : "same"));
   }

   TPaveText *t3 = new TPaveText(0.2, 0.92, 0.8, 1.0, "brNDC"); // left-up   
   t3->SetBorderSize(0);
   t3->SetFillColor(kWhite);
   t3->AddText(Form("Average number of cells in cluster, %s, E > %0.1f, N_{cells} > %d", period.Data(), emin, nCellMin));
   t3->Draw();
  
   cNcell->BuildLegend(0.82, 0.4, 0.9, 0.65);   

   cNcell->Print(Form("figures/Cluster_QA_%s/cNcell_%s.eps", period.Data(), cut.Data()));
   cNcell->Print(Form("figures/Cluster_QA_%s/cNcell_%s.tex", period.Data(), cut.Data()));   
   cNcell->Print(Form("figures/Cluster_QA_%s/cNcell_%s.pdf", period.Data(), cut.Data()));

  TFile fout(Form("ROOT/QA/Cluster_QA_%s_%s.root", period.Data(), cut.Data()), "recreate");
 
 hEvents->Write();

  for (Int_t  iHist = 0; iHist < 5; iHist++) {
     hNClust[iHist]->Write();
     hNClust2[iHist]->Write();
     hNCell[iHist]->Write();
     hEn[iHist]->Write();
  }

  cN    ->Write();
  cN2   ->Write();
  cen   ->Write();
  cNcell->Write();

  fout.Close();

  printf ("\n output stored to the file %s \n", fout.GetName());
}
