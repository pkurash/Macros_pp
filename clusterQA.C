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
	        Int_t nCellMin = 1)
{
   std::vector<Int_t> runNumbers, zeroNumbers;
   Int_t nn;

   ifstream ff;
   ofstream ffzero;

    //printf("%s.txt\n", period.Data());
    ff.open(Form("./datasets/%s.txt", period.Data()));
    while ( ff >> nn) { 
      // cout << nn << endl;
      runNumbers.emplace_back(nn);
    }
    ff.close();


    const Int_t nRuns = runNumbers.size();
    cout << nRuns << " runs:" << endl;

    TFile *f1;
    TList *hData;

    TH2F *hClusterEvsN[5];
    TH1F *hEn[5], *hNClust[5], *hNCell[5], *hNClust2[5], *hClustPt[5], *hcounter;
    Double_t xa[2]={emin, 20.};

    TH1F *hEvents = new TH1F("hEvents", "number of analyzed events;;N_{ev}", nRuns, 0, 1.0*nRuns);

    for (Int_t iHist = 0; iHist < 5; iHist ++ ) {
      TString sMod = iHist == 0 ? "all modules" : Form("module %d", iHist);
      hEn[iHist] =      new TH1F(Form("hEn%d",     iHist), Form("%s;;<E_{cl}>, GeV/c", sMod.Data()), nRuns, 0, nRuns); 
      hNClust[iHist] =  new TH1F(Form("hNClust%d", iHist), Form("%s;;N_{cl}",          sMod.Data()), nRuns, 0, nRuns); 
      hNClust2[iHist] = new TH1F(Form("hNClust2%d",iHist), Form("%s;;N_{cl}",          sMod.Data()), nRuns, 0, nRuns); 
      hNCell[iHist] =   new TH1F(Form("hNCell%d",  iHist), Form("%s;;N_{cell}",        sMod.Data()), nRuns, 0, nRuns); 
    }

    Int_t iFile = 0;
 
    for (auto run : runNumbers) {
      
      //printf("run: %d\n", run);

      TString inFilePath = Form("../MergedResults_Run2/%s/%d/AnalysisResults.root", period.Data(), run);
      
      f1 = TFile::Open(inFilePath);    
      if (!f1) { 
        printf("ERROR: no file for run %d found!\n", run);
        return;
      }

      hData = (TList*)f1->Get("Data");
 
      hcounter =       (TH1F*)hData->FindObject("hSelEvents");
      Int_t nev =hcounter->GetBinContent(8);
      if (nev == 0) {
         printf("ERROR: no events in the file for run %d!\n", run);
	 return;
      }

      hEvents->SetBinContent(iFile + 1, nev);
      hEvents->GetXaxis()->SetBinLabel(iFile + 1, Form("%d", run));
      Double_t en[5], ncl[5], ncl2[5], ncpcl[5];

      for (Int_t i = 0; i < 5; i ++) {
          hClusterEvsN[i]   = (TH2F*)hData->FindObject(Form("hClusterEvsN_%s%s", 
	                                       cut.Data(), i == 0 ? "" : Form("_M%d", i)));
          hClusterEvsN[i]->GetXaxis()->SetRangeUser(emin, 0.01*hClusterEvsN[i]->GetNbinsX());
          hClusterEvsN[i]->GetYaxis()->SetRange(nCellMin, hClusterEvsN[i]->GetNbinsY());
          
	  hClustPt[i] = (TH1F*)hData->FindObject(Form("hCaloPhotonPt_%s%s", cut.Data(),
	                                       i == 0 ? "" : Form("_M%d", i)));
          TH1F *hpx = (TH1F*)hClustPt[i]->Rebin(1, "h", xa);
	  hpx->Scale(1./nev);

          en[i] = hClusterEvsN[i]->GetMean(1);
	  hEn[i]->SetBinContent(iFile + 1, en[i]);
          hEn[i]->SetBinError(iFile + 1, 
	                  hClusterEvsN[i]->GetMean(1)/TMath::Sqrt(hClusterEvsN[i]->Integral()));
          hEn[i]->GetXaxis()->SetBinLabel(iFile + 1, Form("%d", run));

          TH1D *hClusterE = (TH1D*)hClusterEvsN[i]->ProjectionX("hClusterE", 1, hClusterEvsN[i]->GetNbinsY());
	  TH1D *hex = (TH1D*)hClusterE->Rebin(1, "hex", xa);

	  hex->Scale(1./nev);
	  ncl[i] = hex->GetBinContent(1);
          hNClust[i]->SetBinContent(iFile + 1, ncl[i]);
	  hNClust[i]->SetBinError(iFile + 1, hex->GetBinError(1));
	  hNClust[i]->GetXaxis()->SetBinLabel(iFile + 1, Form("%d", run));

	  ncl2[i] = hpx->GetBinContent(1);
	  hNClust2[i]->SetBinContent(iFile + 1, ncl2[i]);
	  hNClust2[i]->SetBinError(iFile + 1, hpx->GetBinError(1));
	  hNClust2[i]->GetXaxis()->SetBinLabel(iFile + 1, Form("%d", run));

          TH1D *hClusterN = (TH1D*)hClusterEvsN[i]->ProjectionY("hClusterN", 1, hClusterEvsN[i]->GetNbinsX());
	  hClusterN->GetXaxis()->SetRangeUser(3, 40);
	  ncpcl[i] = hClusterN->GetMean(); 
	  hNCell[i]->SetBinContent(iFile + 1, ncpcl[i]);
	  hNCell[i]->SetBinError(iFile + 1, ncpcl[i]/TMath::Sqrt(hClusterN->GetEntries()));
          hNCell[i]->GetXaxis()->SetBinLabel(iFile + 1, Form("%d", run));

      }

      if (en[0] == 0 ||  ncl[0] == 0 || ncpcl[0] == 0) {
          printf("run %d: ZERO ALERT!!!\n", run);
	  zeroNumbers.emplace_back(run);
       }
        
      iFile ++;

      delete f1;
      delete hData;
   } 

   ffzero.open(Form("./datasets/%s_bad.txt", period.Data()));
   for (auto zero : zeroNumbers) {
    ffzero << zero << endl;
   }
   
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0); 
   gStyle->SetPaperSize(10.,10.);

   TString figDir = Form("figures/Cluster_QA_%s", period.Data());
   gSystem->Exec(Form("mkdir -p %s", figDir.Data()));
   gSystem->Exec("mkdir -p ROOT/QA/");

   EColor hColor[5] = {kBlack, kBlue, kGreen, kMagenta, kOrange};

   for (Int_t i = 0; i < 5; i ++) {
     hEn[i]->SetLineColor(hColor[i]);
     hEn[i]->SetMarkerColor(hEn[i]->GetLineColor());
     hEn[i]->SetMarkerSize(0.8);
     hEn[i]->SetMarkerStyle(kFullCircle);
     hEn[i]->GetYaxis()->SetRangeUser(0.2,1.02);
     hEn[i]->Sumw2(0);
   }
   
   auto *cen = new TCanvas("cen", "cen", 1200, 600);
   //cen->SetGridx();
   //cen->SetGridy();
   cen->cd();

   for (Int_t iHist = 0; iHist < 5; iHist ++) {
     hEn[iHist]->Draw(Form("%s", iHist == 0 ? "" : "same"));
   }

   auto *t1 = new TPaveText(0.2, 0.92, 0.8, 1.0, "brNDC"); // left-up   
   t1->SetBorderSize(0);
   t1->SetFillColor(kWhite);
   t1->AddText(Form("Average energy of clusters, %s", period.Data()));
   t1->Draw();
  
   cen->BuildLegend(0.82, 0.2, 0.9, 0.45);
 
   TString cenPath = Form("%s/cen_%s", figDir.Data(), cut.Data());
   cen->SaveAs(Form("%s.root", cenPath.Data())); 
   cen->SaveAs(Form("%s.eps",  cenPath.Data())); 
   cen->SaveAs(Form("%s.pdf",  cenPath.Data())); 

   for (Int_t i = 0; i < 5; i ++) {
     hNClust[i]->SetLineColor(hColor[i]);
     hNClust[i]->SetMarkerColor(hNClust[i]->GetLineColor());
     hNClust[i]->SetMarkerSize(0.8);
     hNClust[i]->SetMarkerStyle(kFullCircle);
   }

   auto *cN = new TCanvas("cN", "cN", 1200, 600);
   //cN->SetGridx();
   //cN->SetGridy();
   cN->cd();

   for (Int_t iHist = 0; iHist < 5; iHist ++) {
     hNClust[iHist]->GetYaxis()->SetRangeUser(0, 0.07);
     hNClust[iHist]->Draw(Form("%s", iHist == 0 ? "" : "same"));
   }

   cN->BuildLegend(0.82, 0.4, 0.9, 0.65);

   auto *t2 = new TPaveText(0.2, 0.92, 0.8, 1.0, "brNDC"); // left-up   
   t2->SetBorderSize(0);
   t2->SetFillColor(kWhite);
   t2->AddText(Form("Average number of clusters per event, %s, E > %0.1f, N_{cells} #geq %d", period.Data(), emin, nCellMin));
   t2->Draw();
 
   TString cNPath = Form("%s/cN_%s", figDir.Data(), cut.Data());
   cN->SaveAs(Form("%s.root", cNPath.Data()));   
   cN->SaveAs(Form("%s.eps", cNPath.Data()));   
   cN->SaveAs(Form("%s.pdf", cNPath.Data()));   

   for (Int_t i = 0; i < 5; i ++) {
     hNClust2[i]->SetLineColor(hColor[i]);
     hNClust2[i]->SetMarkerColor(hNClust2[i]->GetLineColor());
     hNClust2[i]->SetMarkerSize(0.8);
     hNClust2[i]->SetMarkerStyle(kFullCircle);
   }

   auto *cN2 = new TCanvas("cN2", "cN2", 1200, 600);
   //cN2->SetGridx();
   //cN2->SetGridy();
   cN2->cd();
 
   for (Int_t iHist = 0; iHist < 5; iHist ++) {
     hNClust2[iHist]->Draw(Form("%s", iHist == 0 ? "" : "same"));
   }

   cN2->BuildLegend(0.82, 0.4, 0.9, 0.65);

   auto *t22 = new TPaveText(0.2, 0.92, 0.8, 1.0, "brNDC"); // left-up   
   t22->SetBorderSize(0);
   t22->SetFillColor(kWhite);
   t22->AddText(Form("Average number of clusters per event, %s, E > %0.1f, N_{cells} #geq %d", period.Data(), emin, nCellMin));
   t22->Draw();
 
   TString cNPath2 = Form("%s/cN2_%s", figDir.Data(), cut.Data());
   cN2->SaveAs(Form("%s.root", cNPath2.Data()));   
   cN2->SaveAs(Form("%s.eps",  cNPath2.Data()));    
   cN2->SaveAs(Form("%s.pdf",  cNPath2.Data()));        
 
   for (Int_t i = 0; i < 5; i ++) {
     hNCell[i]->SetLineColor(hColor[i]);
     hNCell[i]->SetMarkerColor(hNCell[i]->GetLineColor());
     hNCell[i]->SetMarkerSize(0.8);
     hNCell[i]->SetMarkerStyle(kFullCircle);
   }
 
   auto *cNcell = new TCanvas("cNcell", "cNcell", 1200, 600);
   //cNcell->SetGridx();
   //cNcell->SetGridy();
   cNcell->cd();

   for (Int_t iHist = 0; iHist < 5; iHist ++ ) {
     hNCell[iHist]->Draw(Form("%s", iHist == 0 ? "" : "same"));
   }

   auto *t3 = new TPaveText(0.2, 0.92, 0.8, 1.0, "brNDC"); // left-up   
   t3->SetBorderSize(0);
   t3->SetFillColor(kWhite);
   t3->AddText(Form("Average number of cells in cluster, %s, E > %0.1f, N_{cells} > %d", period.Data(), emin, nCellMin));
   t3->Draw();
  
   cNcell->BuildLegend(0.82, 0.4, 0.9, 0.65);   
   
   TString cNcellPath = Form("%s/cNcell_%s.eps", figDir.Data(), cut.Data());
   cNcell->SaveAs(Form("%s.root", cNcellPath.Data()));
   cNcell->SaveAs(Form("%s.eps",  cNcellPath.Data()));   
   cNcell->SaveAs(Form("%s.pdf",  cNcellPath.Data()));

   auto *fout = TFile::Open(Form("ROOT/QA/Cluster_QA_%s_%s.root", period.Data(), cut.Data()), "recreate");
 
   hEvents->Write();

   for (Int_t  iHist = 0; iHist < 5; iHist++) {
      hNClust[iHist]  ->Write();
      hNClust2[iHist] ->Write();
      hNCell[iHist]   ->Write();
      hEn[iHist]      ->Write();
   }

   cN    ->Write();
   cN2   ->Write();
   cen   ->Write();
   cNcell->Write();

   printf ("\n storing output to file %s \n", fout->GetName());
   fout->Close();
}
