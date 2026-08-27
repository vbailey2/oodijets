int GlobalBin(int ix, int iy, int nBinsY)
{
    return (ix - 1) * nBinsY + iy;
}
 
// -----------------------------------------------------------------
// inverse of GlobalBin
// -----------------------------------------------------------------
void InverseGlobalBin(int globalBin, int nBinsY, int &ix, int &iy)
{
    int g = globalBin - 1; // 0-based
    ix = g / nBinsY + 1;
    iy = g % nBinsY + 1;
}


void unfoldData(int nIterations = 2)
{
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  
  TFile *fMC = new TFile("hists/histMC.root","READ");
  TFile *fData = new TFile("hists/hist-full.root","READ");
  TH3F *h_xj_bins = (TH3F*)fMC->Get("h_xj");
  TH3F *h_pt1pt2 = (TH3F*)fMC->Get("h_pt1pt2");
  int cent_N = h_pt1pt2->GetNbinsZ();
  int pt_N = h_pt1pt2->GetNbinsX();
  
  TH2D *hTrue2D[cent_N];
  TH2D *hMeas2D[cent_N];
  TH2D *hUnfolded2D[cent_N];
  TH1D *hTrue1D[cent_N];
  TH1D *hMeas1D[cent_N];
  TH1D *hFake1D[cent_N];
  TH1D *hUnfolded1D[cent_N];

  
  TH3F *hUnfolded3D = (TH3F*)h_pt1pt2->Clone();
  hUnfolded3D->SetName("h_unfoldpt1pt2");
  hUnfolded3D->Reset();
  TH3F *hTrue3D = (TH3F*)h_pt1pt2->Clone();
  hTrue3D->SetName("h_truept1pt2");
  hTrue3D->Reset();
  TH3F *hMeas3D = (TH3F*)h_pt1pt2->Clone();
  hMeas3D->SetName("h_measpt1pt2");
  hMeas3D->Reset();
   
 
  std::vector<RooUnfoldResponse*> response(cent_N);
  std::vector<RooUnfoldBayes>     unfold;
  unfold.reserve(cent_N);

  
  for(int ic = 0; ic < cent_N; ic++)
    {
      hTrue2D[ic] = (TH2D*)fMC->Get(Form("hTrue2D%i",ic));
      hMeas2D[ic] = (TH2D*)fData->Get(Form("hMeas2D%i",ic));
      hTrue1D[ic] = (TH1D*)fMC->Get(Form("hTrue1D%i",ic));
      hMeas1D[ic] = (TH1D*)fData->Get(Form("hMeas1D%i",ic));
      hFake1D[ic] = (TH1D*)fMC->Get(Form("hFake1D%i",ic));
      response[ic] = (RooUnfoldResponse*)fMC->Get(Form("response%i",ic));

      
      const int nBinsY = hTrue2D[ic]->GetNbinsY();
      const int nGlobalBins = hMeas1D[ic]->GetNbinsX();

      //run unfold
      hMeas1D[ic]->Add(hFake1D[ic],-1);
      unfold.emplace_back(response[ic], hMeas1D[ic], nIterations);
      unfold[ic].SetVerbose(1);
      
      hUnfolded1D[ic] = (TH1D*) unfold[ic].Hunfold(RooUnfold::kErrors);
      hUnfolded1D[ic]->SetName(Form("hUnfolded1D_cent%d", ic));
      
      //map back to 2D
      hUnfolded2D[ic] = (TH2D*) hTrue2D[ic]->Clone(Form("hUnfolded2D_cent%d", ic));
      hUnfolded2D[ic]->SetTitle(Form("Unfolded (Bayes), cent %i", ic));
      hUnfolded2D[ic]->SetDirectory(nullptr);
      hUnfolded2D[ic]->Reset();

      //also map 1D truth/reco to 2D
      hTrue2D[ic]->Reset();
      hMeas2D[ic]->Reset();

      
      for (int g = 1; g <= nGlobalBins; ++g) {
	int ix, iy;
	InverseGlobalBin(g, nBinsY, ix, iy);
	hUnfolded2D[ic]->SetBinContent(ix, iy, hUnfolded1D[ic]->GetBinContent(g));
	hUnfolded2D[ic]->SetBinError(ix, iy, hUnfolded1D[ic]->GetBinError(g));
	
	hTrue2D[ic]->SetBinContent(ix, iy, hTrue1D[ic]->GetBinContent(g));
	hTrue2D[ic]->SetBinError(ix, iy, hTrue1D[ic]->GetBinError(g));
	
	hMeas2D[ic]->SetBinContent(ix, iy, hMeas1D[ic]->GetBinContent(g));
        hMeas2D[ic]->SetBinError(ix, iy, hMeas1D[ic]->GetBinError(g));
      }      
    }

  TFile *fout = new TFile("hists/hist-unfoldedData.root","RECREATE");
  h_xj_bins->Write();
  h_pt1pt2->Write();
  for(int ic = 0; ic < cent_N; ic++)
    {
      hTrue2D[ic]->Write();
      hMeas2D[ic]->Write();
      hUnfolded2D[ic]->Write();

      hUnfolded1D[ic]->Write();
      hTrue1D[ic]->Write();
      hMeas1D[ic]->Write();
    }
}
