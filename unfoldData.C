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

// nToys: number of response-matrix toy unfolds for the response-matrix statistical uncertainty
// (added in quadrature to the data statistical uncertainty); nToys <= 1 disables it
void unfoldData(int nIterations = 2, int nToys = 1000)
{
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();

  TFile *fMC = new TFile("hists/histMC_reweight.root", "READ");
  TFile *fData = new TFile("hists/hist-full.root", "READ");
  TH3F *h_xj_bins = (TH3F *)fMC->Get("h_xj");
  TH3F *h_pt1pt2 = (TH3F *)fMC->Get("h_pt1pt2");
  int cent_N = h_pt1pt2->GetNbinsZ();
  int pt_N = h_pt1pt2->GetNbinsX();

  TH2D *hTrue2D[cent_N];
  TH2D *hMeas2D[cent_N];
  TH2D *hUnfolded2D[cent_N];
  TH1D *hTrue1D[cent_N];
  TH1D *hMeas1D[cent_N];
  TH1D *hUnfolded1D[cent_N];
  TH1D *hUnfolded1D_dataStat[cent_N]; // unfolded with data-only stat errors
  TH1D *hUnfolded1D_respStat[cent_N]; // content = response-matrix stat error per bin

  TH3F *hUnfolded3D = (TH3F *)h_pt1pt2->Clone();
  hUnfolded3D->SetName("h_unfoldpt1pt2");
  hUnfolded3D->Reset();
  TH3F *hTrue3D = (TH3F *)h_pt1pt2->Clone();
  hTrue3D->SetName("h_truept1pt2");
  hTrue3D->Reset();
  TH3F *hMeas3D = (TH3F *)h_pt1pt2->Clone();
  hMeas3D->SetName("h_measpt1pt2");
  hMeas3D->Reset();

  std::vector<RooUnfoldResponse *> response(cent_N);
  std::vector<RooUnfoldBayes> unfold;
  unfold.reserve(cent_N);

  for (int ic = 0; ic < cent_N; ic++)
  {
    hTrue2D[ic] = (TH2D *)fMC->Get(Form("hTrue2D%i", ic));
    hMeas2D[ic] = (TH2D *)fData->Get(Form("hMeas2D%i", ic));
    hTrue1D[ic] = (TH1D *)fMC->Get(Form("hTrue1D%i", ic));
    hMeas1D[ic] = (TH1D *)fData->Get(Form("hMeas1D%i", ic));
    response[ic] = (RooUnfoldResponse *)fMC->Get(Form("response%i", ic));

    const int nBinsY = hTrue2D[ic]->GetNbinsY();
    const int nGlobalBins = hMeas1D[ic]->GetNbinsX();

    // run unfold; fakes are not subtracted here -- the response matrix carries them in an extra
    // truth bin (index nGlobalBins+1), which the g <= nGlobalBins loop below drops
    unfold.emplace_back(response[ic], hMeas1D[ic], nIterations);
    unfold[ic].SetVerbose(1);

    hUnfolded1D[ic] = (TH1D *)unfold[ic].Hunfold(RooUnfold::kErrors);
    hUnfolded1D[ic]->SetName(Form("hUnfolded1D_cent%d", ic));

    // keep the data-only statistical error (RooUnfold's propagated measured-data error) for reference
    hUnfolded1D_dataStat[ic] = (TH1D *)hUnfolded1D[ic]->Clone(Form("hUnfolded1D_dataStat_cent%d", ic));
    hUnfolded1D_dataStat[ic]->SetDirectory(nullptr);

    // response-matrix statistical uncertainty: unfold the same data nToys times, each with a
    // response matrix whose bins are fluctuated within their MC statistical uncertainty, and take
    // the per-bin standard deviation of the unfolded results. Each bin is fluctuated as a Poisson
    // draw of its effective number of entries, n_eff = (w/sigma_w)^2, rescaled back to weighted
    // units (no negative contents, sensible for sparse bins). The truth and measured spectra
    // of the toy response are shifted by the change in the response projections, so the misses
    // (truth - projection) are held fixed while efficiency and purity fluctuate with the matrix.
    hUnfolded1D_respStat[ic] = (TH1D *)hUnfolded1D[ic]->Clone(Form("hUnfolded1D_respStat_cent%d", ic));
    hUnfolded1D_respStat[ic]->SetDirectory(nullptr);
    hUnfolded1D_respStat[ic]->Reset();
    if (nToys > 1)
    {
      TH2D *hRespNom = (TH2D *)response[ic]->Hresponse()->Clone(Form("hRespNom_cent%d", ic));
      TH1D *hMeasNom = (TH1D *)response[ic]->Hmeasured()->Clone(Form("hRespMeasNom_cent%d", ic));
      TH1D *hTrueNom = (TH1D *)response[ic]->Htruth()->Clone(Form("hRespTrueNom_cent%d", ic));
      hRespNom->SetDirectory(nullptr);
      hMeasNom->SetDirectory(nullptr);
      hTrueNom->SetDirectory(nullptr);
      TH1D *hProjXNom = hRespNom->ProjectionX(Form("hProjXNom_cent%d", ic), 0, -1);
      TH1D *hProjYNom = hRespNom->ProjectionY(Form("hProjYNom_cent%d", ic), 0, -1);
      hProjXNom->SetDirectory(nullptr);
      hProjYNom->SetDirectory(nullptr);

      const int nRespX = hRespNom->GetNbinsX();
      const int nRespY = hRespNom->GetNbinsY();
      const int nUnf = hUnfolded1D[ic]->GetNbinsX();
      std::vector<double> sum(nUnf + 2, 0.), sum2(nUnf + 2, 0.);

      // fixed seed per centrality: reproducible, and the same toy matrices are used for every
      // iteration count, so the iteration scan compares like with like
      TRandom3 rng(4357 + ic);

      for (int itoy = 0; itoy < nToys; itoy++)
      {
        TH2D *hRespToy = (TH2D *)hRespNom->Clone(Form("hRespToy_cent%d", ic));
        hRespToy->SetDirectory(nullptr);
        for (int ix = 1; ix <= nRespX; ix++)
        {
          for (int iy = 1; iy <= nRespY; iy++)
          {
            double w = hRespNom->GetBinContent(ix, iy);
            double e = hRespNom->GetBinError(ix, iy);
            if (w <= 0 || e <= 0)
              continue;
            double neff = (w / e) * (w / e);
            double toy = rng.Poisson(neff) * (w / neff);
            hRespToy->SetBinContent(ix, iy, toy);
            hRespToy->SetBinError(ix, iy, e);
          }
        }

        TH1D *hProjXToy = hRespToy->ProjectionX(Form("hProjXToy_cent%d", ic), 0, -1);
        TH1D *hProjYToy = hRespToy->ProjectionY(Form("hProjYToy_cent%d", ic), 0, -1);
        hProjXToy->SetDirectory(nullptr);
        hProjYToy->SetDirectory(nullptr);
        TH1D *hMeasToy = (TH1D *)hMeasNom->Clone(Form("hRespMeasToy_cent%d", ic));
        TH1D *hTrueToy = (TH1D *)hTrueNom->Clone(Form("hRespTrueToy_cent%d", ic));
        hMeasToy->SetDirectory(nullptr);
        hTrueToy->SetDirectory(nullptr);
        hMeasToy->Add(hProjXNom, -1);
        hMeasToy->Add(hProjXToy);
        hTrueToy->Add(hProjYNom, -1);
        hTrueToy->Add(hProjYToy);

        RooUnfoldResponse respToy(hMeasToy, hTrueToy, hRespToy, Form("respToy_cent%d", ic), "");
        RooUnfoldBayes unfoldToy(&respToy, hMeas1D[ic], nIterations);
        unfoldToy.SetVerbose(0);
        TH1D *hToy = (TH1D *)unfoldToy.Hunfold(RooUnfold::kNoError);
        for (int g = 1; g <= nUnf; g++)
        {
          double v = hToy->GetBinContent(g);
          sum[g] += v;
          sum2[g] += v * v;
        }
        delete hToy;
        delete hRespToy;
        delete hProjXToy;
        delete hProjYToy;
        delete hMeasToy;
        delete hTrueToy;
      }

      for (int g = 1; g <= nUnf; g++)
      {
        double mean = sum[g] / nToys;
        double var = (sum2[g] - nToys * mean * mean) / (nToys - 1);
        double respErr = var > 0 ? std::sqrt(var) : 0;
        hUnfolded1D_respStat[ic]->SetBinContent(g, respErr);

        // total statistical error = data stat (+) response-matrix stat, in quadrature
        double dataErr = hUnfolded1D_dataStat[ic]->GetBinError(g);
        hUnfolded1D[ic]->SetBinError(g, std::sqrt(dataErr * dataErr + respErr * respErr));
      }

      delete hRespNom;
      delete hMeasNom;
      delete hTrueNom;
      delete hProjXNom;
      delete hProjYNom;
    }

    // map back to 2D
    hUnfolded2D[ic] = (TH2D *)hTrue2D[ic]->Clone(Form("hUnfolded2D_cent%d", ic));
    hUnfolded2D[ic]->SetTitle(Form("Unfolded (Bayes), cent %i", ic));
    hUnfolded2D[ic]->SetDirectory(nullptr);
    hUnfolded2D[ic]->Reset();

    // also map 1D truth/reco to 2D
    hTrue2D[ic]->Reset();
    hMeas2D[ic]->Reset();

    for (int g = 1; g <= nGlobalBins; ++g)
    {
      int ix, iy;
      InverseGlobalBin(g, nBinsY, ix, iy);
      hUnfolded2D[ic]->SetBinContent(ix, iy, hUnfolded1D[ic]->GetBinContent(g));
      hUnfolded2D[ic]->SetBinError(ix, iy, hUnfolded1D[ic]->GetBinError(g));

      hTrue2D[ic]->SetBinContent(ix, iy, hTrue1D[ic]->GetBinContent(g));
      hTrue2D[ic]->SetBinError(ix, iy, hTrue1D[ic]->GetBinError(g));

      hMeas2D[ic]->SetBinContent(ix, iy, hMeas1D[ic]->GetBinContent(g));
      hMeas2D[ic]->SetBinError(ix, iy, hMeas1D[ic]->GetBinError(g));

      // stitch the unfolded 2D slices back into a 3D pt1-pt2-cent histogram,
      // matching the "h_pt1pt2" layout projectDijets() expects as input
      hUnfolded3D->SetBinContent(ix, iy, ic + 1, hUnfolded1D[ic]->GetBinContent(g));
      hUnfolded3D->SetBinError(ix, iy, ic + 1, hUnfolded1D[ic]->GetBinError(g));
    }
  }
  hUnfolded3D->SetName("h_pt1pt2");

  TFile *fout = new TFile(Form("hists/hist-unfoldedData_iter%d.root", nIterations), "RECREATE");
  h_xj_bins->Write();
  hUnfolded3D->Write();
  for (int ic = 0; ic < cent_N; ic++)
  {
    hTrue2D[ic]->Write();
    hMeas2D[ic]->Write();
    hUnfolded2D[ic]->Write();

    hUnfolded1D[ic]->Write();
    hUnfolded1D_dataStat[ic]->Write();
    hUnfolded1D_respStat[ic]->Write();
    hTrue1D[ic]->Write();
    hMeas1D[ic]->Write();
  }
}
