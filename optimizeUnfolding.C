#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

void optimizeUnfolding()
{
  const int nIterMin = 1;
  const int nIterMax = 8;
  const int nIter = nIterMax - nIterMin + 1;

  TFile *f[nIter];
  for (int i = 0; i < nIter; i++)
  {
    int iter = nIterMin + i;
    f[i] = new TFile(Form("hists/hist-unfoldedData_iter%d.root", iter), "READ");
  }

  TH3F *h_pt1pt2 = (TH3F *)f[0]->Get("h_pt1pt2");
  int ncent = h_pt1pt2->GetNbinsZ();

  // 1. sum of per-bin statistical error, and 2. quadrature sum of the bin-by-bin
  // difference from the previous iteration, both per iteration and centrality
  double statErrSum[nIter][ncent];
  double diffQuadSum[nIter][ncent];

  TH2F *hUnfolded2D[nIter][ncent];
  for (int i = 0; i < nIter; i++)
  {
    for (int ic = 0; ic < ncent; ic++)
    {
      hUnfolded2D[i][ic] = (TH2F *)f[i]->Get(Form("hUnfolded2D_cent%d", ic));

      int nx = hUnfolded2D[i][ic]->GetNbinsX();
      int ny = hUnfolded2D[i][ic]->GetNbinsY();

      double errsum = 0;
      for (int ix = 1; ix <= nx; ix++)
      {
        for (int iy = 1; iy <= ny; iy++)
        {
          errsum += hUnfolded2D[i][ic]->GetBinError(ix, iy);
        }
      }
      statErrSum[i][ic] = errsum;

      diffQuadSum[i][ic] = 0;
      if (i > 0)
      {
        double diffsum = 0;
        for (int ix = 11; ix <= nx; ix++) //start at bin 11 since thats our first pt bin
        {
          for (int iy = 11; iy <= ny; iy++)
          {
            double diff = hUnfolded2D[i][ic]->GetBinContent(ix, iy) - hUnfolded2D[i - 1][ic]->GetBinContent(ix, iy);
            diffsum += diff * diff;
          }
        }
        diffQuadSum[i][ic] = std::sqrt(diffsum);
      }
    }
  }

  // 3. plot both quantities, and their quadrature sum, vs iteration number for each centrality
  SetsPhenixStyle();



  int colors[] = {1, 2, 4, kGreen + 2, kViolet, kCyan, kOrange + 2, kMagenta + 2, kAzure - 2};

  for (int ic = 0; ic < ncent; ic++)
  {
    TGraph *gStatErr = new TGraph(nIter);
    gStatErr->SetName(Form("gStatErrSum_cent%d", ic));
    gStatErr->SetTitle(Form("cent %d;iteration;statistical error", ic));

    TGraph *gDiff = new TGraph(nIter - 1);
    gDiff->SetName(Form("gDiffQuadSum_cent%d", ic));
    gDiff->SetTitle(Form("cent %d;iteration;difference from previous iteration", ic));

    TGraph *gQuad = new TGraph(nIter - 1);
    gQuad->SetName(Form("gQuadSum_cent%d", ic));
    gQuad->SetTitle(Form("cent %d;iteration;quadrature sum", ic));

    for (int i = 0; i < nIter; i++)
    {
      int iter = nIterMin + i;
      gStatErr->SetPoint(i, iter, statErrSum[i][ic]);

      if (i > 0)
      {
        double quad = std::sqrt(statErrSum[i][ic] * statErrSum[i][ic] + diffQuadSum[i][ic] * diffQuadSum[i][ic]);
        gDiff->SetPoint(i - 1, iter, diffQuadSum[i][ic]);
        gQuad->SetPoint(i - 1, iter, quad);
      }
    }

    gStatErr->SetMarkerStyle(20);
    gStatErr->SetMarkerColor(colors[0]);
    gStatErr->SetLineColor(colors[0]);

    gDiff->SetMarkerStyle(21);
    gDiff->SetMarkerColor(colors[1]);
    gDiff->SetLineColor(colors[1]);

    gQuad->SetMarkerStyle(22);
    gQuad->SetMarkerColor(colors[2]);
    gQuad->SetLineColor(colors[2]);

    TCanvas *c = new TCanvas(Form("c_cent%d", ic), Form("c_cent%d", ic), 700, 700);

    TMultiGraph *mg = new TMultiGraph();
    mg->SetName(Form("mg_cent%d", ic));
    mg->SetTitle(Form("cent %d;iteration;", ic));
    mg->Add(gStatErr, "LP");
    mg->Add(gDiff, "LP");
    mg->Add(gQuad, "LP");
    mg->Draw("A");

    TLegend *leg = new TLegend(.15, .75, .4, .92);
	  leg->SetFillStyle(0);
	  leg->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
	  leg->AddEntry("", "O+O #sqrt{s_{NN}} = 200 GeV", "");
	  leg->AddEntry("", "anti-#it{k}_{#it{t}} #it{R} = 0.4, |#eta| < 0.7", "");

    TLegend *leg2 = new TLegend(.75, .75, .88, .92);
    leg2->SetFillStyle(0);
    leg2->AddEntry(gStatErr, "#sigma_{stat}", "lp");
    leg2->AddEntry(gDiff, "#sigma_{bin}", "lp");
    leg2->AddEntry(gQuad, "Total", "lp");
    leg->Draw();
    leg2->Draw();

    c->Print(Form("plots/optimizeUnfolding_cent%d.pdf", ic));
  }

  std::cout << "all done" << std::endl;
}
