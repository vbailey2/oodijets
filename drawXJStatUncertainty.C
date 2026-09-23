#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

void drawXJStatUncertainty()
{
	SetsPhenixStyle();
	TH1::SetDefaultSumw2();

	const int finalIter = 3; // chosen from optimizeUnfolding.C's convergence result
	TFile *fproj = new TFile(Form("hists/projections_unfold_iter%d.root", finalIter), "READ");
	TFile *fpp = new TFile("hists/final_plots_pp_r04.root", "READ");

	const int ncent = 4;
	std::string cent_str[] = {"0-20%", "20-40%", "40-60%", "60-80%"};
	int colors[] = {1, 2, 4, kGreen + 2, kViolet};

	// ATLAS O+O/pp reference points (extracted from plots/oodijetatlas.pdf), overlaid on the 0-20% most-central
	// bin (our finest centrality bin, closest match to ATLAS's 0-10%), at their real central values
	TGraph *gAtlas = new TGraph("atlas_oodijet_ratio.txt");
	gAtlas->SetName("gAtlas");
	gAtlas->SetMarkerStyle(24);
	gAtlas->SetMarkerColor(kBlack);
	gAtlas->SetLineColor(kBlack);

	// per-centrality final xj histograms, plus their sum for the centrality-inclusive O+O comparison
	TH2F *h_xj[ncent];
	TH2F *h_xj_inclusive = nullptr;
	for (int ic = 0; ic < ncent; ic++)
	{
		h_xj[ic] = (TH2F *)fproj->Get(Form("h_xj_%i", ic));
		if (ic == 0)
			h_xj_inclusive = (TH2F *)h_xj[ic]->Clone("h_xj_inclusive");
		else
			h_xj_inclusive->Add(h_xj[ic]);
	}

	int nfinal = h_xj[0]->GetNbinsY();

	for (int ipt = 0; ipt < nfinal; ipt++)
	{
		TGraphAsymmErrors *gpp = (TGraphAsymmErrors *)fpp->Get(Form("g_final_xj_statistics_%d_1", ipt));

		double pt1low = h_xj[0]->GetYaxis()->GetBinLowEdge(ipt + 1);
		double pt1high = h_xj[0]->GetYaxis()->GetBinUpEdge(ipt + 1);

		for (int i = 0; i <= ncent; i++) // i = 0..ncent-1: centralities; i = ncent: centrality-inclusive
		{
			bool isInclusive = (i == ncent);
			TH2F *hsrc = isInclusive ? h_xj_inclusive : h_xj[i];
			std::string label = isInclusive ? "Inclusive" : cent_str[i];

			hsrc->GetYaxis()->SetRange(ipt + 1, ipt + 1);
			TH1D *h1d = (TH1D *)hsrc->ProjectionX(Form("h_xj_proj_%d_pt%d", i, ipt));
			h1d->Scale(1. / h1d->Integral(), "width");

			TGraphAsymmErrors *gratio = (TGraphAsymmErrors *)gpp->Clone(Form("g_ratio_%d_pt%d", i, ipt));
			int npts = gratio->GetN();
			for (int ip = 0; ip < npts; ip++)
			{
				double x, ypp;
				gratio->GetPoint(ip, x, ypp);
				int bin = h1d->GetXaxis()->FindBin(x);
				double yhist = h1d->GetBinContent(bin);
				double ehist = h1d->GetBinError(bin);
				double eyppLow = gratio->GetErrorYlow(ip);
				double eyppHigh = gratio->GetErrorYhigh(ip);

				if (ypp == 0 || yhist == 0)
				{
					gratio->SetPoint(ip, x, 1);
					gratio->SetPointEYlow(ip, 0);
					gratio->SetPointEYhigh(ip, 0);
					continue;
				}

				double ratio = yhist / ypp;
				double relErrHist = ehist / yhist;
				double errLow = ratio * std::sqrt(relErrHist * relErrHist + (eyppLow / ypp) * (eyppLow / ypp));
				double errHigh = ratio * std::sqrt(relErrHist * relErrHist + (eyppHigh / ypp) * (eyppHigh / ypp));

				// pin the central value to 1 so only the statistical uncertainty on the ratio remains visible
				gratio->SetPoint(ip, x, 1);
				gratio->SetPointEYlow(ip, errLow);
				gratio->SetPointEYhigh(ip, errHigh);
			}

			gratio->SetTitle("");
			gratio->SetMarkerStyle(20);
			gratio->SetMarkerColor(colors[i]);
			gratio->SetLineColor(colors[i]);
			gratio->GetXaxis()->SetLimits(0.2, 1);
			gratio->GetYaxis()->SetRangeUser(0, 2.5);
			gratio->GetXaxis()->SetTitle("x_{J}");
			gratio->GetYaxis()->SetTitle("O+O/p+p");

			TCanvas *c = new TCanvas(Form("c_xjstat_pt%d_%d", ipt, i), Form("c_xjstat_pt%d_%d", ipt, i), 700, 700);

			TLegend *sphenixLeg = new TLegend(.15, .75, .45, .92);
			sphenixLeg->SetFillStyle(0);
			sphenixLeg->SetTextSize(0.032);
			sphenixLeg->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
			sphenixLeg->AddEntry("", "O+O #sqrt{s_{NN}} = 200 GeV", "");
			sphenixLeg->AddEntry("", "anti-#it{k}_{#it{t}} #it{R} = 0.4, |#eta| < 0.7", "");

			gratio->Draw("AP");

			std::string outname = isInclusive ? Form("plots/xj_statuncertainty_pt%i_inclusive.pdf", ipt)
											   : Form("plots/xj_statuncertainty_pt%i_cent%i.pdf", ipt, i);

			if (i == 0)
			{
				gAtlas->Draw("P SAME");

				TLegend *dataLeg = new TLegend(.15, .61, .92, .74);
				dataLeg->SetFillStyle(0);
				dataLeg->SetTextSize(0.022);
				dataLeg->SetNColumns(2);
				dataLeg->AddEntry(gratio, "sPHENIX statistical reach", "lep");
				dataLeg->AddEntry(gAtlas, "ATLAS result from arXiv:2606.20463", "lep");
				dataLeg->AddEntry("", "#sqrt{s_{NN}} = 200 GeV", "");
				dataLeg->AddEntry("", "#sqrt{s_{NN}} = 5.36 TeV", "");
				dataLeg->AddEntry("", "|#Delta#phi| > 3#pi/4", "");
				dataLeg->AddEntry("", "|#Delta#phi| > 7#pi/8", "");
				dataLeg->AddEntry("", Form("%.1f < p_{T1} < %.1f GeV", pt1low, pt1high), "");
				dataLeg->AddEntry("", "79 < p_{T1} < 89 GeV", "");
				dataLeg->Draw();
			}
			else
			{
				TLegend *dataLeg = new TLegend(.65, .8, .88, .88);
				dataLeg->SetFillStyle(0);
				dataLeg->SetTextSize(0.032);
				dataLeg->AddEntry(gratio, label.c_str(), "lep");
				dataLeg->Draw();
			}

			sphenixLeg->Draw();
			c->Print(outname.c_str());
		}
	}

	std::cout << "all done" << std::endl;
}
