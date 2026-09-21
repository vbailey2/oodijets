#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

void drawXJStatUncertainty()
{
	SetsPhenixStyle();
	TH1::SetDefaultSumw2();
	TH3::SetDefaultSumw2();

	TFile *f = new TFile("hists/hist-full.root", "READ");
	TH3F *h_xj = (TH3F *)f->Get("h_xj");

	int npt = h_xj->GetNbinsY();
	int ncent = h_xj->GetNbinsZ();
	std::string cent_str[] = {"0-20%", "20-40%", "40-60%", "60-80%"};
	int colors[] = {1, 2, 4, kGreen + 2};

	for (int ipt = 0; ipt < npt; ipt++)
	{
		// peripheral (most peripheral centrality bin) reference for this pt bin
		h_xj->GetYaxis()->SetRange(ipt + 1, ipt + 1);
		h_xj->GetZaxis()->SetRange(ncent, ncent);
		TH1D *h_periph = (TH1D *)h_xj->Project3D("x");
		h_periph->SetName(Form("h_xj_periph_pt%i", ipt));
		if (h_periph->Integral() <= 0)
		{
			std::cout << "pt bin " << ipt << ": peripheral slice is empty, skipping" << std::endl;
			continue;
		}
		h_periph->Scale(1. / h_periph->Integral(), "width");

		// first pass: build all centrality ratios for this pt bin and track the largest error
		TH1D *h_ratio[ncent];
		double maxerr = 0;
		for (int icent = 0; icent < ncent; icent++)
		{
			h_xj->GetYaxis()->SetRange(ipt + 1, ipt + 1);
			h_xj->GetZaxis()->SetRange(icent + 1, icent + 1);
			TH1D *h_cent = (TH1D *)h_xj->Project3D("x");
			h_cent->SetName(Form("h_xj_cent%i_pt%i", icent, ipt));
			if (h_cent->Integral() <= 0)
			{
				std::cout << "pt bin " << ipt << ", cent bin " << icent << ": slice is empty, skipping" << std::endl;
				h_ratio[icent] = nullptr;
				continue;
			}
			h_cent->Scale(1. / h_cent->Integral(), "width");

			h_ratio[icent] = (TH1D *)h_cent->Clone(Form("h_xj_ratio_cent%i_pt%i", icent, ipt));
			h_ratio[icent]->Divide(h_periph);

			// zero out the central value so only the statistical uncertainty on the ratio remains visible
			for (int ib = 1; ib <= h_ratio[icent]->GetNbinsX(); ib++)
			{
				double err = h_ratio[icent]->GetBinError(ib);
				h_ratio[icent]->SetBinContent(ib, 0);
				if (h_ratio[icent]->GetXaxis()->GetBinCenter(ib) >= 0.2 && err > maxerr)
					maxerr = err;
			}
		}

		if (maxerr <= 0)
			maxerr = 0.1;

		TCanvas *c = new TCanvas(Form("c_xjstat_pt%i", ipt), Form("c_xjstat_pt%i", ipt), 700, 700);
		TLegend *leg = new TLegend(.65, .7, .88, .88);
		leg->SetFillStyle(0);

		bool drawnany = false;
		for (int icent = 0; icent < ncent; icent++)
		{
			if (!h_ratio[icent])
				continue;

			h_ratio[icent]->SetTitle("");
			h_ratio[icent]->SetMarkerStyle(20);
			h_ratio[icent]->SetMarkerColor(colors[icent]);
			h_ratio[icent]->SetLineColor(colors[icent]);
			h_ratio[icent]->GetXaxis()->SetRangeUser(0.2, 1);
			h_ratio[icent]->GetYaxis()->SetRangeUser(-1.3 * maxerr, 1.3 * maxerr);
			h_ratio[icent]->GetXaxis()->SetTitle("x_{J}");
			h_ratio[icent]->GetYaxis()->SetTitle(Form("stat. uncertainty on ratio to %s", cent_str[ncent - 1].c_str()));

			if (!drawnany)
			{
				h_ratio[icent]->Draw("PE");
				drawnany = true;
			}
			else
				h_ratio[icent]->Draw("PE SAME");
			leg->AddEntry(h_ratio[icent], cent_str[icent].c_str(), "lep");
		}

		leg->Draw();
		c->Print(Form("plots/xj_statuncertainty_pt%i.pdf", ipt));
	}

	std::cout << "all done" << std::endl;
}
