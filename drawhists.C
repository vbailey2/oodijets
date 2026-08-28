#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

void drawhists(int isunfold = 0)
{

	SetsPhenixStyle();
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	TH3::SetDefaultSumw2();

	string fname;
	string projname;
	if (isunfold)
	{
		fname = Form("hists/hist-unfoldedData_iter%d.root", isunfold);
		projname = Form("hists/projections_unfold_iter%d.root", isunfold);
	}
	else
	{
		fname = "hists/hist-full.root";
		projname = "hists/projections.root";
	}
	TFile *f = new TFile(fname.c_str(), "READ");
	TFile *fproj = new TFile(projname.c_str(), "READ");
	TFile *fpp = nullptr;
	if (isunfold)
		fpp = new TFile("hists/final_plots_pp_r04.root", "READ");

	TCanvas *c = new TCanvas("c", "c", 700, 700);

	TCanvas *c2 = new TCanvas("c2", "c2", 700, 700);

	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
	pad1->SetBottomMargin(0.01);
	pad1->SetLeftMargin(0.12);
	pad1->Draw();

	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
	pad2->SetTopMargin(0.03);
	pad2->SetBottomMargin(0.35);
	pad2->SetLeftMargin(0.12);
	pad2->Draw();

	TLegend *leg = new TLegend(.15, .75, .4, .92);
	leg->SetFillStyle(0);
	leg->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
	// leg->AddEntry("","Pythia + HIJING","");
	leg->AddEntry("", "O+O #sqrt{s_{NN}} = 200 GeV", "");
	leg->AddEntry("", "anti-#it{k}_{#it{t}} #it{R} = 0.4, |#eta| < 0.7", "");

	TLegend *cleg = new TLegend(.15, .6, .4, .75);
	cleg->SetFillStyle(0);

	TLegend *hleg = new TLegend(.7, .75, .9, .92);
	hleg->SetFillStyle(0);

	int colors[] = {1, 2, 4, kGreen + 2, kViolet, kCyan, kOrange + 2, kMagenta + 2, kAzure - 2};

	TH3F *h_dphi = (TH3F *)f->Get("h_dphi");

	TH3F *h_xjnoproj3D = (TH3F *)f->Get("h_xj");
	h_xjnoproj3D->SetName("h_xjnoproj3D");
	int npt = h_xjnoproj3D->GetNbinsY();
	int ncent = h_xjnoproj3D->GetNbinsZ();

	std::string cent_str[] = {"0-20%", "20-40%", "40-60%", "60-80%"};

	TH2F *h_xj[ncent];
	TH2F *h_xjnoproj[ncent];
	h_xj[0] = (TH2F *)fproj->Get("h_xj_0");
	const int n_xj = h_xj[0]->GetNbinsX();
	const int n_final = h_xj[0]->GetNbinsY();
	for (int i = 0; i < ncent; i++)
	{
		h_xj[i] = (TH2F *)fproj->Get(Form("h_xj_%i", i));
		if (!isunfold)
		{
			// put the non-projected xj in the new pt bins
			h_xjnoproj[i] = (TH2F *)h_xj[i]->Clone();
			h_xjnoproj[i]->Reset();
			h_xjnoproj[i]->SetName(Form("h_xjnoproj%i", i));
			double bincount[3][19] = {};
			double binerr[3][19] = {};
			for (int ifinal = 0; ifinal < h_xj[i]->GetNbinsY(); ifinal++)
			{
				int finlow = h_xj[i]->GetYaxis()->GetBinLowEdge(ifinal + 1);
				int finhigh = h_xj[i]->GetYaxis()->GetBinLowEdge(ifinal + 2);
				for (int ipt = 0; ipt < h_xjnoproj3D->GetNbinsY(); ipt++)
				{
					int ptlow = h_xjnoproj3D->GetYaxis()->GetBinLowEdge(ipt + 1);
					for (int ix = 0; ix < h_xj[i]->GetNbinsX(); ix++)
					{
						// check this pt bin is in this combined final pt range
						if (ptlow >= finlow && ptlow < finhigh)
						{
							// if it is then add to that bin count
							bincount[ifinal][ix] += h_xjnoproj3D->GetBinContent(ix + 1, ipt + 1, i + 1);
							binerr[ifinal][ix] += h_xjnoproj3D->GetBinError(ix + 1, ipt + 1, i + 1) * h_xjnoproj3D->GetBinError(ix + 1, ipt + 1, i + 1);
						}
					}
				}
				for (int ix = 0; ix < h_xj[i]->GetNbinsX(); ix++)
				{
					h_xjnoproj[i]->SetBinContent(ix + 1, ifinal + 1, bincount[ifinal][ix]);
					h_xjnoproj[i]->SetBinError(ix + 1, ifinal + 1, std::sqrt(binerr[ifinal][ix]));
				}
			}
		}
	}

	TH1F *h_xj1D[ncent][npt];
	TH1F *h_xjnoproj1D[ncent][npt];
	TH1F *h_dphi1D[ncent][npt];

	for (int ipt = 0; ipt < n_final; ipt++)
	{
		if (isunfold)
		{
			c2->cd();
			pad1->cd();
			pad1->Clear();
		}
		else
			c->cd();
		if (!isunfold)
			cleg->AddEntry("", Form("%2.1f < p_{T}^{calib} < %2.1f GeV", h_xj[0]->GetYaxis()->GetBinLowEdge(ipt + 1), h_xj[0]->GetYaxis()->GetBinLowEdge(ipt + 2)), "");
		else
			cleg->AddEntry("", Form("%2.1f < p_{T} < %2.1f GeV", h_xj[0]->GetYaxis()->GetBinLowEdge(ipt + 1), h_xj[0]->GetYaxis()->GetBinLowEdge(ipt + 2)), "");
		for (int icent = 0; icent < ncent; icent++)
		{
			h_xj[icent]->GetYaxis()->SetRange(ipt + 1, ipt + 1);
			h_xj1D[icent][ipt] = (TH1F *)h_xj[icent]->ProjectionX();
			h_xj1D[icent][ipt]->SetName(Form("h_xj_cent%i_pt%i", icent, ipt));
			h_xj1D[icent][ipt]->GetXaxis()->SetRangeUser(0.2, 1);
			h_xj1D[icent][ipt]->Scale(1. / h_xj1D[icent][ipt]->Integral(), "width");

			h_xj1D[icent][ipt]->SetMarkerColor(colors[icent]);
			h_xj1D[icent][ipt]->SetLineColor(colors[icent]);
			h_xj1D[icent][ipt]->GetYaxis()->SetRangeUser(0, 5);
			h_xj1D[icent][ipt]->GetXaxis()->SetTitle("x_{J}");
			if (icent == 0)
				h_xj1D[icent][ipt]->Draw();
			else
				h_xj1D[icent][ipt]->Draw("SAME");
			hleg->AddEntry(h_xj1D[icent][ipt], cent_str[icent].c_str(), "p");
		}

		TGraphAsymmErrors *gpp = nullptr;
		if (isunfold)
		{
			gpp = (TGraphAsymmErrors *)fpp->Get(Form("g_final_xj_statistics_%d_1", ipt));
			gpp->SetMarkerStyle(24);
			gpp->SetMarkerColor(kBlack);
			gpp->SetLineColor(kBlack);
			gpp->Draw("P SAME");
			hleg->AddEntry(gpp, "pp", "p");
		}

		leg->Draw();
		cleg->Draw();
		hleg->Draw();
		string sunfold;
		if (isunfold)
			sunfold = Form("unfold_iter%d", isunfold);
		else
			sunfold = "";

		if (isunfold)
		{
			pad2->cd();
			pad2->Clear();
			for (int icent = 0; icent < ncent; icent++)
			{
				TGraphAsymmErrors *gratio = (TGraphAsymmErrors *)gpp->Clone(Form("g_ratio_cent%i_pt%i", icent, ipt));
				int npts = gratio->GetN();
				for (int ip = 0; ip < npts; ip++)
				{
					double x, ypp;
					gratio->GetPoint(ip, x, ypp);
					int bin = h_xj1D[icent][ipt]->GetXaxis()->FindBin(x);
					double yhist = h_xj1D[icent][ipt]->GetBinContent(bin);
					double ehist = h_xj1D[icent][ipt]->GetBinError(bin);
					double eyppLow = gratio->GetErrorYlow(ip);
					double eyppHigh = gratio->GetErrorYhigh(ip);

					if (ypp == 0 || yhist == 0)
					{
						gratio->SetPoint(ip, x, 0);
						gratio->SetPointEYlow(ip, 0);
						gratio->SetPointEYhigh(ip, 0);
						continue;
					}

					double ratio = yhist / ypp;
					double relErrHist = ehist / yhist;
					double errLow = ratio * std::sqrt(relErrHist * relErrHist + (eyppLow / ypp) * (eyppLow / ypp));
					double errHigh = ratio * std::sqrt(relErrHist * relErrHist + (eyppHigh / ypp) * (eyppHigh / ypp));

					gratio->SetPoint(ip, x, ratio);
					gratio->SetPointEYlow(ip, errLow);
					gratio->SetPointEYhigh(ip, errHigh);
				}

				gratio->SetMarkerColor(colors[icent]);
				gratio->SetLineColor(colors[icent]);
				gratio->SetMarkerStyle(20);

				gratio->SetTitle("");
				gratio->GetYaxis()->SetRangeUser(0, 2);
				gratio->GetXaxis()->SetLimits(0.2, 1);
				gratio->GetXaxis()->SetTitle("x_{J}");
				gratio->GetXaxis()->SetTitleSize(0.12);
				gratio->GetXaxis()->SetTitleOffset(1.0);
				gratio->GetXaxis()->SetLabelSize(0.10);
				gratio->GetXaxis()->SetTickLength(0.07);

				gratio->GetYaxis()->SetTitle("AA/pp");
				gratio->GetYaxis()->SetTitleSize(0.12);
				gratio->GetYaxis()->SetTitleOffset(0.45);
				gratio->GetYaxis()->SetLabelSize(0.10);
				gratio->GetYaxis()->SetNdivisions(505);

				if (icent == 0)
					gratio->Draw("AP");
				else
					gratio->Draw("P SAME");
			}

			TLine *ppline = new TLine(0.2, 1.0, 1.0, 1.0);
			ppline->SetLineColor(kBlack);
			ppline->SetLineStyle(2);
			ppline->Draw("SAME");

			c2->Print(Form("plots/xj_pt%i_%s.pdf", ipt, sunfold.c_str()));
		}
		else
		{
			c->Print(Form("plots/xj_pt%i_%s.pdf", ipt, sunfold.c_str()));
			c->Clear();
		}
		hleg->Clear();
		cleg->Clear();

		if (!isunfold)
		{
			c2->cd();
			for (int icent = 0; icent < ncent; icent++)
			{
				pad1->cd();
				h_xjnoproj[icent]->GetYaxis()->SetRange(ipt + 1, ipt + 1);
				h_xjnoproj1D[icent][ipt] = (TH1F *)h_xjnoproj[icent]->ProjectionX();
				h_xjnoproj1D[icent][ipt]->SetName(Form("h_xjnoproj_cent%i_pt%i", icent, ipt));

				h_xjnoproj1D[icent][ipt]->Scale(1. / h_xjnoproj1D[icent][ipt]->Integral(), "width");

				h_xjnoproj1D[icent][ipt]->SetMarkerColor(colors[icent]);
				h_xjnoproj1D[icent][ipt]->SetMarkerStyle(25);
				h_xjnoproj1D[icent][ipt]->SetLineColor(colors[icent]);
				h_xjnoproj1D[icent][ipt]->GetYaxis()->SetRangeUser(0, 5);
				h_xjnoproj1D[icent][ipt]->GetXaxis()->SetRangeUser(0.2, 1);
				h_xjnoproj1D[icent][ipt]->GetXaxis()->SetTitle("x_{J}");
				h_xj1D[icent][ipt]->Draw();
				hleg->AddEntry(h_xj1D[icent][ipt], Form("Proj. %s", cent_str[icent].c_str()), "p");

				h_xjnoproj1D[icent][ipt]->Draw("SAME");
				hleg->AddEntry(h_xjnoproj1D[icent][ipt], Form("No proj. %s", cent_str[icent].c_str()), "p");
				hleg->Draw();
				cleg->Draw();
				leg->Draw();

				pad2->cd();
				TH1D *hRatio = (TH1D *)h_xj1D[icent][ipt]->Clone("hRatio");
				hRatio->Divide(h_xjnoproj1D[icent][ipt]);
				hRatio->SetTitle("");

				hRatio->SetMinimum(0.);
				hRatio->SetMaximum(2.);
				hRatio->GetXaxis()->SetTitle("x_{J}");
				hRatio->GetXaxis()->SetTitleSize(0.12);
				hRatio->GetXaxis()->SetTitleOffset(1.0);
				hRatio->GetXaxis()->SetLabelSize(0.10);
				hRatio->GetXaxis()->SetTickLength(0.07);

				hRatio->GetYaxis()->SetTitle("Ratio");
				hRatio->GetYaxis()->SetTitleSize(0.12);
				hRatio->GetYaxis()->SetTitleOffset(0.45);
				hRatio->GetYaxis()->SetLabelSize(0.10);
				hRatio->GetYaxis()->SetNdivisions(505);

				hRatio->Draw();

				int firstBin = hRatio->GetXaxis()->GetFirst();
				int lastBin = hRatio->GetXaxis()->GetLast();

				double xmin = hRatio->GetXaxis()->GetBinLowEdge(firstBin);
				double xmax = hRatio->GetXaxis()->GetBinUpEdge(lastBin);
				TLine *line = new TLine(xmin, 1.0,
										xmax, 1.0);
				line->SetLineColor(kBlack);
				line->SetLineStyle(2);
				line->Draw("SAME");

				c2->Print(Form("plots/xj_projclosure_pt%i_cent%i.pdf", ipt, icent));
				hleg->Clear();
			}
			cleg->Clear();
		}
	}
	if (!isunfold)
	{
		c->cd();
		for (int ipt = 0; ipt < npt; ipt++)
		{
			for (int icent = 0; icent < ncent; icent++)
			{
				h_dphi->GetYaxis()->SetRange(ipt + 1, ipt + 1);
				h_dphi->GetZaxis()->SetRange(icent + 1, icent + 1);
				h_dphi1D[icent][ipt] = (TH1F *)h_dphi->Project3D("x");
				h_dphi1D[icent][ipt]->SetName(Form("h_dphi_cent%i_pt%i", icent, ipt));
				h_dphi1D[icent][ipt]->Scale(1. / h_dphi1D[icent][ipt]->Integral());
				h_dphi1D[icent][ipt]->SetMarkerColor(colors[icent]);
				h_dphi1D[icent][ipt]->SetLineColor(colors[icent]);
				h_dphi1D[icent][ipt]->GetYaxis()->SetRangeUser(0, 1);
				h_dphi1D[icent][ipt]->GetXaxis()->SetTitle("#Delta#phi");
				if (icent == 0)
					h_dphi1D[icent][ipt]->Draw();
				else
					h_dphi1D[icent][ipt]->Draw("SAME");
				hleg->AddEntry(h_dphi1D[icent][ipt], cent_str[icent].c_str(), "p");
			}
			leg->Draw();
			cleg->Draw();
			hleg->Draw();
			c->Print(Form("plots/dphi_pt%i.pdf", ipt));
			c->Clear();
			cleg->Clear();
			hleg->Clear();
		}
	}
}
