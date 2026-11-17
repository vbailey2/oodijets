#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

// Truth (pT1, pT2) prior reweighting of the MC, following danjlis/dijet-unfolding (createResponse_*: the
// "PRIMER" pass and h_flatreweight_pt1pt2).
//
// Inputs: a data unfold done with the response matrix *without* the prior weight (the "primer" unfold),
// which contains both the unfolded data (hUnfolded2D_cent%d) and the MC truth used to build that
// response (hTrue2D%d), both symmetric in (pT1, pT2).
//
// Per centrality and per truth cell (i, j), with k = max(i, j) the leading-pT bin of the cell:
//   U_ij -> U_ij / sum_{cells with max = k} U     (only when perLeadingBin is true; otherwise / sum_all U)
//   T_ij -> T_ij / sum_{cells with max = k} T
//   w_ij = U_ij / T_ij, clamped to [wmin, wmax]
// With perLeadingBin (the reference's choice) each leading-pT row is normalized separately, so only the
// pT2 | pT1 shape (the xJ shape) of the MC prior is corrected, not its leading-jet spectrum.
// Cells with no unfolded or truth entries take the weight of the nearest filled cell in the same
// leading-pT row, walking from the diagonal (xJ = 1) down to low pT2 (as in the reference); cells
// still without a weight are set to 1.
//
// Weights are only derived inside the reco acceptance (leading truth bin >= minLeadBin, subleading truth
// bin >= minSubBin; defaults = measure tier pt_bins[10], pt_bins[4] -> 1-based bins 11, 5, which is now the
// reco selection for both data and MC in getDijets.C). Truth cells below
// the reco cuts are reached only through migrations, so the unfold there is barely constrained and the
// primer ratio is noise (e.g. it hit the wmax clamp at pT2 ~ 6 GeV, which then ran away with iterations).
// Those cells keep weight 1 and are excluded from the row normalizations.
void getPriorReweighting(const char *unfoldfile = "hists/hist-unfoldedData_primer.root",
						 const char *outfile = "hists/prior_reweight.root",
						 bool perLeadingBin = true, double wmin = 0.1, double wmax = 5.,
						 int minLeadBin = 11, int minSubBin = 5)
{
	TFile *f = TFile::Open(unfoldfile, "READ");
	if (!f || f->IsZombie())
	{
		std::cerr << "ERROR: cannot open " << unfoldfile << std::endl;
		return;
	}

	SetsPhenixStyle();
	gStyle->SetPalette(kBird);

	std::vector<TH2D *> hw;
	for (int ic = 0;; ic++)
	{
		TH2D *hU = (TH2D *)f->Get(Form("hUnfolded2D_cent%d", ic));
		TH2D *hT = (TH2D *)f->Get(Form("hTrue2D%d", ic));
		if (!hU || !hT)
			break;

		const int n = hU->GetNbinsX();
		TH2D *h = (TH2D *)hU->Clone(Form("h_prior_reweight_cent%d", ic));
		h->SetDirectory(nullptr);
		h->Reset();
		h->SetTitle(";#it{p}_{T,1}^{truth} [GeV];#it{p}_{T,2}^{truth} [GeV];prior weight");

		// normalizations: per leading-pT row k = max(i, j), or one overall sum
		std::vector<double> sumU(n + 1, 0.), sumT(n + 1, 0.);
		double totU = 0, totT = 0;
		for (int i = 1; i <= n; i++)
			for (int j = 1; j <= n; j++)
			{
				int k = std::max(i, j);
				if (k < minLeadBin || std::min(i, j) < minSubBin)
					continue;
				double u = hU->GetBinContent(i, j), t = hT->GetBinContent(i, j);
				if (u > 0)
				{
					sumU[k] += u;
					totU += u;
				}
				if (t > 0)
				{
					sumT[k] += t;
					totT += t;
				}
			}

		for (int i = 1; i <= n; i++)
			for (int j = 1; j <= n; j++)
			{
				int k = std::max(i, j);
				if (k < minLeadBin || std::min(i, j) < minSubBin)
					continue;
				double nu = perLeadingBin ? sumU[k] : totU;
				double nt = perLeadingBin ? sumT[k] : totT;
				double u = hU->GetBinContent(i, j), t = hT->GetBinContent(i, j);
				if (nu <= 0 || nt <= 0 || u <= 0 || t <= 0)
					continue;
				double w = (u / nu) / (t / nt);
				h->SetBinContent(i, j, std::min(std::max(w, wmin), wmax));
			}

		// fill empty cells from the nearest filled cell in the same leading-pT row (diagonal -> low pT2),
		// keeping the matrix symmetric; anything left over gets weight 1
		for (int k = 1; k <= n; k++)
		{
			double last = 0;
			for (int j = k; j >= 1; j--)
			{
				if (k < minLeadBin || j < minSubBin)
				{
					h->SetBinContent(k, j, 1.);
					h->SetBinContent(j, k, 1.);
					continue;
				}
				double w = h->GetBinContent(k, j);
				if (w <= 0)
					w = h->GetBinContent(j, k);
				if (w <= 0)
					w = last;
				else
					last = w;
				if (w <= 0)
					w = 1.;
				h->SetBinContent(k, j, w);
				h->SetBinContent(j, k, w);
			}
		}
		hw.push_back(h);

		TCanvas *c = new TCanvas(Form("c_prior_cent%d", ic), "", 700, 600);
		c->SetRightMargin(0.18);
		c->SetLogx();
		c->SetLogy();
		h->SetMinimum(wmin);
		h->SetMaximum(wmax);
		h->Draw("colz");
		TLegend *leg = new TLegend(.17, .78, .5, .92);
		leg->SetFillStyle(0);
		leg->AddEntry("", "#it{#bf{sPHENIX}} Simulation Internal", "");
		leg->AddEntry("", Form("Prior reweighting, cent %d", ic), "");
		leg->Draw();
		c->Print(Form("plots/prior_reweight_cent%d.pdf", ic));
	}
	f->Close();

	if (hw.empty())
	{
		std::cerr << "ERROR: no hUnfolded2D_cent%d / hTrue2D%d found in " << unfoldfile << std::endl;
		return;
	}

	TFile *fout = new TFile(outfile, "RECREATE");
	for (auto *h : hw)
		h->Write();
	fout->Close();
	std::cout << "wrote " << hw.size() << " prior weight maps to " << outfile << std::endl;
}
