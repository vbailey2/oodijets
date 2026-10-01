#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

void drawResponseMatrix()
{
	SetsPhenixStyle();
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	TH3::SetDefaultSumw2();

	TFile *f = new TFile("hists/histMC_reweight.root", "READ");

	// the data/MC weights applied on top of the MC to produce this reweighted response matrix
	TFile *frw = new TFile("hists/centrality_mbd_reweight.root", "READ");
	TH1D *h_totalcalo_et_ratio = (TH1D *)frw->Get("h_totalcalo_et_ratio");
	TH1D *h_vz_ratio = (TH1D *)frw->Get("h_vz_ratio");

	h_totalcalo_et_ratio->SetTitle(";totalcalo_et;Data/MC weight");
	h_totalcalo_et_ratio->SetMarkerStyle(20);
	h_totalcalo_et_ratio->SetMarkerColor(kBlack);
	h_totalcalo_et_ratio->SetLineColor(kBlack);

	h_vz_ratio->SetTitle(";v_{z} [cm];Data/MC weight");
	h_vz_ratio->SetMarkerStyle(20);
	h_vz_ratio->SetMarkerColor(kBlack);
	h_vz_ratio->SetLineColor(kBlack);

	TCanvas *cTotalCaloEt = new TCanvas("c_totalcalo_et_weight", "c_totalcalo_et_weight", 700, 700);
	h_totalcalo_et_ratio->Draw("PE");
	cTotalCaloEt->Print("plots/weight_totalcalo_et.pdf");

	TCanvas *cVz = new TCanvas("c_vz_weight", "c_vz_weight", 700, 700);
	h_vz_ratio->Draw("PE");
	cVz->Print("plots/weight_vz.pdf");

	TH3F *h_pt1pt2 = (TH3F *)f->Get("h_pt1pt2");
	int cent_N = h_pt1pt2->GetNbinsZ();
	int pt_N = h_pt1pt2->GetNbinsX();

	// fakes: reco dijets with no matching truth dijet (recorded separately, since
	// response[centbin].Fake() is never called in getDijets.C)
	TH3F *h_pt1pt2fake = (TH3F *)f->Get("h_pt1pt2fake");

	std::string cent_str[] = {"0-20%", "20-40%", "40-60%", "60-80%"};

	for (int ic = 0; ic < cent_N; ic++)
	{
		RooUnfoldResponse *response = (RooUnfoldResponse *)f->Get(Form("response%i", ic));

		TH2D *h_response = (TH2D *)response->Hresponse();
		h_response->SetName(Form("h_response_flat_cent%i", ic));
		h_response->SetTitle(Form("cent %s;measured global bin (p_{T,1}^{reco}#times%d + p_{T,2}^{reco});true global bin (p_{T,1}^{truth}#times%d + p_{T,2}^{truth})",
								   cent_str[ic].c_str(), pt_N, pt_N));

		TCanvas *c = new TCanvas(Form("c_response_cent%i", ic), Form("c_response_cent%i", ic), 700, 700);
		c->SetRightMargin(0.15);
		c->SetLogz();

		h_response->Draw("colz");

		TLegend *sphenixLeg = new TLegend(.15, .75, .45, .92);
		sphenixLeg->SetFillStyle(0);
		sphenixLeg->SetTextSize(0.032);
		sphenixLeg->AddEntry("", "#it{#bf{sPHENIX}} Simulation Internal", "");
		sphenixLeg->AddEntry("", "Pythia8 + HIJING O+O #sqrt{s_{NN}} = 200 GeV", "");
		sphenixLeg->AddEntry("", "anti-#it{k}_{#it{t}} #it{R} = 0.4, |#eta| < 0.7", "");
		sphenixLeg->AddEntry("", Form("%s",cent_str[ic].c_str()), "");
		sphenixLeg->Draw();

		c->Print(Form("plots/response_matrix_cent%i.pdf", ic));

		// fakes, as a 2D (p_{T,1}^{reco}, p_{T,2}^{reco}) histogram
		h_pt1pt2fake->GetZaxis()->SetRange(ic + 1, ic + 1);
		TH2D *h_fake2D = (TH2D *)h_pt1pt2fake->Project3D("yx");
		h_fake2D->SetName(Form("h_fake2D_cent%i", ic));
		h_fake2D->SetTitle(Form("cent %s fakes;p_{T,1}^{reco} [GeV];p_{T,2}^{reco} [GeV]", cent_str[ic].c_str()));

		TCanvas *cFake = new TCanvas(Form("c_fake_cent%i", ic), Form("c_fake_cent%i", ic), 700, 700);
		cFake->SetRightMargin(0.15);
		cFake->SetLogz();
		h_fake2D->Draw("colz");

		TLegend *fakeLeg = (TLegend *)sphenixLeg->Clone();
		fakeLeg->Draw();

		cFake->Print(Form("plots/response_fakes_cent%i.pdf", ic));

		// misses, as a 2D (p_{T,1}^{truth}, p_{T,2}^{truth}) histogram: there's no standalone
		// miss-only 2D histogram written to disk, since response[centbin].Miss() only feeds the
		// flattened hTrue1D/the response object, and hTrue2D%i (2D p_{T,1},p_{T,2} truth) holds
		// matched+missed truth combined. So unflatten the response matrix's matched-truth axis
		// (its Y projection, over the flattened global-bin truth axis) back into hTrue2D%i's
		// (p_{T,1},p_{T,2}) binning using the same GlobalBin mapping as getDijets.C
		// (gTrue = (ix-1)*pt_N + iy, ix = p_{T,1} bin, iy = p_{T,2} bin), then subtract that
		// matched-truth 2D histogram from the total truth to leave the misses alone.
		TH2D *h_true2D = (TH2D *)f->Get(Form("hTrue2D%i", ic));
		TH1D *h_matchedtrue1D = (TH1D *)h_response->ProjectionY(Form("h_matchedtrue1D_cent%i", ic));
		TH2D *h_matchedtrue2D = (TH2D *)h_true2D->Clone(Form("h_matchedtrue2D_cent%i", ic));
		h_matchedtrue2D->Reset();
		for (int g = 1; g <= h_matchedtrue1D->GetNbinsX(); g++)
		{
			int ix = (g - 1) / pt_N + 1;
			int iy = (g - 1) % pt_N + 1;
			h_matchedtrue2D->SetBinContent(ix, iy, h_matchedtrue1D->GetBinContent(g));
			h_matchedtrue2D->SetBinError(ix, iy, h_matchedtrue1D->GetBinError(g));
		}
		TH2D *h_miss2D = (TH2D *)h_true2D->Clone(Form("h_miss2D_cent%i", ic));
		h_miss2D->Add(h_matchedtrue2D, -1);
		h_miss2D->SetTitle(Form("cent %s misses;p_{T,1}^{truth} [GeV];p_{T,2}^{truth} [GeV]", cent_str[ic].c_str()));

		TCanvas *cMiss = new TCanvas(Form("c_miss_cent%i", ic), Form("c_miss_cent%i", ic), 700, 700);
		cMiss->SetRightMargin(0.15);
		cMiss->SetLogz();
		h_miss2D->Draw("colz");

		TLegend *missLeg = (TLegend *)sphenixLeg->Clone();
		missLeg->Draw();

		cMiss->Print(Form("plots/response_misses_cent%i.pdf", ic));
	}

	std::cout << "all done" << std::endl;
}
