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

		// fakes: reco dijets with no matching truth dijet (recorded separately, since
		// response[centbin].Fake() is never called in getDijets.C), flattened over the
		// same (p_{T,1}#times%d + p_{T,2}) global-bin axis as the response matrix's measured axis
		TH1D *h_fake = (TH1D *)f->Get(Form("hFake1D%i", ic));
		h_fake->SetName(Form("h_fake_cent%i", ic));
		h_fake->SetTitle(Form("cent %s;measured global bin (p_{T,1}^{reco}#times%d + p_{T,2}^{reco});fakes",
							   cent_str[ic].c_str(), pt_N));
		h_fake->SetLineColor(kBlack);
		h_fake->SetMarkerColor(kBlack);
		h_fake->SetMarkerStyle(20);

		TCanvas *cFake = new TCanvas(Form("c_fake_cent%i", ic), Form("c_fake_cent%i", ic), 700, 700);
		cFake->SetLogy();
		h_fake->Draw("PE");

		TLegend *fakeLeg = (TLegend *)sphenixLeg->Clone();
		fakeLeg->Draw();

		cFake->Print(Form("plots/response_fakes_cent%i.pdf", ic));

		// misses: truth dijets with no matching reco dijet, i.e. the truth entries filled via
		// response[centbin].Miss() rather than response[centbin].Fill(); computed as the total
		// truth distribution (hTrue1D, matched+missed) minus the matched-truth projection of
		// the response matrix (its true-axis, i.e. Y, projection)
		TH1D *h_true = (TH1D *)f->Get(Form("hTrue1D%i", ic));
		TH1D *h_matchedtrue = (TH1D *)h_response->ProjectionY(Form("h_matchedtrue_cent%i", ic));
		TH1D *h_miss = (TH1D *)h_true->Clone(Form("h_miss_cent%i", ic));
		h_miss->Add(h_matchedtrue, -1);
		h_miss->SetTitle(Form("cent %s;true global bin (p_{T,1}^{truth}#times%d + p_{T,2}^{truth});misses",
							   cent_str[ic].c_str(), pt_N));
		h_miss->SetLineColor(kBlack);
		h_miss->SetMarkerColor(kBlack);
		h_miss->SetMarkerStyle(20);

		TCanvas *cMiss = new TCanvas(Form("c_miss_cent%i", ic), Form("c_miss_cent%i", ic), 700, 700);
		cMiss->SetLogy();
		h_miss->Draw("PE");

		TLegend *missLeg = (TLegend *)sphenixLeg->Clone();
		missLeg->Draw();

		cMiss->Print(Form("plots/response_misses_cent%i.pdf", ic));
	}

	std::cout << "all done" << std::endl;
}
