#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

void drawResponseMatrix()
{
	SetsPhenixStyle();
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	TH3::SetDefaultSumw2();

	TFile *f = new TFile("hists/histMC.root", "READ");

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
		c->Print(Form("plots/response_matrix_cent%i.pdf", ic));
	}

	std::cout << "all done" << std::endl;
}
