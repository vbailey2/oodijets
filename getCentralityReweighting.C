#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

void drawDataMCRatio(TH1D *h_data, TH1D *h_mc, TH1D *h_ratio, const char *xtitle, const char *outname)
{
	TCanvas *c = new TCanvas(outname, outname, 700, 700);

	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
	pad1->SetBottomMargin(0.01);
	pad1->SetLeftMargin(0.12);
	pad1->Draw();

	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
	pad2->SetTopMargin(0.03);
	pad2->SetBottomMargin(0.35);
	pad2->SetLeftMargin(0.12);
	pad2->Draw();

	pad1->cd();
	h_data->SetTitle("");
	h_data->SetMarkerStyle(20);
	h_data->SetMarkerColor(kBlack);
	h_data->SetLineColor(kBlack);
	h_data->GetXaxis()->SetTitle(xtitle);
	h_data->Draw("p");

	h_mc->SetMarkerStyle(21);
	h_mc->SetMarkerColor(kRed);
	h_mc->SetLineColor(kRed);
	h_mc->Draw("p same");

	TLegend *leg = new TLegend(.65, .7, .88, .88);
	leg->SetFillStyle(0);
	leg->AddEntry(h_data, "Data", "p");
	leg->AddEntry(h_mc, "MC", "p");
	leg->Draw();

	pad2->cd();
	h_ratio->SetTitle("");
	h_ratio->SetMinimum(0.);
	h_ratio->SetMaximum(2.);
	h_ratio->SetMarkerStyle(20);
	h_ratio->SetMarkerColor(kBlack);
	h_ratio->SetLineColor(kBlack);

	h_ratio->GetXaxis()->SetTitle(xtitle);
	h_ratio->GetXaxis()->SetTitleSize(0.12);
	h_ratio->GetXaxis()->SetTitleOffset(1.0);
	h_ratio->GetXaxis()->SetLabelSize(0.10);
	h_ratio->GetXaxis()->SetTickLength(0.07);

	h_ratio->GetYaxis()->SetTitle("Data/MC");
	h_ratio->GetYaxis()->SetTitleSize(0.12);
	h_ratio->GetYaxis()->SetTitleOffset(0.45);
	h_ratio->GetYaxis()->SetLabelSize(0.10);
	h_ratio->GetYaxis()->SetNdivisions(505);

	h_ratio->Draw("p");

	int firstBin = h_ratio->GetXaxis()->GetFirst();
	int lastBin = h_ratio->GetXaxis()->GetLast();
	double xmin = h_ratio->GetXaxis()->GetBinLowEdge(firstBin);
	double xmax = h_ratio->GetXaxis()->GetBinUpEdge(lastBin);
	TLine *line = new TLine(xmin, 1.0, xmax, 1.0);
	line->SetLineColor(kBlack);
	line->SetLineStyle(2);
	line->Draw("SAME");

	c->Print(Form("plots/%s.pdf", outname));
}

void getCentralityReweighting()
{
	SetsPhenixStyle();
	TH1::SetDefaultSumw2();

	TFile *fData = new TFile("hists/hist-full.root", "READ");
	TFile *fMC = new TFile("hists/histMC.root", "READ");

	TH1D *h_centrality_data = (TH1D *)fData->Get("h_centrality");
	h_centrality_data->SetName("h_centrality_data");
	TH1D *h_centrality_mc = (TH1D *)fMC->Get("h_centrality");
	h_centrality_mc->SetName("h_centrality_mc");
	TH1D *h_mbd_charge_sum_data = (TH1D *)fData->Get("h_mbd_charge_sum");
	h_mbd_charge_sum_data->SetName("h_mbd_charge_sum_data");
	TH1D *h_mbd_charge_sum_mc = (TH1D *)fMC->Get("h_mbd_charge_sum");
	h_mbd_charge_sum_mc->SetName("h_mbd_charge_sum_mc");
	TH1D *h_totalcalo_et_data = (TH1D *)fData->Get("h_totalcalo_et");
	h_totalcalo_et_data->SetName("h_totalcalo_et_data");
	TH1D *h_totalcalo_et_mc = (TH1D *)fMC->Get("h_totalcalo_et");
	h_totalcalo_et_mc->SetName("h_totalcalo_et_mc");

	h_centrality_data->Scale(1. / h_centrality_data->Integral(), "width");
	h_centrality_mc->Scale(1. / h_centrality_mc->Integral(), "width");
	h_mbd_charge_sum_data->Scale(1. / h_mbd_charge_sum_data->Integral(), "width");
	h_mbd_charge_sum_mc->Scale(1. / h_mbd_charge_sum_mc->Integral(), "width");
	h_totalcalo_et_data->Scale(1. / h_totalcalo_et_data->Integral(), "width");
	h_totalcalo_et_mc->Scale(1. / h_totalcalo_et_mc->Integral(), "width");

	TH1D *h_centrality_ratio = (TH1D *)h_centrality_data->Clone("h_centrality_ratio");
	h_centrality_ratio->Divide(h_centrality_mc);

	TH1D *h_mbd_charge_sum_ratio = (TH1D *)h_mbd_charge_sum_data->Clone("h_mbd_charge_sum_ratio");
	h_mbd_charge_sum_ratio->Divide(h_mbd_charge_sum_mc);

	TH1D *h_totalcalo_et_ratio = (TH1D *)h_totalcalo_et_data->Clone("h_totalcalo_et_ratio");
	h_totalcalo_et_ratio->Divide(h_totalcalo_et_mc);

	drawDataMCRatio(h_centrality_data, h_centrality_mc, h_centrality_ratio, "Centrality bin", "datamc_centrality");
	drawDataMCRatio(h_mbd_charge_sum_data, h_mbd_charge_sum_mc, h_mbd_charge_sum_ratio, "mbd_charge_sum", "datamc_mbd_charge_sum");
	drawDataMCRatio(h_totalcalo_et_data, h_totalcalo_et_mc, h_totalcalo_et_ratio, "totalcalo_et", "datamc_totalcalo_et");

	TFile *fout = new TFile("hists/centrality_mbd_reweight.root", "RECREATE");
	h_centrality_ratio->Write();
	h_mbd_charge_sum_ratio->Write();
	h_totalcalo_et_ratio->Write();
	std::cout << "all done" << std::endl;
}
