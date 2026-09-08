void getCentralityReweighting()
{
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

	h_centrality_data->Scale(1. / h_centrality_data->Integral(), "width");
	h_centrality_mc->Scale(1. / h_centrality_mc->Integral(), "width");
	h_mbd_charge_sum_data->Scale(1. / h_mbd_charge_sum_data->Integral(), "width");
	h_mbd_charge_sum_mc->Scale(1. / h_mbd_charge_sum_mc->Integral(), "width");

	TH1D *h_centrality_ratio = (TH1D *)h_centrality_data->Clone("h_centrality_ratio");
	h_centrality_ratio->Divide(h_centrality_mc);

	TH1D *h_mbd_charge_sum_ratio = (TH1D *)h_mbd_charge_sum_data->Clone("h_mbd_charge_sum_ratio");
	h_mbd_charge_sum_ratio->Divide(h_mbd_charge_sum_mc);

	TFile *fout = new TFile("hists/centrality_mbd_reweight.root", "RECREATE");
	h_centrality_ratio->Write();
	h_mbd_charge_sum_ratio->Write();
	std::cout << "all done" << std::endl;
}
