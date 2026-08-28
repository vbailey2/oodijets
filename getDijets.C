#include <cmath>

const int nCentBins = 4;
double centbins[nCentBins + 1] = {0, 20, 40, 60, 80};
double mbdcharge_boundaries[nCentBins + 1] = {200, 80.3746, 47.2647, 26.1444, 12.8285};
double mbdcharge_boundaries_mc[nCentBins + 1] = {200, 80.4934, 47.2058, 25.9975, 12.3982};
std::string cent_str[nCentBins] = {"0-20%", "20-40%", "40-60%", "60-80%"};

inline const char *DefaultMBMbdScaleFile()
{
	return "/gpfs/mnt/gpfs02/sphenix/user/jpark4/Analysis/Oxygen/sPHENIX-OxygenAnalysis/Corrections/mb_mbd_scale_fit.root";
}

inline double GetMbdChargeSumScaleMC(const char *scale_file = DefaultMBMbdScaleFile())
{
	static TString cached_scale_file = "";
	static double cached_scale = 1.;

	if (cached_scale_file != scale_file)
	{
		cached_scale_file = scale_file;
		cached_scale = 1.;

		TFile *fscale = TFile::Open(scale_file, "READ");
		if (!fscale || fscale->IsZombie())
		{
			std::cerr << "ERROR: cannot open MBD scale file " << scale_file << std::endl;
			if (fscale)
				fscale->Close();
			std::exit(1);
		}

		TParameter<double> *best_scale =
			dynamic_cast<TParameter<double> *>(fscale->Get("best_scale_value"));
		if (!best_scale)
		{
			std::cerr << "ERROR: missing best_scale_value in " << scale_file << std::endl;
			fscale->Close();
			std::exit(1);
		}

		cached_scale = best_scale->GetVal();
		fscale->Close();
	}

	return cached_scale;
}

auto FindCentBinData = [](float mbd_charge_sum) -> int
{
	for (int i = 0; i < nCentBins; i++)
	{
		if (mbd_charge_sum < mbdcharge_boundaries[i] && mbd_charge_sum >= mbdcharge_boundaries[i + 1])
			return i;
	}
	return -1;
};

auto FindCentBinMC = [](float mbd_charge_sum) -> int
{
	const double mbd_charge_sum_scaled = mbd_charge_sum * GetMbdChargeSumScaleMC();
	for (int i = 0; i < nCentBins; i++)
	{
		if (mbd_charge_sum_scaled < mbdcharge_boundaries_mc[i] && mbd_charge_sum_scaled >= mbdcharge_boundaries_mc[i + 1])
			return i;
	}
	return -1;
};

// for putting 2D response into 1D & back
int GlobalBin(int ix, int iy, int nBinsY)
{
	return (ix - 1) * nBinsY + iy;
}

void InverseGlobalBin(int globalBin, int nBinsY, int &ix, int &iy)
{
	int g = globalBin - 1; // 0-based
	ix = g / nBinsY + 1;
	iy = g % nBinsY + 1;
}

void getDijets(string infile = "/sphenix/tg/tg01/jets/jpark4/Run25OO/TTrees/SkimmedTrees/outree_skimjet_MC_inclusive_merged.root", string outfile = "output.root", bool ismc = true)
{
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	TH3::SetDefaultSumw2();

	TFile *f = new TFile(infile.c_str());
	TTree *t = (TTree *)f->Get("tree");

	Short_t njetsdata;
	Int_t njetsmc;
	Float_t eta[10];
	Float_t phi[10];
	Float_t ptcalib[10];
	Float_t e[10];
	Float_t pt[10];
	Float_t time[10];
	Float_t mbdtime;
	Float_t mbdcharge;

	Int_t ntruthjets;
	Float_t trutheta[10];
	Float_t truthphi[10];
	Float_t truthe[10];
	Float_t truthpt[10];
	Double_t weight = 1.;

	t->SetBranchAddress("mbd_mean_time", &mbdtime);
	t->SetBranchAddress("mbd_charge_sum", &mbdcharge);

	if (!ismc)
	{
		t->SetBranchAddress("nJetsSub_r04", &njetsdata);
		t->SetBranchAddress("jet_sub_r04_eta", &eta);
		t->SetBranchAddress("jet_sub_r04_phi", &phi);
		t->SetBranchAddress("jet_sub_r04_e", &e);
		t->SetBranchAddress("jet_sub_r04_pt", &pt);
		t->SetBranchAddress("jet_sub_r04_calib_pt", &ptcalib);
		t->SetBranchAddress("jet_sub_r04_time", &time);
	}

	if (ismc)
	{
		t->SetBranchAddress("nJetsReco_r04", &njetsmc);
		t->SetBranchAddress("jet_reco_r04_eta", &eta);
		t->SetBranchAddress("jet_reco_r04_phi", &phi);
		t->SetBranchAddress("jet_reco_r04_e", &e);
		t->SetBranchAddress("jet_reco_r04_pt", &pt);
		t->SetBranchAddress("jet_reco_r04_pt_calib", &ptcalib);

		t->SetBranchAddress("nJetsTruth_r04", &ntruthjets);
		t->SetBranchAddress("jet_truth_r04_eta", &trutheta);
		t->SetBranchAddress("jet_truth_r04_phi", &truthphi);
		t->SetBranchAddress("jet_truth_r04_pt", &truthpt);
		t->SetBranchAddress("weight", &weight);
	}

	// make hists

	// log pt binning
	/*const int pt_N = 19;
	Float_t pt_bins[pt_N+1];
	float minVal = 5.5;
	float maxVal = 68;
	float alpha = pow(maxVal/minVal,1./pt_N);
	for(int i = 0; i < pt_N+1; i++){
	  pt_bins[i] = minVal*pow(alpha,i);
	  std::cout<<"Bin "<<i<<" "<<pt_bins[i]<<std::endl;
	  }*/

	// Dan binning
	const int pt_N = 19;
	Float_t pt_bins[pt_N + 1];
	float minVal = 5.5;
	float fixedBin = 14.; // this is a fixed pt1 bin we want
	int beforebins = 7;	  // this is the number of bins below the fixed bin value
	float alpha = pow(fixedBin / minVal, 1. / beforebins);
	float maxVal = minVal * TMath::Power(alpha, pt_N);

	for (int i = 0; i < pt_N + 1; i++)
	{
		pt_bins[i] = minVal * pow(alpha, i);
		std::cout << "Bin " << i << " " << pt_bins[i] << std::endl;
	}

	// truth dijet / response matrix phase space, aligned to bin edges 9 and 3 (~18.2836, ~8.20847)
	const float leadpttrue_min = pt_bins[9];
	const float subpttrue_min = pt_bins[3];

	// cent binning
	int cent_N = 4;
	Float_t cent_bins[] = {0, 1, 2, 3, 4};

	// xj binning
	int xj_N = pt_N;
	Float_t xj_bins[xj_N + 1];

	for (int i = 0; i < xj_N + 1; i++)
	{
		xj_bins[i] = pow(alpha, i - xj_N);
		std::cout << "xj bin " << xj_bins[i] << std::endl;
	}

	// dphi binning
	int dphi_N = 20;
	Float_t dphi_bins[dphi_N + 1];
	for (int i = 0; i < dphi_N + 1; i++)
	{
		dphi_bins[i] = i * (TMath::Pi() / dphi_N);
	}

	TH3F *h_dphi = new TH3F("h_dphi", "", dphi_N, dphi_bins, pt_N, pt_bins, cent_N, cent_bins);
	TH3F *h_xj = new TH3F("h_xj", "", xj_N, xj_bins, pt_N, pt_bins, cent_N, cent_bins);
	TH3F *h_pt1pt2 = new TH3F("h_pt1pt2", "", pt_N, pt_bins, pt_N, pt_bins, cent_N, cent_bins);

	// only filled in MC
	TH2D *hTrue2D[cent_N];
	TH2D *hMeas2D[cent_N];
	TH1D *hTrue1D[cent_N];
	TH1D *hMeas1D[cent_N];
	TH2D *hTrue2DA[cent_N];
	TH2D *hMeas2DA[cent_N];
	TH1D *hTrue1DA[cent_N];
	TH1D *hMeas1DA[cent_N];
	TH2D *hTrue2DB[cent_N];
	TH2D *hMeas2DB[cent_N];
	TH1D *hTrue1DB[cent_N];
	TH1D *hMeas1DB[cent_N];
	TH1D *hFake1D[cent_N];
	TH1D *hFake1DA[cent_N];
	TH1D *hFake1DB[cent_N];
	TH1D *h_counting[cent_N];
	TH3F *h_pt1pt2fake = new TH3F("h_pt1pt2fake", "", pt_N, pt_bins, pt_N, pt_bins, cent_N, cent_bins);

	std::vector<RooUnfoldResponse> response;
	response.reserve(nCentBins);
	std::vector<RooUnfoldResponse> responseA;
	responseA.reserve(nCentBins);
	std::vector<RooUnfoldResponse> responseB;
	responseB.reserve(nCentBins);

	int nGlobalBins = pt_N * pt_N;
	for (int i = 0; i < cent_N; i++)
	{
		hTrue2D[i] = new TH2D(Form("hTrue2D%i", i), "Truth;x;y", pt_N, pt_bins, pt_N, pt_bins);
		hMeas2D[i] = new TH2D(Form("hMeas2D%i", i), "Measured;x;y", pt_N, pt_bins, pt_N, pt_bins);

		hTrue2DA[i] = new TH2D(Form("hTrue2DA%i", i), "Truth;x;y", pt_N, pt_bins, pt_N, pt_bins);
		hMeas2DA[i] = new TH2D(Form("hMeas2DA%i", i), "Measured;x;y", pt_N, pt_bins, pt_N, pt_bins);

		hTrue2DB[i] = new TH2D(Form("hTrue2DB%i", i), "Truth;x;y", pt_N, pt_bins, pt_N, pt_bins);
		hMeas2DB[i] = new TH2D(Form("hMeas2DB%i", i), "Measured;x;y", pt_N, pt_bins, pt_N, pt_bins);

		hTrue1D[i] = new TH1D(Form("hTrue1D%i", i), "Truth (flattened);global bin", nGlobalBins, 0.5, nGlobalBins + 0.5);
		hTrue1DA[i] = new TH1D(Form("hTrue1DA%i", i), "Truth (flattened);global bin", nGlobalBins, 0.5, nGlobalBins + 0.5);
		hTrue1DB[i] = new TH1D(Form("hTrue1DB%i", i), "Truth (flattened);global bin", nGlobalBins, 0.5, nGlobalBins + 0.5);

		hMeas1D[i] = new TH1D(Form("hMeas1D%i", i), "Measured (flattened);global bin", nGlobalBins, 0.5, nGlobalBins + 0.5);
		hMeas1DA[i] = new TH1D(Form("hMeas1DA%i", i), "Measured (flattened);global bin", nGlobalBins, 0.5, nGlobalBins + 0.5);
		hMeas1DB[i] = new TH1D(Form("hMeas1DB%i", i), "Measured (flattened);global bin", nGlobalBins, 0.5, nGlobalBins + 0.5);

		hFake1D[i] = new TH1D(Form("hFake1D%i", i), "Fakes (flattened);global bin", nGlobalBins, 0.5, nGlobalBins + 0.5);
		hFake1DA[i] = new TH1D(Form("hFake1DA%i", i), "Fakes (flattened);global bin", nGlobalBins, 0.5, nGlobalBins + 0.5);
		hFake1DB[i] = new TH1D(Form("hFake1DB%i", i), "Fakes (flattened);global bin", nGlobalBins, 0.5, nGlobalBins + 0.5);

		response.emplace_back(hMeas1D[i], hTrue1D[i], Form("response%i", i), "2D response (flattened)");
		responseA.emplace_back(hMeas1DA[i], hTrue1DA[i], Form("responseA%i", i), "2D response (flattened)");
		responseB.emplace_back(hMeas1DB[i], hTrue1DB[i], Form("responseB%i", i), "2D response (flattened)");

		h_counting[i] = new TH1D(Form("h_counting%i", i), "", 7, 0.5, 7.5); // total events, pass cent, reco dijet, truth dijet, match dijet, unmatch reco, unmatch truth
	}

	int nentries = t->GetEntries();
	for (int i = 0; i < nentries; i++)
	{
		t->GetEntry(i);

		int half = i % 2;

		// event level cuts
		int centbin = -999;
		if (!ismc)
			centbin = FindCentBinData(mbdcharge);
		else
			centbin = FindCentBinMC(mbdcharge);
		if (centbin < 0)
			continue;
		int njets = 0;
		if (ismc)
			njets = njetsmc;
		else
			njets = njetsdata;

		h_counting[centbin]->Fill(1);

		// loop through jets to get lead & sublead
		float leadpt = 0;
		float leadphi = 0;
		float leadtime = 999;
		float leadeta = 0;
		float subpt = 0;
		float subphi = 0;
		float subeta = 0;
		float subtime = 999;
		bool recodijet = false;

		if (njets >= 2)
		{
			for (int j = 0; j < njets; j++)
			{
				// jet level cuts

				if (ptcalib[j] > leadpt)
				{
					subpt = leadpt;
					subphi = leadphi;
					subeta = leadeta;
					subtime = leadtime;
					leadpt = ptcalib[j];
					leadphi = phi[j];
					leadeta = eta[j];
					leadtime = time[j];
				}
				else if (ptcalib[j] > subpt)
				{
					subpt = ptcalib[j];
					subphi = phi[j];
					subeta = eta[j];
					subtime = time[j];
				}
			}

			// jet time cuts for data only
			if (!(!ismc && (fabs(leadtime) > 6 || fabs(leadtime - subtime) > 3)))
			{

				// jet pt cuts
				if (subpt >= 12)
				{

					// check if back to back
					float dPhi = leadphi - subphi;
					while (dPhi > TMath::Pi())
						dPhi -= 2 * TMath::Pi();
					while (dPhi < -TMath::Pi())
						dPhi += 2 * TMath::Pi();
					h_dphi->Fill(dPhi, leadpt, centbin, weight);
					if (fabs(dPhi) > 3 * TMath::Pi() / 4) // backtoback
					{
						h_xj->Fill(subpt / leadpt, leadpt, centbin, weight);
						h_pt1pt2->Fill(subpt, leadpt, centbin, 0.5 * weight);
						h_pt1pt2->Fill(leadpt, subpt, centbin, 0.5 * weight);
						int ixReco = hMeas2D[centbin]->GetXaxis()->FindBin(leadpt);
						int iyReco = hMeas2D[centbin]->GetYaxis()->FindBin(subpt);
						int gReco = GlobalBin(ixReco, iyReco, pt_N);
						int gRecosym = GlobalBin(iyReco, ixReco, pt_N);

						if (!ismc) // fill data here, otherwise fill later in mc bit
						{
							hMeas2D[centbin]->Fill(leadpt, subpt, weight * 0.5);
							hMeas2D[centbin]->Fill(subpt, leadpt, weight * 0.5);
							hMeas1D[centbin]->Fill(gReco, weight * 0.5);
							hMeas1D[centbin]->Fill(gRecosym, weight * 0.5);
						}
						recodijet = true;
						h_counting[centbin]->Fill(2);
					}
				}
			}
		}

		// if MC then also run this for truth jets and see if those are matched to the reco jets
		if (ismc)
		{
			bool ismatch = false;
			bool swapmatch = false;
			bool truthdijet = false;
			float leadpttrue = 0;
			float leadphitrue = 0;
			float leadetatrue = 0;
			float subpttrue = 0;
			float subphitrue = 0;
			float subetatrue = 0;
			if (ntruthjets >= 2)
			{
				for (int j = 0; j < ntruthjets; j++)
				{
					if (truthpt[j] > leadpttrue)
					{
						subpttrue = leadpttrue;
						subphitrue = leadphitrue;
						subetatrue = leadetatrue;
						leadpttrue = truthpt[j];
						leadphitrue = truthphi[j];
						leadetatrue = trutheta[j];
					}
					else if (truthpt[j] > subpttrue)
					{
						subpttrue = truthpt[j];
						subphitrue = truthphi[j];
						subetatrue = trutheta[j];
					}
				}
				// truth jet cuts
				if (subpttrue >= subpttrue_min && leadpttrue >= leadpttrue_min && fabs(leadetatrue) < 0.7 && fabs(subetatrue) < 0.7)
				{
					// check if truth is back to back
					float dPhi = leadphi - subphi;
					while (dPhi > TMath::Pi())
						dPhi -= 2 * TMath::Pi();
					while (dPhi < -TMath::Pi())
						dPhi += 2 * TMath::Pi();

					if (fabs(dPhi) > 3 * TMath::Pi() / 4) // backtoback
					{
						// truth dijet
						truthdijet = true;
						h_counting[centbin]->Fill(3);

						// then check if matched to reco
						if (recodijet) // only match if the reco passes dijet cuts
						{
							// direct matching: lead reco - lead truth, sublead reco - sublead truth
							dPhi = leadphi - leadphitrue;
							while (dPhi > TMath::Pi())
								dPhi -= 2 * TMath::Pi();
							while (dPhi < -TMath::Pi())
								dPhi += 2 * TMath::Pi();
							float deta = leadeta - leadetatrue;
							float dR = sqrt(dPhi * dPhi + deta * deta);
							bool leadmatchlead = (dR < 0.3);

							dPhi = subphi - subphitrue;
							while (dPhi > TMath::Pi())
								dPhi -= 2 * TMath::Pi();
							while (dPhi < -TMath::Pi())
								dPhi += 2 * TMath::Pi();
							deta = subeta - subetatrue;
							dR = sqrt(dPhi * dPhi + deta * deta);
							bool submatchsub = (dR < 0.3);

							// swapped matching: lead reco - sublead truth, sublead reco - lead truth
							dPhi = leadphi - subphitrue;
							while (dPhi > TMath::Pi())
								dPhi -= 2 * TMath::Pi();
							while (dPhi < -TMath::Pi())
								dPhi += 2 * TMath::Pi();
							deta = leadeta - subetatrue;
							dR = sqrt(dPhi * dPhi + deta * deta);
							bool leadmatchsub = (dR < 0.3);

							dPhi = subphi - leadphitrue;
							while (dPhi > TMath::Pi())
								dPhi -= 2 * TMath::Pi();
							while (dPhi < -TMath::Pi())
								dPhi += 2 * TMath::Pi();
							deta = subeta - leadetatrue;
							dR = sqrt(dPhi * dPhi + deta * deta);
							bool submatchlead = (dR < 0.3);

							if (leadmatchlead && submatchsub)
							{
								ismatch = true; // lead reco matches lead truth, sublead reco matches sublead truth
								h_counting[centbin]->Fill(4);
							}
							else if (leadmatchsub && submatchlead)
							{
								ismatch = true; // lead reco matches sublead truth, sublead reco matches lead truth
								swapmatch = true;
								h_counting[centbin]->Fill(4);
							}
						}
					}
				}
			}
			// now depending on if we have reco/truth/matched dijets fill response

			// get global binning
			int ixReco = hMeas2D[centbin]->GetXaxis()->FindBin(leadpt);
			int iyReco = hMeas2D[centbin]->GetYaxis()->FindBin(subpt);
			int gReco = GlobalBin(ixReco, iyReco, pt_N);
			int gRecosym = GlobalBin(iyReco, ixReco, pt_N);

			int ixTrue = hTrue2D[centbin]->GetXaxis()->FindBin(leadpttrue);
			int iyTrue = hTrue2D[centbin]->GetYaxis()->FindBin(subpttrue);
			int gTrue = GlobalBin(ixTrue, iyTrue, pt_N);
			int gTruesym = GlobalBin(iyTrue, ixTrue, pt_N);

			// std::cout<<"is reco "<<recodijet<<" is truth "<<truthdijet<<" Match? "<<ismatch<<std::endl;
			if (ismatch)
			{
				hMeas2D[centbin]->Fill(leadpt, subpt, weight * 0.5);
				hMeas2D[centbin]->Fill(subpt, leadpt, weight * 0.5); // symmetrize
				hTrue2D[centbin]->Fill(leadpttrue, subpttrue, weight * 0.5);
				hTrue2D[centbin]->Fill(subpttrue, leadpttrue, weight * 0.5);

				hMeas1D[centbin]->Fill(gReco, weight * 0.5);
				hMeas1D[centbin]->Fill(gRecosym, weight * 0.5);
				hTrue1D[centbin]->Fill(gTrue, weight * 0.5);
				hTrue1D[centbin]->Fill(gTruesym, weight * 0.5);
				// if swapmatch, the reco lead/sublead pt bin actually corresponds to the truth sublead/lead pt bin
				if (!swapmatch)
				{
					response[centbin].Fill(gReco, gTrue, weight * 0.5);
					response[centbin].Fill(gRecosym, gTruesym, weight * 0.5);
				}
				else
				{
					response[centbin].Fill(gReco, gTruesym, weight * 0.5);
					response[centbin].Fill(gRecosym, gTrue, weight * 0.5);
				}

				if (half == 0)
				{
					hMeas2DA[centbin]->Fill(leadpt, subpt, weight * 0.5);
					hMeas2DA[centbin]->Fill(subpt, leadpt, weight * 0.5);

					hTrue2DA[centbin]->Fill(leadpttrue, subpttrue, weight * 0.5);
					hTrue2DA[centbin]->Fill(subpttrue, leadpttrue, weight * 0.5);

					hMeas1DA[centbin]->Fill(gReco, weight * 0.5);
					hMeas1DA[centbin]->Fill(gRecosym, weight * 0.5);
					hTrue1DA[centbin]->Fill(gTrue, weight * 0.5);
					hTrue1DA[centbin]->Fill(gTruesym, weight * 0.5);
					if (!swapmatch)
					{
						responseA[centbin].Fill(gReco, gTrue, weight * 0.5);
						responseA[centbin].Fill(gRecosym, gTruesym, weight * 0.5);
					}
					else
					{
						responseA[centbin].Fill(gReco, gTruesym, weight * 0.5);
						responseA[centbin].Fill(gRecosym, gTrue, weight * 0.5);
					}
				}
				else
				{
					hMeas2DB[centbin]->Fill(leadpt, subpt, weight * 0.5);
					hMeas2DB[centbin]->Fill(subpt, leadpt, weight * 0.5);

					hTrue2DB[centbin]->Fill(leadpttrue, subpttrue, weight * 0.5);
					hTrue2DB[centbin]->Fill(subpttrue, leadpttrue, weight * 0.5);

					hMeas1DB[centbin]->Fill(gReco, weight * 0.5);
					hMeas1DB[centbin]->Fill(gRecosym, weight * 0.5);
					hTrue1DB[centbin]->Fill(gTrue, weight * 0.5);
					hTrue1DB[centbin]->Fill(gTruesym, weight * 0.5);
					if (!swapmatch)
					{
						responseB[centbin].Fill(gReco, gTrue, weight * 0.5);
						responseB[centbin].Fill(gRecosym, gTruesym, weight * 0.5);
					}
					else
					{
						responseB[centbin].Fill(gReco, gTruesym, weight * 0.5);
						responseB[centbin].Fill(gRecosym, gTrue, weight * 0.5);
					}
				}

				// std::cout<<"Filling response"<<std::endl;
			}

			else if (truthdijet) // umatched truth dijet is miss
			{
				hTrue2D[centbin]->Fill(leadpttrue, subpttrue, weight * 0.5);
				hTrue2D[centbin]->Fill(subpttrue, leadpttrue, weight * 0.5);
				hTrue1D[centbin]->Fill(gTrue, weight * 0.5);
				hTrue1D[centbin]->Fill(gTruesym, weight * 0.5);
				response[centbin].Miss(gTrue, weight * 0.5);
				response[centbin].Miss(gTruesym, weight * 0.5);
				// std::cout<<"Filling miss"<<std::endl;
				h_counting[centbin]->Fill(6);

				if (half == 0)
				{
					hTrue2DA[centbin]->Fill(leadpttrue, subpttrue, weight * 0.5);
					hTrue2DA[centbin]->Fill(subpttrue, leadpttrue, weight * 0.5);
					hTrue1DA[centbin]->Fill(gTrue, weight * 0.5);
					hTrue1DA[centbin]->Fill(gTruesym, weight * 0.5);
					responseA[centbin].Miss(gTrue, weight * 0.5);
					responseA[centbin].Miss(gTruesym, weight * 0.5);
				}
				else
				{
					hTrue2DB[centbin]->Fill(leadpttrue, subpttrue, weight * 0.5);
					hTrue2DB[centbin]->Fill(subpttrue, leadpttrue, weight * 0.5);
					hTrue1DB[centbin]->Fill(gTrue, weight * 0.5);
					hTrue1DB[centbin]->Fill(gTruesym, weight * 0.5);
					responseB[centbin].Miss(gTrue, weight * 0.5);
					responseB[centbin].Miss(gTruesym, weight * 0.5);
				}

				if (recodijet) // if there's also a reco dijet and they are just unmatched this is also a fake
				{
					hMeas2D[centbin]->Fill(leadpt, subpt, weight * 0.5);
					hMeas2D[centbin]->Fill(subpt, leadpt, weight * 0.5);
					hMeas1D[centbin]->Fill(gReco, weight * 0.5);
					hMeas1D[centbin]->Fill(gRecosym, weight * 0.5);
					// this just doesn't work in our version of roounfold >:(
					// response[centbin].Fake(gReco,weight*0.5);
					// response[centbin].Fake(gRecosym,weight*0.5);
					hFake1D[centbin]->Fill(gReco, weight * 0.5);
					hFake1D[centbin]->Fill(gRecosym, weight * 0.5);

					h_pt1pt2fake->Fill(leadpt, subpt, centbin, weight * 0.5);
					h_pt1pt2fake->Fill(subpt, leadpt, centbin, weight * 0.5);
					// std::cout<<"Filling fake & miss"<<std::endl;
					h_counting[centbin]->Fill(5);

					if (half == 0)
					{
						hMeas2DA[centbin]->Fill(leadpt, subpt, weight * 0.5);
						hMeas2DA[centbin]->Fill(subpt, leadpt, weight * 0.5);
						hMeas1DA[centbin]->Fill(gReco, weight * 0.5);
						hMeas1DA[centbin]->Fill(gRecosym, weight * 0.5);
						// responseA[centbin].Fake(gReco,weight*0.5);
						// responseA[centbin].Fake(gRecosym,weight*0.5);
						hFake1DA[centbin]->Fill(gReco, weight * 0.5);
						hFake1DA[centbin]->Fill(gRecosym, weight * 0.5);
					}
					else
					{
						hMeas2DB[centbin]->Fill(leadpt, subpt, weight * 0.5);
						hMeas2DB[centbin]->Fill(subpt, leadpt, weight * 0.5);
						hMeas1DB[centbin]->Fill(gReco, weight * 0.5);
						hMeas1DB[centbin]->Fill(gRecosym, weight * 0.5);
						// responseB[centbin].Fake(gReco,weight*0.5);
						// responseB[centbin].Fake(gRecosym,weight*0.5);
						hFake1DB[centbin]->Fill(gReco, weight * 0.5);
						hFake1DB[centbin]->Fill(gRecosym, weight * 0.5);
					}
				}
			}
			else if (recodijet) // no truth dijet but there is a reco, this is a fake
			{
				hMeas2D[centbin]->Fill(leadpt, subpt, weight * 0.5);
				hMeas2D[centbin]->Fill(subpt, leadpt, weight * 0.5);
				hMeas1D[centbin]->Fill(gReco, weight * 0.5);
				hMeas1D[centbin]->Fill(gRecosym, weight * 0.5);
				// response[centbin].Fake(gReco,weight*0.5);
				// response[centbin].Fake(gRecosym,weight*0.5);
				hFake1D[centbin]->Fill(gReco, weight * 0.5);
				hFake1D[centbin]->Fill(gRecosym, weight * 0.5);
				h_pt1pt2fake->Fill(leadpt, subpt, centbin, weight * 0.5);
				h_pt1pt2fake->Fill(subpt, leadpt, centbin, weight * 0.5);
				h_counting[centbin]->Fill(5);

				if (half == 0)
				{
					hMeas2DA[centbin]->Fill(leadpt, subpt, weight * 0.5);
					hMeas2DA[centbin]->Fill(subpt, leadpt, weight * 0.5);
					hMeas1DA[centbin]->Fill(gReco, weight * 0.5);
					hMeas1DA[centbin]->Fill(gRecosym, weight * 0.5);
					hFake1DA[centbin]->Fill(gReco, weight * 0.5);
					hFake1DA[centbin]->Fill(gRecosym, weight * 0.5);
					// responseA[centbin].Fake(gReco,weight*0.5);
					// responseA[centbin].Fake(gRecosym,weight*0.5);
				}
				else
				{
					hMeas2DB[centbin]->Fill(leadpt, subpt, weight * 0.5);
					hMeas2DB[centbin]->Fill(subpt, leadpt, weight * 0.5);
					hMeas1DB[centbin]->Fill(gReco, weight * 0.5);
					hMeas1DB[centbin]->Fill(gRecosym, weight * 0.5);
					hFake1DB[centbin]->Fill(gReco, weight * 0.5);
					hFake1DB[centbin]->Fill(gRecosym, weight * 0.5);
					// responseB[centbin].Fake(gReco,weight*0.5);
					// responseB[centbin].Fake(gRecosym,weight*0.5);
				}
			}
		}
	}

	TFile *fout = new TFile(outfile.c_str(), "RECREATE");
	h_dphi->Write();
	h_xj->Write();
	h_pt1pt2->Write();

	for (int i = 0; i < cent_N; i++)
	{
		hTrue2D[i]->Write();
		hMeas2D[i]->Write();

		hTrue1D[i]->Write();
		hMeas1D[i]->Write();
		hFake1D[i]->Write();

		response[i].Write();
		h_counting[i]->Write();

		hTrue2DA[i]->Write();
		hMeas2DA[i]->Write();

		hTrue1DA[i]->Write();
		hMeas1DA[i]->Write();
		hFake1DA[i]->Write();

		responseA[i].Write();

		hTrue2DB[i]->Write();
		hMeas2DB[i]->Write();

		hTrue1DB[i]->Write();
		hMeas1DB[i]->Write();
		hFake1DB[i]->Write();

		responseB[i].Write();
	}

	h_pt1pt2fake->Write();
}
