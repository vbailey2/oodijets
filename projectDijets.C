void projectDijets()
{
	TFile *f = new TFile("hists/hist-full.root", "READ");
	TH3F *h_pt1pt2 = (TH3F *)f->Get("h_pt1pt2");
	int ncent = h_pt1pt2->GetNbinsZ();
	TH2F *h_2D[ncent];
	const int npt = h_pt1pt2->GetNbinsX();
	const int nfinal = 3;				  // number of pt bins in final measurement
	int finalptbins[] = {11, 13, 15, 18}; // pt bin indices of those
										  // corresponds to bin low edges of 23.9,31.2,40.7,60.8)
	TH2F *h_2D_unsym[ncent];
	TH2F *h_xj[ncent];
	TH3F *h_xj_bins = (TH3F *)f->Get("h_xj");
	const TArrayD *xjaxis = h_xj_bins->GetXaxis()->GetXbins();
	Double_t final_pt_binning[nfinal + 1];
	for (int i = 0; i < nfinal + 1; i++)
	{
		final_pt_binning[i] = h_xj_bins->GetYaxis()->GetBinLowEdge(finalptbins[i] + 1);
		// std::cout<<"bin "<<i <<" is "<<final_pt_binning[i]<<std::endl;
	}

	for (int i = 0; i < ncent; i++)
	{
		h_pt1pt2->GetZaxis()->SetRange(i + 1, i + 1);
		h_2D[i] = (TH2F *)h_pt1pt2->Project3D("yx");
		h_2D[i]->SetName(Form("h_2D_%i", i));

		h_2D_unsym[i] = (TH2F *)h_2D[i]->Clone();
		h_2D_unsym[i]->Reset();
		h_2D_unsym[i]->SetName(Form("h_2D_unsym_%i", i));

		// unsymmetrize
		for (int ipt2 = 0; ipt2 < npt; ipt2++)
		{
			for (int ipt1 = ipt2; ipt1 < npt; ipt1++)
			{
				// std::cout<<"doing bin "<<ipt1<<" , "<<ipt2<<" and "<<ipt2<<" , "<<ipt1<<std::endl;
				float bincount = 0;
				bincount += h_2D[i]->GetBinContent(ipt1 + 1, ipt2 + 1);
				if (ipt1 != ipt2)
					bincount += h_2D[i]->GetBinContent(ipt2 + 1, ipt1 + 1);
				float binerr = 0;
				binerr += h_2D[i]->GetBinError(ipt1 + 1, ipt2 + 1) * h_2D[i]->GetBinError(ipt1 + 1, ipt2 + 1);
				if (ipt1 != ipt2)
					binerr += h_2D[i]->GetBinError(ipt2 + 1, ipt1 + 1) * h_2D[i]->GetBinError(ipt2 + 1, ipt1 + 1);
				h_2D_unsym[i]->SetBinContent(ipt1 + 1, ipt2 + 1, bincount);
				h_2D_unsym[i]->SetBinError(ipt1 + 1, ipt2 + 1, std::sqrt(binerr));
			}
		}

		// project
		// h_xj[i] = (TH2F*)h_xj_bins->Project3D("yx");
		h_xj[i] = new TH2F(Form("h_xj_%i", i), "", h_xj_bins->GetNbinsX(), xjaxis->GetArray(), nfinal, final_pt_binning);
		h_xj[i]->SetName(Form("h_xj_%i", i));
		float binvals[nfinal][19] = {};
		float binerrs[nfinal][19] = {};
		for (int ipt2 = 0; ipt2 < npt; ipt2++)
		{
			for (int ipt1 = ipt2; ipt1 < npt; ipt1++)
			{
				// map this bin to one of the final bins
				int finalbin = -1;
				for (int ibin = 0; ibin < nfinal; ibin++)
				{
					if (ipt1 >= finalptbins[ibin] && ipt1 < finalptbins[ibin + 1])
						finalbin = ibin;
				}
				// std::cout<<"This is bin "<<ipt1<<" which maps to final bin "<<finalbin<<std::endl;
				if (finalbin < 0)
					continue;
				double binval = h_2D_unsym[i]->GetBinContent(ipt1 + 1, ipt2 + 1);
				double binerr = h_2D_unsym[i]->GetBinError(ipt1 + 1, ipt2 + 1);
				// std::cout<<binval<<std::endl;
				if (ipt1 == ipt2)
				{
					binvals[finalbin][npt - 1] += binval;
					binerrs[finalbin][npt - 1] += binerr * binerr;
				}
				else
				{
					int diff = ipt1 - ipt2;
					int bin = npt - diff;
					// std::cout<<"pt "<<ipt1<<" "<<ipt2<<" is bin "<<bin<<std::endl;
					binvals[finalbin][bin] += binval / 2.;
					binvals[finalbin][bin - 1] += binval / 2.;
					binerrs[finalbin][bin] += (binerr / 2.) * (binerr / 2.);
					binerrs[finalbin][bin - 1] += (binerr / 2.) * (binerr / 2.);
				}
			}
		}

		// fill xj hist
		for (int ix = 0; ix < npt; ix++)
		{
			for (int ipt = 0; ipt < nfinal; ipt++)
			{
				h_xj[i]->SetBinContent(ix + 1, ipt + 1, binvals[ipt][ix]);
				h_xj[i]->SetBinError(ix + 1, ipt + 1, std::sqrt(binerrs[ipt][ix]));
			}
		}

	} // cent loop

	TFile *fout = new TFile("hists/projections.root", "RECREATE");
	for (int i = 0; i < ncent; i++)
	{
		h_2D[i]->Write();
		h_2D_unsym[i]->Write();
		h_xj[i]->Write();
	}
	std::cout << "all done" << std::endl;
}
