void projectUnfoldClosure(bool ishalf = 0)
{
  //TFile *f = new TFile("hists/hist-unfoldedMC.root","READ");
  TFile *f = new TFile("hists/hist-unfoldedData.root","READ");
  TH3F *h_pt1pt2 = (TH3F*)f->Get("h_pt1pt2");
  int ncent = h_pt1pt2->GetNbinsZ();
  TH2F *h_2Dunfold[ncent];
  TH2F *h_2Dtrue[ncent];
  TH2F *h_2Dmeas[ncent];
  const int npt = h_pt1pt2->GetNbinsX();
  const int nfinal = 3; //number of pt bins in final measurement
  int finalptbins[] = {11,13,15,18}; //pt bin indices of those
                                     //corresponds to bin low edges of 23.9,31.2,40.7,60.8)
  TH2F *h_2D_unsymunfold[ncent];
  TH2F *h_2D_unsymmeas[ncent];
  TH2F *h_2D_unsymtrue[ncent];
  TH2F *h_xjunfold[ncent];
  TH2F *h_xjtrue[ncent];
  TH2F *h_xjmeas[ncent];
  TH3F *h_xj_bins = (TH3F*)f->Get("h_xj");
  const TArrayD *xjaxis = h_xj_bins->GetXaxis()->GetXbins();
  Double_t final_pt_binning[nfinal+1];
  for(int i = 0; i < nfinal+1; i++)
    {
      final_pt_binning[i] = h_xj_bins->GetYaxis()->GetBinLowEdge(finalptbins[i]+1);
      //std::cout<<"bin "<<i <<" is "<<final_pt_binning[i]<<std::endl;
    }
  
  
  for(int i = 0; i < ncent; i++)
    {
      if(ishalf) h_2Dunfold[i] =(TH2F*) f->Get(Form("hUnfolded2DHalf_cent%i",i));
      else h_2Dunfold[i] =(TH2F*) f->Get(Form("hUnfolded2D_cent%i",i));
      h_2D_unsymunfold[i] = (TH2F*)h_2Dunfold[i]->Clone();
      h_2D_unsymunfold[i]->SetName(Form("h_2D_unsymunfold_%i",i));
      h_2D_unsymunfold[i]->Reset();
 
      if(ishalf)h_2Dtrue[i] =(TH2F*) f->Get(Form("hTrue2DB%i",i));
      else h_2Dtrue[i] =(TH2F*) f->Get(Form("hTrue2D%i",i));
      h_2D_unsymtrue[i] = (TH2F*)h_2Dtrue[i]->Clone();
      h_2D_unsymtrue[i]->SetName(Form("h_2D_unsymtrue_%i",i));
      h_2D_unsymtrue[i]->Reset();
      
      if(ishalf)h_2Dmeas[i] =(TH2F*) f->Get(Form("hMeas2DB%i",i));
      else h_2Dmeas[i] =(TH2F*) f->Get(Form("hMeas2D%i",i));
      h_2D_unsymmeas[i] = (TH2F*)h_2Dmeas[i]->Clone();
      h_2D_unsymmeas[i]->SetName(Form("h_2D_unsymmeas_%i",i));
      h_2D_unsymmeas[i]->Reset();
      
      //unsymmetrize 
      for(int ipt2 = 0; ipt2 < npt; ipt2++)
	{
	  for(int ipt1 = ipt2; ipt1 < npt; ipt1++)
	    {
	      //unfold
	      float bincount = 0;
	      bincount += h_2Dunfold[i]->GetBinContent(ipt1+1,ipt2+1);
	      bincount += h_2Dunfold[i]->GetBinContent(ipt2+1,ipt1+1);
	      float binerr = 0;
	      binerr += h_2Dunfold[i]->GetBinError(ipt1+1,ipt2+1)*h_2Dunfold[i]->GetBinError(ipt1+1,ipt2+1);
	      binerr += h_2Dunfold[i]->GetBinError(ipt2+1,ipt1+1)*h_2Dunfold[i]->GetBinError(ipt2+1,ipt1+1);
	      h_2D_unsymunfold[i]->SetBinContent(ipt1+1,ipt2+1, bincount);
	      h_2D_unsymunfold[i]->SetBinError(ipt1+1,ipt2+1, std::sqrt(binerr));

	      //measured
	      float bincountmeas = 0;
              bincountmeas += h_2Dmeas[i]->GetBinContent(ipt1+1,ipt2+1);
              bincountmeas += h_2Dmeas[i]->GetBinContent(ipt2+1,ipt1+1);
              float binerrmeas = 0;
              binerrmeas += h_2Dmeas[i]->GetBinError(ipt1+1,ipt2+1)*h_2Dmeas[i]->GetBinError(ipt1+1,ipt2+1);
              binerrmeas += h_2Dmeas[i]->GetBinError(ipt2+1,ipt1+1)*h_2Dmeas[i]->GetBinError(ipt2+1,ipt1+1);
              h_2D_unsymmeas[i]->SetBinContent(ipt1+1,ipt2+1, bincountmeas);
              h_2D_unsymmeas[i]->SetBinError(ipt1+1,ipt2+1, std::sqrt(binerrmeas));

	      //true
	      float bincounttrue = 0;
              bincounttrue += h_2Dtrue[i]->GetBinContent(ipt1+1,ipt2+1);
              bincounttrue += h_2Dtrue[i]->GetBinContent(ipt2+1,ipt1+1);
              float binerrtrue = 0;
              binerrtrue += h_2Dtrue[i]->GetBinError(ipt1+1,ipt2+1)*h_2Dtrue[i]->GetBinError(ipt1+1,ipt2+1);
              binerrtrue += h_2Dtrue[i]->GetBinError(ipt2+1,ipt1+1)*h_2Dtrue[i]->GetBinError(ipt2+1,ipt1+1);
              h_2D_unsymtrue[i]->SetBinContent(ipt1+1,ipt2+1, bincounttrue);
              h_2D_unsymtrue[i]->SetBinError(ipt1+1,ipt2+1, std::sqrt(binerrtrue));
	    }
	}

      //project
      h_xjunfold[i] = new TH2F(Form("h_xjunfold_%i",i),"",h_xj_bins->GetNbinsX(), xjaxis->GetArray(),nfinal,final_pt_binning);
      h_xjmeas[i] = new TH2F(Form("h_xjmeas_%i",i),"",h_xj_bins->GetNbinsX(), xjaxis->GetArray(),nfinal,final_pt_binning);
      h_xjtrue[i] = new TH2F(Form("h_xjtrue_%i",i),"",h_xj_bins->GetNbinsX(), xjaxis->GetArray(),nfinal,final_pt_binning);
      
      float binvals[nfinal][19] = {};
      float binerrs[nfinal][19] = {};
      float binvalsmeas[nfinal][19] = {};
      float binerrsmeas[nfinal][19] = {};
      float binvalstrue[nfinal][19] = {};
      float binerrstrue[nfinal][19] = {};
      for(int ipt2 = 0; ipt2 < npt; ipt2++)
        {
          for(int ipt1 = ipt2; ipt1 < npt; ipt1++)
            {
	      //map this bin to one of the final bins
	      int finalbin = -1;
	      for(int ibin = 0; ibin < nfinal; ibin++)
		{
		  if(ipt1 >= finalptbins[ibin] && ipt1 < finalptbins[ibin+1]) finalbin = ibin;
		}

	      if(finalbin < 0) continue;
	      double binval = h_2D_unsymunfold[i]->GetBinContent(ipt1+1,ipt2+1);
	      double binerr = h_2D_unsymunfold[i]->GetBinError(ipt1+1,ipt2+1);
	      double binvalmeas = h_2D_unsymmeas[i]->GetBinContent(ipt1+1,ipt2+1);
              double binerrmeas = h_2D_unsymmeas[i]->GetBinError(ipt1+1,ipt2+1);
	      double binvaltrue = h_2D_unsymtrue[i]->GetBinContent(ipt1+1,ipt2+1);
              double binerrtrue = h_2D_unsymtrue[i]->GetBinError(ipt1+1,ipt2+1);
	      
	      if(ipt1 == ipt2)
		{
		  binvals[finalbin][npt-1] += binval;
		  binerrs[finalbin][npt-1] += binerr*binerr;
		  binvalsmeas[finalbin][npt-1] += binvalmeas;
                  binerrsmeas[finalbin][npt-1] += binerrmeas*binerrmeas;
		  binvalstrue[finalbin][npt-1] += binvaltrue;
                  binerrstrue[finalbin][npt-1] += binerrtrue*binerrtrue;
		}
	      else
		{
		  int diff = ipt1-ipt2;
		  int bin = npt-diff;
		  binvals[finalbin][bin] += binval/2.;
		  binvals[finalbin][bin-1] += binval/2.;
		  binerrs[finalbin][bin] += (binerr/2.)*(binerr/2.);
		  binerrs[finalbin][bin-1] += (binerr/2.)*(binerr/2.);

		  binvalsmeas[finalbin][bin] += binvalmeas/2.;
                  binvalsmeas[finalbin][bin-1] += binvalmeas/2.;
                  binerrsmeas[finalbin][bin] += (binerrmeas/2.)*(binerrmeas/2.);
                  binerrsmeas[finalbin][bin-1] += (binerrmeas/2.)*(binerrmeas/2.);

		  binvalstrue[finalbin][bin] += binvaltrue/2.;
                  binvalstrue[finalbin][bin-1] += binvaltrue/2.;
                  binerrstrue[finalbin][bin] += (binerrtrue/2.)*(binerrtrue/2.);
                  binerrstrue[finalbin][bin-1] += (binerrtrue/2.)*(binerrtrue/2.);
		}
	    }
	}

      //fill xj hist
      for(int ix = 0; ix < npt; ix++)
	{
	  for(int ipt = 0; ipt < nfinal; ipt++)
	    {
	      h_xjunfold[i]->SetBinContent(ix+1,ipt+1,binvals[ipt][ix]);	      
	      h_xjunfold[i]->SetBinError(ix+1,ipt+1,std::sqrt(binerrs[ipt][ix]));
	      h_xjmeas[i]->SetBinContent(ix+1,ipt+1,binvalsmeas[ipt][ix]);
              h_xjmeas[i]->SetBinError(ix+1,ipt+1,std::sqrt(binerrsmeas[ipt][ix]));
	      h_xjtrue[i]->SetBinContent(ix+1,ipt+1,binvalstrue[ipt][ix]);
              h_xjtrue[i]->SetBinError(ix+1,ipt+1,std::sqrt(binerrstrue[ipt][ix]));
	    }
	}
      
    }//cent loop

  string outname;
  if(ishalf) outname = "hists/projections_unfoldMC_Half.root";
  else outname = "hists/projections_unfoldData.root";
  TFile *fout = new TFile(outname.c_str(),"RECREATE");
  for(int i = 0; i < ncent; i++)
    {
      h_2Dunfold[i]->Write();
      h_2D_unsymunfold[i]->Write();
      h_xjunfold[i]->Write();
      h_xjmeas[i]->Write();
      h_xjtrue[i]->Write();
    }
  std::cout<<"all done"<<std::endl;
}
