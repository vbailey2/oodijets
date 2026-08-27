#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

void draw_dataMC()
{

  SetsPhenixStyle();
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();

  TFile *fdata = new TFile("hists/projections.root","READ");
  TFile *fMC = new TFile("hists/projectionsMC.root","READ");
  
  TCanvas *c = new TCanvas("c","c",700,700);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.01);   // small gap to bottom pad
  pad1->SetLeftMargin(0.12);
  pad1->Draw();
  
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
  pad2->SetTopMargin(0.03);
  pad2->SetBottomMargin(0.35);
  pad2->SetLeftMargin(0.12);
  pad2->Draw();
  
  TLegend *leg = new TLegend(.15,.75,.4,.92);
  leg->SetFillStyle(0);
  leg->AddEntry("","#it{#bf{sPHENIX}} Internal","");
  leg->AddEntry("","O+O #sqrt{s_{NN}} = 200 GeV","");
  leg->AddEntry("","anti-#it{k}_{#it{t}} #it{R} = 0.4, |#eta| < 0.7","");


  TLegend *cleg = new TLegend(.15,.5,.4,.75);
  cleg->SetFillStyle(0);
  
  TLegend *hleg = new TLegend(.6,.75,.9,.92);
  hleg->SetFillStyle(0);


  int colors[] = {1,2,4,kGreen+2, kViolet,kCyan,kOrange+2,kMagenta+2,kAzure-2};
  
 
  std::string cent_str[] = {"0-20%","20-40%","40-60%","60-80%"};
  int ncent = 4;
  
  TH2F *h_xjMC[ncent];
  TH2F *h_xjdata[ncent];
 
  h_xjMC[0] = (TH2F*)fMC->Get("h_xj_0");
  const int n_xj = h_xjMC[0]->GetNbinsX();
  const int npt = h_xjMC[0]->GetNbinsY();
  for(int i = 0; i < ncent; i++)
    {
      h_xjMC[i] = (TH2F*)fMC->Get(Form("h_xj_%i",i));
      h_xjMC[i]->SetName(Form("h_xjMC_%i",i));
      h_xjdata[i] = (TH2F*)fdata->Get(Form("h_xj_%i",i));
      h_xjdata[i]->SetName(Form("h_xjdata_%i",i));
      
    }
  
  TH1F *h_xj1DMC[ncent][npt];
  TH1F *h_xj1Ddata[ncent][npt];
  
  for(int ipt = 0; ipt < npt; ipt++)
    {

      for(int icent = 0; icent < ncent; icent++)
	{
	  pad1->cd();
	  h_xjMC[icent]->GetYaxis()->SetRange(ipt+1,ipt+1);	  
	  h_xj1DMC[icent][ipt] = (TH1F*)h_xjMC[icent]->ProjectionX();
	  h_xj1DMC[icent][ipt]->SetName(Form("h_xj_cent%i_pt%i_MC",icent,ipt));
	  h_xj1DMC[icent][ipt]->Scale(1./h_xj1DMC[icent][ipt]->Integral(),"width");
	  
	  h_xj1DMC[icent][ipt]->SetMarkerColor(colors[icent]);
	  h_xj1DMC[icent][ipt]->SetMarkerStyle(25);
	  h_xj1DMC[icent][ipt]->SetLineColor(colors[icent]);
	  h_xj1DMC[icent][ipt]->GetYaxis()->SetRangeUser(0,5);
	  h_xj1DMC[icent][ipt]->GetXaxis()->SetRangeUser(0.2,1);
	  h_xj1DMC[icent][ipt]->GetXaxis()->SetTitle("x_{J}");
	  h_xj1DMC[icent][ipt]->GetXaxis()->SetLabelSize(0);   // hide x labels (shown in ratio pad)
	  h_xj1DMC[icent][ipt]->GetXaxis()->SetTitle("");
	  h_xj1DMC[icent][ipt]->Draw();
	  hleg->AddEntry(h_xj1DMC[icent][ipt],"Pythia+Hijing","p");

	  h_xjdata[icent]->GetYaxis()->SetRange(ipt+1,ipt+1);
          h_xj1Ddata[icent][ipt] = (TH1F*)h_xjdata[icent]->ProjectionX();
          h_xj1Ddata[icent][ipt]->SetName(Form("h_xj_cent%i_pt%i_data",icent,ipt));
          h_xj1Ddata[icent][ipt]->Scale(1./h_xj1Ddata[icent][ipt]->Integral(),"width");

          h_xj1Ddata[icent][ipt]->SetMarkerColor(colors[icent]);
          h_xj1Ddata[icent][ipt]->SetLineColor(colors[icent]);
          h_xj1Ddata[icent][ipt]->GetYaxis()->SetRangeUser(0,5);
          h_xj1Ddata[icent][ipt]->GetXaxis()->SetRangeUser(0.2,1);
          h_xj1Ddata[icent][ipt]->GetXaxis()->SetTitle("x_{J}");
          
          h_xj1Ddata[icent][ipt]->Draw("SAME");
          hleg->AddEntry(h_xj1Ddata[icent][ipt],"Data","p");

	  cleg->AddEntry("",Form("%2.1f < p_{T,1}^{calib} < %2.1f GeV",h_xjMC[0]->GetYaxis()->GetBinLowEdge(ipt+1),h_xjMC[0]->GetYaxis()->GetBinLowEdge(ipt+2)),"");
	  cleg->AddEntry("","p_{T,2}^{calib} > 12 GeV","");
	  cleg->AddEntry("",cent_str[icent].c_str(),"");

	  leg->Draw();
	  cleg->Draw();
	  hleg->Draw();

	  pad2->cd();
	  TH1D *hRatio = (TH1D*)h_xj1Ddata[icent][ipt]->Clone("hRatio");
	  hRatio->Divide(h_xj1DMC[icent][ipt]);
	  hRatio->SetMinimum(0.0);
	  hRatio->SetMaximum(2.0);
	  hRatio->GetXaxis()->SetTitle("x_{J}");
	  hRatio->GetXaxis()->SetTitleSize(0.12);
	  hRatio->GetXaxis()->SetTitleOffset(1.0);
	  hRatio->GetXaxis()->SetLabelSize(0.10);
	  hRatio->GetXaxis()->SetTickLength(0.07);

	  hRatio->GetYaxis()->SetTitle("Data/MC");
	  hRatio->GetYaxis()->SetTitleSize(0.12);
	  hRatio->GetYaxis()->SetTitleOffset(0.45);
	  hRatio->GetYaxis()->SetLabelSize(0.10);
	  hRatio->GetYaxis()->SetNdivisions(505);
    
	  hRatio->Draw();

	  
	  int firstBin = hRatio->GetXaxis()->GetFirst();
	  int lastBin  = hRatio->GetXaxis()->GetLast();

	  double xmin = hRatio->GetXaxis()->GetBinLowEdge(firstBin);
	  double xmax = hRatio->GetXaxis()->GetBinUpEdge(lastBin);
	  TLine *line = new TLine(xmin, 1.0,
                             xmax, 1.0);
	  line->SetLineColor(kBlack);
	  line->SetLineStyle(2);
	  line->Draw("SAME");

	  c->cd();
	  c->Update();
	  
	  c->Print(Form("plots/xj_pt%i_cent%i_dataMC.pdf",ipt,icent));

	  hleg->Clear();
	  cleg->Clear();
	}
    }
}
