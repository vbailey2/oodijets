#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

void drawunfoldclosure(bool ishalf = 1)
{

  SetsPhenixStyle();
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();

  string infile;
  if(ishalf) infile = "hists/projections_unfoldMC_Half.root";
  else infile = "hists/projections_unfoldMC.root";
  TFile *f = new TFile("hists/hist-unfoldedMC.root","READ");
  TFile *fproj = new TFile(infile.c_str(),"READ");
  
  
  TCanvas *c = new TCanvas("c","c",700,700);

  TCanvas *c2 = new TCanvas("c2","c2",700,700);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.01);
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
  leg->AddEntry("","Pythia + HIJING","");
  leg->AddEntry("","O+O #sqrt{s_{NN}} = 200 GeV","");
  leg->AddEntry("","anti-#it{k}_{#it{t}} #it{R} = 0.4, |#eta| < 0.7","");


  TLegend *cleg = new TLegend(.15,.6,.4,.75);
  cleg->SetFillStyle(0);
  
  TLegend *hleg = new TLegend(.6,.75,.9,.92);
  hleg->SetFillStyle(0);


  int colors[] = {1,2,4,kGreen+2, kViolet,kCyan,kOrange+2,kMagenta+2,kAzure-2};
  
 
  TH3F* h_bins = (TH3F*)f->Get("h_pt1pt2");
  int npt = h_bins->GetNbinsY();
  int ncent = h_bins->GetNbinsZ();


  std::string cent_str[] = {"0-20%","20-40%","40-60%","60-80%"};
  
  TH2F *h_xj[ncent];
  h_xj[0] = (TH2F*)fproj->Get("h_xjunfold_0");
  TH2D *hTrue2D[ncent];
  TH2D *hMeas2D[ncent];
  
  const int n_xj = h_xj[0]->GetNbinsX();
  const int n_final = h_xj[0]->GetNbinsY();
  
  for(int i = 0; i < ncent; i++)
    {
      h_xj[i] = (TH2F*)fproj->Get(Form("h_xjunfold_%i",i));
      hTrue2D[i] = (TH2D*)fproj->Get(Form("h_xjtrue_%i",i));
      hMeas2D[i] = (TH2D*)fproj->Get(Form("h_xjmeas_%i",i));
    }

  TH1F *h_xj1D[ncent][npt];
  TH1D *hTrue1D[ncent][npt];
  TH1D *hMeas1D[ncent][npt];
  
  c2->cd();
  for(int ipt = 0; ipt < n_final; ipt++)
    {
      for(int icent = 0; icent < ncent; icent++)
        {
	  pad1->cd();

	  //draw truth
	  hTrue2D[icent]->GetYaxis()->SetRange(ipt+1,ipt+1);
	  hTrue1D[icent][ipt] = (TH1D*)hTrue2D[icent]->ProjectionX();
	  hTrue1D[icent][ipt]->SetName(Form("h_xjtrue_cent%i_pt%i",icent,ipt));
          //hTrue1D[icent][ipt]->Scale(1./hTrue1D[icent][ipt]->Integral(),"width");
	  hTrue1D[icent][ipt]->SetMarkerColor(colors[icent]);
          hTrue1D[icent][ipt]->SetLineColor(colors[icent]);
	  //hTrue1D[icent][ipt]->GetYaxis()->SetRangeUser(0,5);
          hTrue1D[icent][ipt]->GetXaxis()->SetRangeUser(0.3,1);
	  hTrue1D[icent][ipt]->GetXaxis()->SetTitle("x_{J}");
	  hTrue1D[icent][ipt]->Draw();
	  hleg->AddEntry(hTrue1D[icent][ipt],"Truth","p");
	  
	  //draw unfold
	  h_xj[icent]->GetYaxis()->SetRange(ipt+1,ipt+1);
          h_xj1D[icent][ipt] = (TH1F*)h_xj[icent]->ProjectionX();
          h_xj1D[icent][ipt]->SetName(Form("h_xj_cent%i_pt%i",icent,ipt));
          //h_xj1D[icent][ipt]->Scale(1./h_xj1D[icent][ipt]->Integral(),"width");

          h_xj1D[icent][ipt]->SetMarkerColor(colors[icent]);
	  h_xj1D[icent][ipt]->SetMarkerStyle(25);
          h_xj1D[icent][ipt]->SetLineColor(colors[icent]);
          h_xj1D[icent][ipt]->GetYaxis()->SetRangeUser(0,5);
          h_xj1D[icent][ipt]->GetXaxis()->SetRangeUser(0.3,1);
          h_xj1D[icent][ipt]->GetXaxis()->SetTitle("x_{J}");
	  h_xj1D[icent][ipt]->Draw("SAME");
	  
	  hleg->AddEntry(h_xj1D[icent][ipt],"Unfold","p");
	  cleg->AddEntry("",Form("%2.1f < p_{T} < %2.1f GeV",h_xj[0]->GetYaxis()->GetBinLowEdge(ipt+1),	h_xj[0]->GetYaxis()->GetBinLowEdge(ipt+2)),"");
	  cleg->AddEntry("",cent_str[icent].c_str(),"");
	  
	  hleg->Draw();
	  cleg->Draw();
	  leg->Draw();
	  
	  pad2->cd();
	  TH1D *hRatio = (TH1D*)h_xj1D[icent][ipt]->Clone("hRatio");
	  hRatio->Divide(hTrue1D[icent][ipt]);
	  hRatio->SetTitle("");
	  
	  hRatio->SetMinimum(0.0);
	  hRatio->SetMaximum(2.0);
	  hRatio->GetXaxis()->SetTitle("x_{J}");
	  hRatio->GetXaxis()->SetTitleSize(0.12);
	  hRatio->GetXaxis()->SetTitleOffset(1.0);
	  hRatio->GetXaxis()->SetLabelSize(0.10);
	  hRatio->GetXaxis()->SetTickLength(0.07);

	  hRatio->GetYaxis()->SetTitle("Unfold/Truth");
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

	  if(ishalf) c2->Print(Form("plots/xj_unfoldhalfclosure_pt%i_cent%i.pdf",ipt,icent));
	  else c2->Print(Form("plots/xj_unfoldfullclosure_pt%i_cent%i.pdf",ipt,icent));
	  
	  hleg->Clear();
	  cleg->Clear();
	}
    }
}
