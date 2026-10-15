#ifndef truth_sumet_weight_h
#define truth_sumet_weight_h

// Truth-ΣE_T → data/MC weight helper for MC analyses.
//
// For each MC event with a known truth_sumet, this returns the data/MC weight
// to apply, derived from the response + resolution fits in
// TruthETWeight/truth_et_response_fit.root and the best calo scale +
// data/MC ratio in Corrections/mb_calo_scale_fit.root.
//
// Per event the procedure is:
//   1. mu  = f_sumet_response(truth_sumet)
//   2. sig = f_sumet_resolution_calo(truth_sumet) * mu
//   3. r   = TRandom3->Gaus(mu, sig)
//   4. reco_smear_scaled = best_calo_scale * truth_sumet * r
//   5. weight = r_best_calo->GetBinContent(FindBin(reco_smear_scaled))
//   6. return (weight > 0) ? weight : 1.0
//
// Out-of-range (truth_sumet ≤ 0 or large): the TF1s are evaluated freely
// (extrapolation). Caller is responsible for any pre-filtering.
//
// RNG: a single TRandom3 with default seed 20260526 is used. Override via
// SetTruthSumETSmearSeed(seed) once before the event loop if needed.
//
// Usage:
//   #include "/sphenix/user/jpark4/Analysis/Oxygen/sPHENIX-OxygenAnalysis/Headers/truth_sumet_weight.h"
//   ...
//   const double w = GetTruthSumETSmearWeight(truth_sumet);
//   h->Fill(some_observable, w);

#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TParameter.h>
#include <TRandom3.h>
#include <cstdlib>
#include <iostream>

inline const char *DefaultTruthSumETResponseFitFile()
{
   return "/gpfs/mnt/gpfs02/sphenix/user/jpark4/Analysis/Oxygen/sPHENIX-OxygenAnalysis/TruthETWeight/truth_et_response_fit.root";
}

inline const char *DefaultTruthSumETBestScaleFile()
{
   return "/gpfs/mnt/gpfs02/sphenix/user/jpark4/Analysis/Oxygen/sPHENIX-OxygenAnalysis/Corrections/mb_calo_scale_fit.root";
}

inline TRandom3 &TruthSumETSmearRng_()
{
   static TRandom3 rng(20260526u);
   return rng;
}

inline void SetTruthSumETSmearSeed(unsigned int seed)
{
   TruthSumETSmearRng_().SetSeed(seed);
}

inline double GetTruthSumETSmearWeight(double truth_sumet)
{
   static bool   initialized      = false;
   static TF1   *f_resp           = nullptr;
   static TF1   *f_res            = nullptr;
   static TH1D  *r_best_calo      = nullptr;
   static double best_calo_scale  = 1.;

   if (!initialized) {
      // Response + resolution TF1s
      const char *fit_path = DefaultTruthSumETResponseFitFile();
      TFile *ff = TFile::Open(fit_path, "READ");
      if (!ff || ff->IsZombie()) {
         std::cerr << "ERROR: cannot open " << fit_path << std::endl;
         std::exit(1);
      }
      TF1 *r_in = dynamic_cast<TF1 *>(ff->Get("f_sumet_response"));
      TF1 *e_in = dynamic_cast<TF1 *>(ff->Get("f_sumet_resolution_calo"));
      if (!r_in || !e_in) {
         std::cerr << "ERROR: missing f_sumet_response / f_sumet_resolution_calo in "
                   << fit_path << std::endl;
         ff->Close();
         std::exit(1);
      }
      f_resp = static_cast<TF1 *>(r_in->Clone("f_sumet_response_cached"));
      f_res  = static_cast<TF1 *>(e_in->Clone("f_sumet_resolution_cached"));
      ff->Close();

      // best calo scale + data/MC ratio
      const char *scale_path = DefaultTruthSumETBestScaleFile();
      TFile *fs = TFile::Open(scale_path, "READ");
      if (!fs || fs->IsZombie()) {
         std::cerr << "ERROR: cannot open " << scale_path << std::endl;
         std::exit(1);
      }
      TParameter<double> *bp = dynamic_cast<TParameter<double> *>(fs->Get("best_scale_value"));
      if (!bp) {
         std::cerr << "ERROR: missing best_scale_value in " << scale_path << std::endl;
         fs->Close();
         std::exit(1);
      }
      best_calo_scale = bp->GetVal();
      TH1D *rb = dynamic_cast<TH1D *>(fs->Get("r_best_calo"));
      if (!rb) {
         std::cerr << "ERROR: missing r_best_calo in " << scale_path << std::endl;
         fs->Close();
         std::exit(1);
      }
      r_best_calo = static_cast<TH1D *>(rb->Clone("r_best_calo_cached"));
      r_best_calo->SetDirectory(nullptr);
      fs->Close();

      initialized = true;
      std::cout << "TruthSumETWeight initialized: best_calo_scale = "
                << best_calo_scale << std::endl;
   }

   const double mu  = f_resp->Eval(truth_sumet);
   const double sig = f_res ->Eval(truth_sumet) * mu;
   const double r   = TruthSumETSmearRng_().Gaus(mu, sig);
   const double reco_smear_scaled = best_calo_scale * truth_sumet * r;

   const int bin = r_best_calo->GetXaxis()->FindBin(reco_smear_scaled);
   const double w = r_best_calo->GetBinContent(bin);
   return (w > 0.) ? w : 1.;
}

#endif
