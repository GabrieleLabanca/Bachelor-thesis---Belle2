#include "v1742.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
using namespace std;

#include <assert.h>
#include "TMath.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH1.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TF1.h"
#include <TFitResultPtr.h>

// classes by Gabriele Labanca
#include "classes/digit_utilities.h"


int prog07_fit_final(TString fname) {
  ifstream is(fname+".dat");
  V1742 v(4);
  
  ////////////////
  // HISTOGRAMS //
  //////////////// 
  TH1D * hTrigStart = new TH1D("hTrigStart", "t0 times", 1025, -1., 1024.);
  TH2D * hScopeAll = new TH2D("hScopeAll", "scope histogram", 1250, 0., 1250., 1800, -1200., 600.);
  TH2D * hScopeT860 = new TH2D("hScopeT860", "scope histogram, events at 860 peak of time of arrival", 1250, 0., 1250., 1800, -1200., 600.);
  TH2D * hScopeTemp = new TH2D("hScopeTemp", "temporary, scope histogram", 1250, 0., 1250., 1800, -1200., 600.);
  TH2D * hScopeT860_verypho = new TH2D("hScopeT860_verypho", "scope histogram, events at 860 peak of time of arrival, max>75",1250, 0., 1250., 1800, -1200., 600.);
  TH2D * hScopeT860_pho = new TH2D("hScopeT860_pho", "scope histogram, events at 860 peak of time of arrival, 50<max<75",1250, 0., 1250., 1800, -1200., 600.);
  TH2D * hScopeT860_phobkg = new TH2D("hScopeT860_phobkg", "scope histogram, events at 860 peak of time of arrival, 25<max<50",1250, 0., 1250., 1800, -1200., 600.);
  TH2D * hScopeT860_bkg = new TH2D("hScopeT860_bkg", "scope histogram, events at 860 peak of time of arrival, max<25",1250, 0., 1250., 1800, -1200., 600.);
  TH2D * hScope_badmatch = new TH2D("hScope_badmatch","Scope histogram, events mismatched by amplitude-integral(sum) cross-check",1250, 0., 1250., 1800, -1200., 600.);


  ///////////////////////
  // FIT INTEGRAL HISTOGRAMS
  //  - length of integration
  //  - from real minimum or from +50
  //  - avg or not avg
  //  - no errors or errors-not_avg or errors/sqrt(N)-avg
  //
  //  NAMES: hFitIntegral + Bare/Avg + Min/Retarded (+ Errors)
  //////////////////////
  TH1D* hFitIntgrBareMin500 = new TH1D("hFitIntgrBareMin500", 
                                "Integral of the fitting function from 0 to 500 bins after the minimum", 
                                1500, -10000., 1000000.);
  // start point: minimum+50bins
  TH1D* hFitIntgrBareRetarded500 = new TH1D("hFitIntgrBareRetarded500",
                                "Integral of the fitting function from 50 to 550 bins after the minimum",
                                1500, -10000., 1000000.);
  // averaged points
  TH1D* hFitIntgrAvgMin500 = new TH1D("hFitIntgrAvgMin500",
                                "Integral of the fitting function from 0 to 500 bins after the minimum, averaged points",
                                1500, -10000., 1000000.);
  // errors
  TH1D* hFitIntgrBareMinErrors500 = new TH1D("hFitIntgrBareMinErrors500",
                                "Integral of the fitting function from 0 to 500 bins after the minimum, fit with errors",
                                1500, -10000., 1000000.);
  TH1D* hFitIntgrAvgMinErrors500 = new TH1D("hFitIntgrAvgMinErrors500",
                                "Integral of the fitting function from 0 to 500 bins after the minimum, averaged points, fit with errors on the averages",
                                1500, -10000., 1000000.);



  ///////////////
  // VARIABLES //
  ///////////////
  // time when trigger starts
  double t0=-999;
  // bin start signal
  int sign_start;
  // extimated distance between trigger start and signal start
  int central_bin = 290;
  // count good/bad events
  int counter=0; int counter_good=0; int counter_bad=0;
  
  
  // FIT //
  // number of fit points
  int const n_fit = 50;
  // arrays for fit
  //   amplitudes
  double fit_points[n_fit];
  //   bins
  double fit_bins[n_fit];
  //   errors
  double fit_errorsy[n_fit];
  double fit_errorsx[n_fit];
  double fit_zeros[n_fit]={0}; // errors on bin
  // length of integration
  int intgr_length = 500;
  // #bins over which the mean is calculated = 2*n_mean+1
  int n_mean = 2;

  ///////////////
  // ROOT FILE //
  /////////////// TODO
  TFile *f_data = new TFile(fname+"_07.root","recreate");
  cout << "Writing on file ''" << fname << "_07.03.root''" << endl;
  TTree *main_data = new TTree("main_data","input_data_fit");
  
  // fit trees
  TTree * fit_noerr_noavg = new TTree("fit_noerr_noavg","fit_noerr_noavg");
  TTree * fit_err_noavg   = new TTree("fit_err_noavg"  ,"fit_err_noavg"  );
  TTree * fit_err_avg     = new TTree("fit_err_avg"    , "fit_err_avg"   );
 


  
  ////////////
  // LEAVES //
  ////////////
  // fit parameters: [0]+exp([1]+[2]*x)
  double ampli_t0_290=-999;
  double Background=-999;
  double Bkg=-999;
  double bkg_rms = -999.;
  double bkg_hampli = -999.;
  double_t time_min = -999; // time of arrival after the trigger 
  double min = -999;
  main_data->Branch("ampli_t0_290",&ampli_t0_290,"ampli_t0_290/D");
  fit_noerr_noavg->Branch("ampli_t0_290",&ampli_t0_290,"ampli_t0_290/D");
  fit_err_noavg->Branch("ampli_t0_290",&ampli_t0_290,"ampli_t0_290/D");
  fit_err_avg->Branch("ampli_t0_290",&ampli_t0_290,"ampli_t0_290/D");
  
  main_data->Branch("Background",&Background,"Background/D");
  main_data->Branch("Bkg",&Bkg,"Bkg/D");
  main_data->Branch("time_min",&time_min,"time_min/D");
  main_data->Branch("min",&min,"min/D");
  fit_noerr_noavg->Branch("Background",&Background,"Background/D");
  fit_noerr_noavg->Branch("Bkg",&Bkg,"Bkg/D");
  fit_noerr_noavg->Branch("time_min",&time_min,"time_min/D");
  fit_noerr_noavg->Branch("min",&min,"min/D");
  fit_err_noavg->Branch("Background",&Background,"Background/D");
  fit_err_noavg->Branch("Bkg",&Bkg,"Bkg/D");
  fit_err_noavg->Branch("time_min",&time_min,"time_min/D");
  fit_err_noavg->Branch("min",&min,"min/D");
  fit_err_avg->Branch("Background",&Background,"Background/D");
  fit_err_avg->Branch("Bkg",&Bkg,"Bkg/D");
  fit_err_avg->Branch("time_min",&time_min,"time_min/D");
  fit_err_avg->Branch("min",&min,"min/D");

  main_data->Branch("bkg_rms",&bkg_rms,"bkg_rms/D");
  main_data->Branch("bkg_hampli",&bkg_hampli,"bkg_hampli/D");
  fit_noerr_noavg->Branch("bkg_rms",&bkg_rms,"bkg_rms/D");
  fit_noerr_noavg->Branch("bkg_hampli",&bkg_hampli,"bkg_hampli/D");
  fit_err_noavg->Branch("bkg_rms",&bkg_rms,"bkg_rms/D");
  fit_err_noavg->Branch("bkg_hampli",&bkg_hampli,"bkg_hampli/D");
  fit_err_avg->Branch("bkg_rms",&bkg_rms,"bkg_rms/D");
  fit_err_avg->Branch("bkg_hampli",&bkg_hampli,"bkg_hampli/D");

  
  main_data->Branch("t0",&t0,"t0/D");

  double fit_A = -999;
  double fit_k  = -999;
  int fit_stat = -999;
  int fit_chi2 = -999;
  double integral_min_300;
  double integral_min_400;
  double integral_min_500;


  int ampli_min_860 = -999;
  int integral_all = -999;
  int integral_bins_400 = -999;
  int integral_from_820 = -999;
  main_data->Branch("ampli_min_860",&ampli_min_860,"ampli_min_860/I");
  main_data->Branch("integral_all",&integral_all,"integral_all/I");
  main_data->Branch("integral_bins_400",&integral_bins_400,"integral_bins_400/I");
  main_data->Branch("integral_from_820",&integral_from_820,"integral_from_820/I");
  fit_noerr_noavg->Branch("ampli_min_860",&ampli_min_860,"ampli_min_860/I");
  fit_noerr_noavg->Branch("integral_all",&integral_all,"integral_all/I");
  fit_noerr_noavg->Branch("integral_bins_400",&integral_bins_400,"integral_bins_400/I");
  fit_noerr_noavg->Branch("integral_from_820",&integral_from_820,"integral_from_820/I");
  fit_err_noavg->Branch("ampli_min_860",&ampli_min_860,"ampli_min_860/I");
  fit_err_noavg->Branch("integral_all",&integral_all,"integral_all/I");
  fit_err_noavg->Branch("integral_bins_400",&integral_bins_400,"integral_bins_400/I");
  fit_err_noavg->Branch("integral_from_820",&integral_from_820,"integral_from_820/I");
  fit_err_avg->Branch("ampli_min_860",&ampli_min_860,"ampli_min_860/I");
  fit_err_avg->Branch("integral_all",&integral_all,"integral_all/I");
  fit_err_avg->Branch("integral_bins_400",&integral_bins_400,"integral_bins_400/I");
  fit_err_avg->Branch("integral_from_820",&integral_from_820,"integral_from_820/I");

  fit_noerr_noavg->Branch("fit_k",&fit_k,"fit_k/D");
  fit_noerr_noavg->Branch("fit_A",&fit_A,"fit_A/D");  
  fit_noerr_noavg->Branch("fit_stat",&fit_stat,"fit_stat/I");
  fit_noerr_noavg->Branch("fit_chi2",&fit_chi2,"fit_chi2/I");
  fit_noerr_noavg->Branch("integral_min_300",&integral_min_300,"integral_min_300/I");
  fit_noerr_noavg->Branch("integral_min_400",&integral_min_400,"integral_min_400/I");
int integral_min_500_NN = -999;
  fit_noerr_noavg->Branch("integral_min_500",&integral_min_500,"integral_min_500/I");

  fit_err_noavg->Branch("fit_k",&fit_k,"fit_k/D");
  fit_err_noavg->Branch("fit_A",&fit_A,"fit_A/D");  
  fit_err_noavg->Branch("fit_stat",&fit_stat,"fit_stat/I");
  fit_err_noavg->Branch("fit_chi2",&fit_chi2,"fit_chi2/I");
  fit_err_noavg->Branch("integral_min_300",&integral_min_300,"integral_min_300/I");
  fit_err_noavg->Branch("integral_min_400",&integral_min_400,"integral_min_400/I");
  fit_err_noavg->Branch("integral_min_500",&integral_min_500,"integral_min_500/I");

  fit_err_avg->Branch("fit_k",&fit_k,"fit_k/D");
  fit_err_avg->Branch("fit_A",&fit_A,"fit_A/D");  
  fit_err_avg->Branch("fit_stat",&fit_stat,"fit_stat/I");
  fit_err_avg->Branch("fit_chi2",&fit_chi2,"fit_chi2/I");
  fit_err_avg->Branch("integral_min_300",&integral_min_300,"integral_min_300/D");
  fit_err_avg->Branch("integral_min_400",&integral_min_400,"integral_min_400/D");
  fit_err_avg->Branch("integral_min_500",&integral_min_500,"integral_min_500/D");


  double zero_value_fit = -999;
  fit_noerr_noavg->Branch("zero_value_fit",&zero_value_fit,"zero_value_fit/D");
  fit_err_noavg->Branch("zero_value_fit",&zero_value_fit,"zero_value_fit/D");
  fit_err_avg->Branch("zero_value_fit",&zero_value_fit,"zero_value_fit/D");




  // while read from file - loops over the events
  while ((is >> v))  { //TODO
    // start time of trigger
    t0 = v.triggercft(8);
    sign_start = t0 + central_bin;

    
    // declare digit: vector of entries
    vector<float> digit = v.waveform(0); 
    counter++;
    
    /////////////////
    // good events //
    /////////////////
    if (digit.size()!=0){
      
      // creates a digit_utilities object
      digit_utilities mydig (&digit);


      // BKG mean first 100 bins //
      counter_good++;
      double av = 0.; // average of the background
      double av_rms = 0.; // standard deviation of the background
      double av_hamp = 0.; // half of the maximum (amp_+ - amp_-)
      
      for (int i=10; i<110; ++i) {
        av += digit[i];
      }
      av *= 0.01; 
       
      for (int i=100; i<200; ++i) {
        av_rms += (av - digit[i]) * (av - digit[i]);
      }
      av_rms = sqrt(av_rms)/10;
      bkg_rms = av_rms;      
      double hampP = 0.;
      double hampN = 0.;
      double temp_rel_amp;
      for (int i=10; i<110; ++i) {
        temp_rel_amp = digit[i]-av;
        if(temp_rel_amp > (hampP)) hampP = temp_rel_amp;
        if(temp_rel_amp < (hampN)) hampN = temp_rel_amp;
      }
      av_hamp = (hampP - hampN)/2.;
      bkg_hampli = av_hamp;


      integral_all = mydig.integral(0,1024,av);
      integral_bins_400 = mydig.integral(sign_start,400,av);
      integral_from_820 = mydig.integral(820,1024-820,av);



      min = 100000.;
      int mindex=0;
      // fill hScopeAll
      //   x: ranges over digit.size scaled to the same start
      //   y: the corresponding digit value - average subtracted
      for (int i=0; i<digit.size(); ++i) {
	hScopeAll->Fill(double(i)-t0+200., digit[i]-av);
	// looks for absolute minimum 
        if (digit[i]<min){
	  min=digit[i]; // min
	  mindex=i;     // bin of min
	}
      }
      min = av - min;
      // time of arrival after the trigger
      time_min = mindex - t0;
      // fills the TH2D with events at the peak of mintime at 860
      for (int i=0; i<digit.size(); ++i) {
        if(time_min>835 && time_min<865) {
          hScopeT860->Fill(double(i)-t0+200., digit[i]-av);   

          ampli_min_860 = mydig.min(t0+280,330-280);
          if     ( (av-ampli_min_860) > 75)                        hScopeT860_verypho->Fill(double(i)-t0+200., digit[i]-av);
          else if( (av-ampli_min_860) > 50 && (av-ampli_min_860) < 75) hScopeT860_pho->Fill(double(i)-t0+200., digit[i]-av);
          else if( (av-ampli_min_860) > 25 && (av-ampli_min_860)<50)   hScopeT860_phobkg->Fill(double(i)-t0+200., digit[i]-av);
          else if( (av-ampli_min_860) < 25)                        hScopeT860_bkg->Fill(double(i)-t0+200., digit[i]-av);
        }
      }

      // events mismatched by amplitude-integral(sum) cross-check
      for (int i=0; i<digit.size(); ++i) {
        if(ampli_t0_290<50 && integral_bins_400>3000) hScope_badmatch->Fill(double(i)-t0+200., digit[i]-av);
      }


      // cft TODO tindex
      int tindex=0;
      for (int j=mindex; 
           j>mindex-20 and mindex>0 and mindex<digit.size(); 
          --j){ 	
      	if (digit[j]-av>(min-av)*0.5) {
	  tindex=j;
	  break;
	}
      }
      
      double time = tindex + (digit[tindex]-min*0.5)/(digit[tindex+1]-digit[tindex]);

   


      // fill hAmpli
      ampli_t0_290=-digit[t0+central_bin]+av;/////Questa e' l'ampiezza che mi interessa






      // FIT //
      for (int i=0; i<n_fit; i++) fit_bins[i] = (i*500./n_fit);
      for (int i=0; i<n_fit; i++) fit_errorsx[i] = 1./sqrt(6); // triangular error
      TFitResultPtr r;
      // no error, not averaged // 
      // extracts fit coordinates
      for(int i=0; i<n_fit; i++) fit_points[i] = (av - digit.at(sign_start + i*500./n_fit) );
      // graph
      TGraph* gFitNeNa = new TGraph(n_fit,fit_bins,fit_points);
      // function
      TF1* fFitNeNa = new TF1("fFitNeNa","exp([0]-[1]*x)",0,500.0);
      r = ( gFitNeNa->Fit(fFitNeNa,"S:Q") );  // this to read the status of the fit: 0=ok, 4=bad
      fit_A = fFitNeNa->GetParameter(0);
      fit_k = fFitNeNa->GetParameter(1);
      fit_stat = (Int_t) r->Status();
      fit_chi2 = fFitNeNa->GetChisquare();
      integral_min_300 = fFitNeNa->Integral(0,300);
      integral_min_400 = fFitNeNa->Integral(0,400);
      integral_min_500 = fFitNeNa->Integral(0,500);
      zero_value_fit = fFitNeNa->Eval(0);
      fit_noerr_noavg->Fill();

  
      


      // error, not averaged // 
      // extracts fit coordinates
      for(int i=0; i<n_fit; i++) fit_points[i] = (av - digit.at(sign_start + i*500./n_fit) );
      for(int i=0; i<n_fit; i++) fit_errorsy[i] = bkg_rms;
      // graph
      TGraphErrors* gFitENa = new TGraphErrors(n_fit,fit_bins,fit_points,fit_errorsx,fit_errorsy);
      // function
      TF1* fFitENa = new TF1("fFitENa","exp([0]-[1]*x)",0,500.0);
      r = ( gFitENa->Fit(fFitENa,"S:Q") );  // this to read the status of the fit: 0=ok, 4=bad
      fit_A = fFitENa->GetParameter(0);
      fit_k = fFitENa->GetParameter(1);
      fit_stat = (Int_t) r->Status();
      fit_chi2 = fFitENa->GetChisquare();
      integral_min_300 = fFitENa->Integral(0,300);
      integral_min_400 = fFitENa->Integral(0,400);
      integral_min_500 = fFitENa->Integral(0,500);
      zero_value_fit = fFitENa->Eval(0);
      fit_err_noavg->Fill();






      // error, averaged // 
      // extracts fit coordinates
      for (int i=0; i<n_fit; i++) fit_errorsx[i] = 2./sqrt(6); // triangular error
      for(int i=0; i<n_fit; i++) fit_points[i] = (av - mydig.mean(sign_start + i*500./n_fit,n_mean) );
      for(int i=0; i<n_fit; i++) fit_errorsy[i] = bkg_rms/sqrt(double(2*n_mean+1));
      // graph
      TGraphErrors* gFitEA = new TGraphErrors(n_fit,fit_bins,fit_points,fit_errorsx,fit_errorsy);
      // function
      TF1* fFitEA = new TF1("fFitEA","exp([0]-[1]*x)",0,500.0);
      r = ( gFitEA->Fit(fFitEA,"S:Q") );  // this to read the status of the fit: 0=ok, 4=bad
      fit_A = fFitEA->GetParameter(0);
      fit_k = fFitEA->GetParameter(1);
      fit_stat = (Int_t) r->Status();
      fit_chi2 = fFitEA->GetChisquare();
      integral_min_300 = fFitEA->Integral(0,300);
      integral_min_400 = fFitEA->Integral(0,400);
      integral_min_500 = fFitEA->Integral(0,500);
      zero_value_fit = fFitEA->Eval(0);
      fit_err_avg->Fill();








      Background = -digit[50]+av;
      Bkg = -digit[t0-50]+av;
      
      int indx = v.startindex(31); //returns begin signal of trigger
      main_data->Fill();
      if(counter%1000 == 0) cerr << "##### NUMBER " << counter << " #####" << endl;
    } else {counter_bad++;}
  }
  
 

  cout<<"N tot  : "<<counter<<endl;
  cout<<"N good : "<<counter_good<<endl;
  cout<<"N bad  : "<<counter_bad<<endl;







  //////////////////////////////
  // write on the file f_data //
  //////////////////////////////
  f_data->cd();

  main_data->Write();

  fit_noerr_noavg->Write();
  fit_err_noavg->Write();
  fit_err_avg->Write();

  hScopeAll->Write();
  hScopeT860->Write(); 
  hScopeT860_verypho->Write(); 
  hScopeT860_pho->Write(); 
  hScopeT860_phobkg->Write(); 
  hScopeT860_bkg->Write(); 
  hScope_badmatch->Write();

  
  // amplitudes and integrals
  TH1D * h_integral_bins_400 = new TH1D("h_integral_bins_400","Integral summing over 400 bins",1000,-100,200000);
  main_data->Draw("integral_bins_400>>h_integral_bins_400");
  h_integral_bins_400->Write();
  TH2D * h_ampli_integral = new TH2D("h_ampli_integral","Amplitude at 290(sum integral)",500,0,220000,500,0,1600);
  main_data->Draw("ampli_t0_290:integral_bins_400>>h_ampli_integral","","colz");
  h_ampli_integral->Write(); 
  



  // NO ERRORS, NO AVERAGE
  f_data->mkdir("noerr_noavg");
  f_data->cd("noerr_noavg");
  // A
  TH1D * hNN_A = new TH1D("hNN_A","Fit no-errors, no-average: A",100,-40,40);
  fit_noerr_noavg->Draw("fit_A>>hNN_A","fit_A>-40&&fit_A<40");
  TH1D * hNN_A_pho = new TH1D("hNN_A_pho","Fit no-errors, no-average: A; photon events",100,0,10);
  fit_noerr_noavg->Draw("fit_A>>hNN_A_pho","fit_A>-40&&fit_A<40 && ampli_t0_290 > 25");
  TH1D * hNN_A_bkg = new TH1D("hNN_A_bkg","Fit no-errors, no-average: A; background events",100,-40,40);
  fit_noerr_noavg->Draw("fit_A>>hNN_A_bkg","fit_A>-40&&fit_A<40 && ampli_t0_290 < 25");
  TH1D * hNN_A_good = new TH1D("hNN_A_good","Fit no-errors, no-average: A; good fit",100,-40,40);
  fit_noerr_noavg->Draw("fit_A>>hNN_A_good","fit_A>-40&&fit_A<40 && fit_stat == 0");
  TH1D * hNN_A_bad = new TH1D("hNN_A_bad","Fit no-errors, no-average: A; bad fit",100,-40,40);
  fit_noerr_noavg->Draw("fit_A>>hNN_A_bad","fit_A>-40&&fit_A<40 && fit_stat == 4");
  TH1D * hNN_A_PG = new TH1D("hNN_A_PG","Fit no-errors, no-average: A; photon events, good fit",100,0,10);
  fit_noerr_noavg->Draw("fit_A>>hNN_A_PG","fit_A>-40&&fit_A<40 && ampli_t0_290 > 25 && fit_stat == 0");
  TH1D * hNN_A_PB = new TH1D("hNN_A_PB","Fit no-errors, no-average: A; photon events, bad fit",100,-40,40);
  fit_noerr_noavg->Draw("fit_A>>hNN_A_PB","fit_A>-40&&fit_A<40 && ampli_t0_290 > 25 && fit_stat == 4");
  TH1D * hNN_A_BG = new TH1D("hNN_A_BG","Fit no-errors, no-average: A; background events, good fit",100,-40,40);
  fit_noerr_noavg->Draw("fit_A>>hNN_A_BG","fit_A>-40&&fit_A<40 && ampli_t0_290 < 25 && fit_stat == 0");
  TH1D * hNN_A_BB = new TH1D("hNN_A_BB","Fit no-errors, no-average: A; background events, bad fit",100,-40,40);
  fit_noerr_noavg->Draw("fit_A>>hNN_A_BB","fit_A>-40&&fit_A<40 && ampli_t0_290 < 25 && fit_stat == 4");
  // k
  TH1D * hNN_k = new TH1D("hNN_k","Fit no-errors, no-average: k",100,-0.05,0.2);
  fit_noerr_noavg->Draw("fit_k>>hNN_k","fit_k>-0.05&&fit_k<0.2");
  TH1D * hNN_k_pho = new TH1D("hNN_k_pho","Fit no-errors, no-average: k; photon events",100,-0.05,0.2);
  fit_noerr_noavg->Draw("fit_k>>hNN_k_pho","fit_k>-0.05&&fit_k<0.2 && ampli_t0_290 > 25");
  TH1D * hNN_k_bkg = new TH1D("hNN_k_bkg","Fit no-errors, no-average: k; background events",100,-0.05,0.2);
  fit_noerr_noavg->Draw("fit_k>>hNN_k_bkg","fit_k>-0.05&&fit_k<0.2 && ampli_t0_290 < 25");
  TH1D * hNN_k_good = new TH1D("hNN_k_good","Fit no-errors, no-average: k; good fit",100,-0.05,0.2);
  fit_noerr_noavg->Draw("fit_k>>hNN_k_good","fit_k>-0.05&&fit_k<0.2 && fit_stat == 0");
  TH1D * hNN_k_bad = new TH1D("hNN_k_bad","Fit no-errors, no-average: k; bad fit",100,-0.05,0.2);
  fit_noerr_noavg->Draw("fit_k>>hNN_k_bad","fit_k>-0.05&&fit_k<0.2 && fit_stat == 4");
  TH1D * hNN_k_PG = new TH1D("hNN_k_PG","Fit no-errors, no-average: k; photon events, good fit",100,-0.05,0.2);
  fit_noerr_noavg->Draw("fit_k>>hNN_k_PG","fit_k>-0.05&&fit_k<0.2 && ampli_t0_290 > 25 && fit_stat == 0");
  TH1D * hNN_k_PB = new TH1D("hNN_k_PB","Fit no-errors, no-average: k; photon events, bad fit",100,-0.05,0.2);
  fit_noerr_noavg->Draw("fit_k>>hNN_k_PB","fit_k>-0.05&&fit_k<0.2 && ampli_t0_290 > 25 && fit_stat == 4");
  TH1D * hNN_k_BG = new TH1D("hNN_k_BG","Fit no-errors, no-average: k; background events, good fit",100,-0.05,0.2);
  fit_noerr_noavg->Draw("fit_k>>hNN_k_BG","fit_k>-0.05&&fit_k<0.2 && ampli_t0_290 < 25 && fit_stat == 0");
  TH1D * hNN_k_BB = new TH1D("hNN_k_BB","Fit no-errors, no-average: k; background events, bad fit",100,-0.05,0.2);
  fit_noerr_noavg->Draw("fit_k>>hNN_k_BB","fit_k>-0.05&&fit_k<0.2 && ampli_t0_290 < 25 && fit_stat == 4");
  // k(A)
  TH2D * hNN_kA = new TH2D("hNN_kA","Fit no-errors, no-average: k(A)",1000,-40,40,1000,-0.05,0.2);
  fit_noerr_noavg->Draw("fit_k:fit_A>>hNN_kA","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40","colz");
  TH2D * hNN_kA_pho = new TH2D("hNN_kA_pho", "Fit no-errors, no-average: k(A); photon events",1000,-40,40,1000,-0.05,0.2);
  fit_noerr_noavg->Draw("fit_k:fit_A>>hNN_kA_pho","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && ampli_t0_290 > 25","colz");
  TH2D * hNN_kA_bkg = new TH2D("hNN_kA_bkg", "Fit no-errors, no-average: k(A); background events",1000,-40,40,1000,-0.05,0.2);
  fit_noerr_noavg->Draw("fit_k:fit_A>>hNN_kA_bkg","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && ampli_t0_290 < 25","colz");
  TH2D * hNN_kA_good = new TH2D("hNN_kA_good", "Fit no-errors, no-average: k(A); good fit",1000,-40,40,1000,-0.05,0.2);
  fit_noerr_noavg->Draw("fit_k:fit_A>>hNN_kA_good","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && fit_stat == 0","colz");
  TH2D * hNN_kA_bad = new TH2D("hNN_kA_bad","Fit no-errors, no-average: k(A); bad  fit",1000,-40,40,1000,-0.05,0.2);
  fit_noerr_noavg->Draw("fit_k:fit_A>>hNN_kA_bad","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && fit_stat == 4","colz");
  TH2D * hNN_kA_PG = new TH2D("hNN_kA_PG","Fit no-errors, no-average: k(A); photon events, good fit",1000,-40,40,1000,-0.05,0.2);
  fit_noerr_noavg->Draw("fit_k:fit_A>>hNN_kA_PG","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && ampli_t0_290 > 25&& fit_stat == 0","colz");
  TH2D * hNN_kA_PB = new TH2D("hNN_kA_PB","Fit no-errors, no-average: k(A); photon events, bad  fit",1000,-40,40,1000,-0.05,0.2);
  fit_noerr_noavg->Draw("fit_k:fit_A>>hNN_kA_PB","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && ampli_t0_290 > 25&& fit_stat == 4","colz");
  TH2D * hNN_kA_BG = new TH2D("hNN_kA_BG","Fit no-errors, no-average: k(A); background events, good fit",1000,-40,40,1000,-0.05,0.2);
  fit_noerr_noavg->Draw("fit_k:fit_A>>hNN_kA_BG","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && ampli_t0_290 < 25 && fit_stat == 0","colz");
  TH2D * hNN_kA_BB = new TH2D("hNN_kA_BB","Fit no-errors, no-average: k(A); background events, bad  fit",1000,-40,40,1000,-0.05,0.2);
  fit_noerr_noavg->Draw("fit_k:fit_A>>hNN_kA_BB","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && ampli_t0_290 < 25 && fit_stat == 4","colz");
  // chi2
  TH1D * hNN_chi2 = new TH1D("hNN_chi2","Fit no-errors, no-average: chi2",1000,0,300000);
  fit_noerr_noavg->Draw("fit_chi2>>hNN_chi2","fit_chi2<300000");
  // integrals
  TH1D * hNN_intgr_300 = new TH1D("hNN_intgr_300","Fit no-errors, no-average: integral 300 bins",1000,0,100000);
  fit_noerr_noavg->Draw("integral_min_300>>hNN_intgr_300");
  TH1D * hNN_intgr_400 = new TH1D("hNN_intgr_400","Fit no-errors, no-average: integral 400 bins",1000,0,100000);
  fit_noerr_noavg->Draw("integral_min_400>>hNN_intgr_400");
  TH1D * hNN_intgr_500 = new TH1D("hNN_intgr_500","Fit no-errors, no-average: integral 500 bins",1000,0,100000);
  fit_noerr_noavg->Draw("integral_min_500>>hNN_intgr_500");
  // amplitude at 290 vs. amplitude of fit
  TH2D * hNN_ampli_fitampli_290 = new TH2D("hNN_ampli_fitampli_290","Amplitude at zero of fit in function of amplitude at 290 of signal",1000,-300,1000,1000,0,1300);
  fit_noerr_noavg->Draw("zero_value_fit:ampli_t0_290>>hNN_ampli_fitampli_290","","colz");
  // integral of fit vs. amplitude of fit
  TH2D * hNN_fitampli_integral_300 = new TH2D("hNN_fitampli_integral_300","Amplitude at zero of fit in function of the integral of the fit (300 bins)",500,0,220000,500,0,1600);
  fit_noerr_noavg->Draw("zero_value_fit:integral_min_300>>hNN_fitampli_integral_300","","colz");
  TH2D * hNN_fitampli_integral_400 = new TH2D("hNN_fitampli_integral_400","Amplitude at zero of fit in function of the integral of the fit (400 bins)",500,0,220000,500,0,1600);
  fit_noerr_noavg->Draw("zero_value_fit:integral_min_400>>hNN_fitampli_integral_400","","colz");
  TH2D * hNN_fitampli_integral_500 = new TH2D("hNN_fitampli_integral_500","Amplitude at zero of fit in function of the integral of the fit (500 bins)",500,0,220000,500,0,1600);
  fit_noerr_noavg->Draw("zero_value_fit:integral_min_500>>hNN_fitampli_integral_500","","colz");
  // integral of fit vs. amplitude of fit WHEN ampli_t0_290 doesn't match integral (sum over bins)
  TH2D * hNN_fitampli_integral_300_badmatch = new TH2D("hNN_fitampli_integral_300_badmatch","Amplitude at zero of fit in function of the integral of the fit (300 bins); amplitude of signal doesn't match integral (sum over bins)",500,0,220000,500,0,1600);
  fit_noerr_noavg->Draw("zero_value_fit:integral_min_300>>hNN_fitampli_integral_300_badmatch","ampli_t0_290<50&&integral_bins_400>5000","colz");
  TH2D * hNN_fitampli_integral_400_badmatch = new TH2D("hNN_fitampli_integral_400_badmatch","Amplitude at zero of fit in function of the integral of the fit (400 bins); amplitude of signal doesn't match integral (sum over bins)",500,0,220000,500,0,1600);
  fit_noerr_noavg->Draw("zero_value_fit:integral_min_400>>hNN_fitampli_integral_400_badmatch","ampli_t0_290<50&&integral_bins_400>5000","colz");
  TH2D * hNN_fitampli_integral_500_badmatch = new TH2D("hNN_fitampli_integral_500_badmatch","Amplitude at zero of fit in function of the integral of the fit (500 bins); amplitude of signal doesn't match integral (sum over bins)",500,0,220000,500,0,1600);
  fit_noerr_noavg->Draw("zero_value_fit:integral_min_500>>hNN_fitampli_integral_500_badmatch","ampli_t0_290<50&&integral_bins_400>5000","colz");
  // amplitude at 290 vs. integral of fit
  TH2D * hNN_ampli_integral_300 = new TH2D("hNN_ampli_integral_300","Amplitude at 290 of signal in function of the integral of the fit (300 bins)",500,0,220000,500,0,1600);
  fit_noerr_noavg->Draw("ampli_t0_290:integral_min_300>>hNN_ampli_integral_300","","colz");
  TH2D * hNN_ampli_integral_400 = new TH2D("hNN_ampli_integral_400","Amplitude at 290 of signal in function of the integral of the fit (400 bins)",500,0,220000,500,0,1600);
  fit_noerr_noavg->Draw("ampli_t0_290:integral_min_400>>hNN_ampli_integral_400","","colz");
  TH2D * hNN_ampli_integral_500 = new TH2D("hNN_ampli_integral_500","Amplitude at 290 of signal in function of the integral of the fit (500 bins)",500,0,220000,500,0,1600);
  fit_noerr_noavg->Draw("ampli_t0_290:integral_min_500>>hNN_ampli_integral_500","","colz");


  
  

  // write on file
  hNN_A->Write();
  hNN_A_pho->Write(); hNN_A_bkg->Write(); hNN_A_good->Write(); hNN_A_bad->Write();
  hNN_A_PG->Write();  hNN_A_PB->Write();  hNN_A_BG->Write();   hNN_A_BB->Write();
  hNN_k->Write();
  hNN_k_pho->Write(); hNN_k_bkg->Write(); hNN_k_good->Write(); hNN_k_bad->Write();
  hNN_k_PG->Write();  hNN_k_PB->Write();  hNN_k_BG->Write();   hNN_k_BB->Write();
  hNN_kA->Write();
  hNN_kA_pho->Write(); hNN_kA_bkg->Write(); hNN_kA_good->Write(); hNN_kA_bad->Write(); 
  hNN_kA_PG->Write();  hNN_kA_PB->Write();  hNN_kA_BG->Write();   hNN_kA_BB->Write();
  hNN_chi2->Write();
  hNN_intgr_300->Write(); hNN_intgr_400->Write(); hNN_intgr_500->Write();
  hNN_ampli_fitampli_290->Write();
  hNN_fitampli_integral_300->Write(); hNN_fitampli_integral_400->Write(); hNN_fitampli_integral_500->Write();
  hNN_fitampli_integral_300_badmatch->Write(); hNN_fitampli_integral_400_badmatch->Write(); hNN_fitampli_integral_500_badmatch->Write();
  hNN_ampli_integral_300->Write(); hNN_ampli_integral_400->Write(); hNN_ampli_integral_500->Write();



  //  ERRORS, NO AVERAGE
  f_data->mkdir("err_noavg");
  f_data->cd("err_noavg");
  // A
  TH1D * hEN_A = new TH1D("hEN_A","Fit errors, no-average: A",100,-40,40);
  fit_err_noavg->Draw("fit_A>>hEN_A","fit_A>-40&&fit_A<40");
  TH1D * hEN_A_pho = new TH1D("hEN_A_pho","Fit errors, no-average: A; photon events",100,0,10);
  fit_err_noavg->Draw("fit_A>>hEN_A_pho","fit_A>-40&&fit_A<40 && ampli_t0_290 > 25");
  TH1D * hEN_A_bkg = new TH1D("hEN_A_bkg","Fit errors, no-average: A; background events",100,-40,40);
  fit_err_noavg->Draw("fit_A>>hEN_A_bkg","fit_A>-40&&fit_A<40 && ampli_t0_290 < 25");
  TH1D * hEN_A_good = new TH1D("hEN_A_good","Fit errors, no-average: A; good fit",100,-40,40);
  fit_err_noavg->Draw("fit_A>>hEN_A_good","fit_A>-40&&fit_A<40 && fit_stat == 0");
  TH1D * hEN_A_bad = new TH1D("hEN_A_bad","Fit errors, no-average: A; bad fit",100,-40,40);
  fit_err_noavg->Draw("fit_A>>hEN_A_bad","fit_A>-40&&fit_A<40 && fit_stat == 4");
  TH1D * hEN_A_PG = new TH1D("hEN_A_PG","Fit errors, no-average: A; photon events, good fit",100,0,10);
  fit_err_noavg->Draw("fit_A>>hEN_A_PG","fit_A>-40&&fit_A<40 && ampli_t0_290 > 25 && fit_stat == 0");
  TH1D * hEN_A_PB = new TH1D("hEN_A_PB","Fit errors, no-average: A; photon events, bad fit",100,-40,40);
  fit_err_noavg->Draw("fit_A>>hEN_A_PB","fit_A>-40&&fit_A<40 && ampli_t0_290 > 25 && fit_stat == 4");
  TH1D * hEN_A_BG = new TH1D("hEN_A_BG","Fit errors, no-average: A; background events, good fit",100,-40,40);
  fit_err_noavg->Draw("fit_A>>hEN_A_BG","fit_A>-40&&fit_A<40 && ampli_t0_290 < 25 && fit_stat == 0");
  TH1D * hEN_A_BB = new TH1D("hEN_A_BB","Fit errors, no-average: A; background events, bad fit",100,-40,40);
  fit_err_noavg->Draw("fit_A>>hEN_A_BB","fit_A>-40&&fit_A<40 && ampli_t0_290 < 25 && fit_stat == 4");
  // k
  TH1D * hEN_k = new TH1D("hEN_k","Fit errors, no-average: k",100,-0.05,0.2);
  fit_err_noavg->Draw("fit_k>>hEN_k","fit_k>-0.05&&fit_k<0.2");
  TH1D * hEN_k_pho = new TH1D("hEN_k_pho","Fit errors, no-average: k; photon events",100,-0.05,0.2);
  fit_err_noavg->Draw("fit_k>>hEN_k_pho","fit_k>-0.05&&fit_k<0.2 && ampli_t0_290 > 25");
  TH1D * hEN_k_bkg = new TH1D("hEN_k_bkg","Fit errors, no-average: k; background events",100,-0.05,0.2);
  fit_err_noavg->Draw("fit_k>>hEN_k_bkg","fit_k>-0.05&&fit_k<0.2 && ampli_t0_290 < 25");
  TH1D * hEN_k_good = new TH1D("hEN_k_good","Fit errors, no-average: k; good fit",100,-0.05,0.2);
  fit_err_noavg->Draw("fit_k>>hEN_k_good","fit_k>-0.05&&fit_k<0.2 && fit_stat == 0");
  TH1D * hEN_k_bad = new TH1D("hEN_k_bad","Fit errors, no-average: k; bad fit",100,-0.05,0.2);
  fit_err_noavg->Draw("fit_k>>hEN_k_bad","fit_k>-0.05&&fit_k<0.2 && fit_stat == 4");
  TH1D * hEN_k_PG = new TH1D("hEN_k_PG","Fit errors, no-average: k; photon events, good fit",100,-0.05,0.2);
  fit_err_noavg->Draw("fit_k>>hEN_k_PG","fit_k>-0.05&&fit_k<0.2 && ampli_t0_290 > 25 && fit_stat == 0");
  TH1D * hEN_k_PB = new TH1D("hEN_k_PB","Fit errors, no-average: k; photon events, bad fit",100,-0.05,0.2);
  fit_err_noavg->Draw("fit_k>>hEN_k_PB","fit_k>-0.05&&fit_k<0.2 && ampli_t0_290 > 25 && fit_stat == 4");
  TH1D * hEN_k_BG = new TH1D("hEN_k_BG","Fit errors, no-average: k; background events, good fit",100,-0.05,0.2);
  fit_err_noavg->Draw("fit_k>>hEN_k_BG","fit_k>-0.05&&fit_k<0.2 && ampli_t0_290 < 25 && fit_stat == 0");
  TH1D * hEN_k_BB = new TH1D("hEN_k_BB","Fit errors, no-average: k; background events, bad fit",100,-0.05,0.2);
  fit_err_noavg->Draw("fit_k>>hEN_k_BB","fit_k>-0.05&&fit_k<0.2 && ampli_t0_290 < 25 && fit_stat == 4");
  // k(A)
  TH2D * hEN_kA = new TH2D("hEN_kA","Fit errors, no-average: k(A)",1000,-40,40,1000,-0.05,0.2);
  fit_err_noavg->Draw("fit_k:fit_A>>hEN_kA","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40","colz");
  TH2D * hEN_kA_pho = new TH2D("hEN_kA_pho", "Fit errors, no-average: k(A); photon events",1000,-40,40,1000,-0.05,0.2);
  fit_err_noavg->Draw("fit_k:fit_A>>hEN_kA_pho","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && ampli_t0_290 > 25","colz");
  TH2D * hEN_kA_bkg = new TH2D("hEN_kA_bkg", "Fit errors, no-average: k(A); background events",1000,-40,40,1000,-0.05,0.2);
  fit_err_noavg->Draw("fit_k:fit_A>>hEN_kA_bkg","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && ampli_t0_290 < 25","colz");
  TH2D * hEN_kA_good = new TH2D("hEN_kA_good", "Fit errors, no-average: k(A); good fit",1000,-40,40,1000,-0.05,0.2);
  fit_err_noavg->Draw("fit_k:fit_A>>hEN_kA_good","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && fit_stat == 0","colz");
  TH2D * hEN_kA_bad = new TH2D("hEN_kA_bad","Fit errors, no-average: k(A); bad  fit",1000,-40,40,1000,-0.05,0.2);
  fit_err_noavg->Draw("fit_k:fit_A>>hEN_kA_bad","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && fit_stat == 4","colz");
  TH2D * hEN_kA_PG = new TH2D("hEN_kA_PG","Fit errors, no-average: k(A); photon events, good fit",1000,-40,40,1000,-0.05,0.2);
  fit_err_noavg->Draw("fit_k:fit_A>>hEN_kA_PG","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && ampli_t0_290 > 25&& fit_stat == 0","colz");
  TH2D * hEN_kA_PB = new TH2D("hEN_kA_PB","Fit errors, no-average: k(A); photon events, bad  fit",1000,-40,40,1000,-0.05,0.2);
  fit_err_noavg->Draw("fit_k:fit_A>>hEN_kA_PB","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && ampli_t0_290 > 25&& fit_stat == 4","colz");
  TH2D * hEN_kA_BG = new TH2D("hEN_kA_BG","Fit errors, no-average: k(A); background events, good fit",1000,-40,40,1000,-0.05,0.2);
  fit_err_noavg->Draw("fit_k:fit_A>>hEN_kA_BG","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && ampli_t0_290 < 25 && fit_stat == 0","colz");
  TH2D * hEN_kA_BB = new TH2D("hEN_kA_BB","Fit errors, no-average: k(A); background events, bad  fit",1000,-40,40,1000,-0.05,0.2);
  fit_err_noavg->Draw("fit_k:fit_A>>hEN_kA_BB","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && ampli_t0_290 < 25 && fit_stat == 4","colz");
  // chi2
  TH1D * hEN_chi2 = new TH1D("hEN_chi2","Fit errors, no-average: chi2",1000,0,20000);
  fit_err_noavg->Draw("fit_chi2>>hEN_chi2","fit_chi2<300000");
  // integrals
  TH1D * hEN_intgr_300 = new TH1D("hEN_intgr_300","Fit errors, no-average: integral 300 bins",1000,0,100000);
  fit_err_noavg->Draw("integral_min_300>>hEN_intgr_300");
  TH1D * hEN_intgr_400 = new TH1D("hEN_intgr_400","Fit errors, no-average: integral 400 bins",1000,0,100000);
  fit_err_noavg->Draw("integral_min_400>>hEN_intgr_400");
  TH1D * hEN_intgr_500 = new TH1D("hEN_intgr_500","Fit errors, no-average: integral 500 bins",1000,0,100000);
  fit_err_noavg->Draw("integral_min_500>>hEN_intgr_500");
  // amplitude at 290 vs. amplitude of fit
  TH2D * hEN_ampli_fitampli_290 = new TH2D("hEN_ampli_fitampli_290","Amplitude at zero of fit in function of amplitude at 290 of signal",1000,-300,1000,1000,0,1300);
  fit_err_noavg->Draw("zero_value_fit:ampli_t0_290>>hEN_ampli_fitampli_290","","colz");
  // integral of fit vs. amplitude of fit
  TH2D * hEN_fitampli_integral_300 = new TH2D("hEN_fitampli_integral_300","Amplitude at zero of fit in function of the integral of the fit (300 bins)",500,0,220000,500,0,1600);
  fit_err_noavg->Draw("zero_value_fit:integral_min_300>>hEN_fitampli_integral_300","","colz");
  TH2D * hEN_fitampli_integral_400 = new TH2D("hEN_fitampli_integral_400","Amplitude at zero of fit in function of the integral of the fit (400 bins)",500,0,220000,500,0,1600);
  fit_err_noavg->Draw("zero_value_fit:integral_min_400>>hEN_fitampli_integral_400","","colz");
  TH2D * hEN_fitampli_integral_500 = new TH2D("hEN_fitampli_integral_500","Amplitude at zero of fit in function of the integral of the fit (500 bins)",500,0,220000,500,0,1600);
  fit_err_noavg->Draw("zero_value_fit:integral_min_500>>hEN_fitampli_integral_500","","colz");
  // integral of fit vs. amplitude of fit WHEN ampli_t0_290 doesn't match integral (sum over bins)
  TH2D * hEN_fitampli_integral_300_badmatch = new TH2D("hEN_fitampli_integral_300_badmatch","Amplitude at zero of fit in function of the integral of the fit (300 bins); amplitude of signal doesn't match integral (sum over bins)",500,0,220000,500,0,1600);
  fit_err_noavg->Draw("zero_value_fit:integral_min_300>>hEN_fitampli_integral_300_badmatch","ampli_t0_290<50&&integral_bins_400>5000","colz");
  TH2D * hEN_fitampli_integral_400_badmatch = new TH2D("hEN_fitampli_integral_400_badmatch","Amplitude at zero of fit in function of the integral of the fit (400 bins); amplitude of signal doesn't match integral (sum over bins)",500,0,220000,500,0,1600);
  fit_err_noavg->Draw("zero_value_fit:integral_min_400>>hEN_fitampli_integral_400_badmatch","ampli_t0_290<50&&integral_bins_400>5000","colz");
  TH2D * hEN_fitampli_integral_500_badmatch = new TH2D("hEN_fitampli_integral_500_badmatch","Amplitude at zero of fit in function of the integral of the fit (500 bins); amplitude of signal doesn't match integral (sum over bins)",500,0,220000,500,0,1600);
  fit_err_noavg->Draw("zero_value_fit:integral_min_500>>hEN_fitampli_integral_500_badmatch","ampli_t0_290<50&&integral_bins_400>5000","colz");
  // amplitude at 290 vs. integral of fit
  TH2D * hEN_ampli_integral_300 = new TH2D("hEN_ampli_integral_300","Amplitude at 290 of signal in function of the integral of the fit (300 bins)",500,0,220000,500,0,1600);
  fit_err_noavg->Draw("ampli_t0_290:integral_min_300>>hEN_ampli_integral_300","","colz");
  TH2D * hEN_ampli_integral_400 = new TH2D("hEN_ampli_integral_400","Amplitude at 290 of signal in function of the integral of the fit (400 bins)",500,0,220000,500,0,1600);
  fit_err_noavg->Draw("ampli_t0_290:integral_min_400>>hEN_ampli_integral_400","","colz");
  TH2D * hEN_ampli_integral_500 = new TH2D("hEN_ampli_integral_500","Amplitude at 290 of signal in function of the integral of the fit (500 bins)",500,0,220000,500,0,1600);
  fit_err_noavg->Draw("ampli_t0_290:integral_min_500>>hEN_ampli_integral_500","","colz");
  // write on file
  hEN_A->Write(); 
  hEN_A_pho->Write(); hEN_A_bkg->Write(); hEN_A_good->Write(); hEN_A_bad->Write();
  hEN_A_PG->Write();  hEN_A_PB->Write();  hEN_A_BG->Write();   hEN_A_BB->Write(); hEN_k->Write();
  hEN_k_pho->Write(); hEN_k_bkg->Write(); hEN_k_good->Write(); hEN_k_bad->Write();
  hEN_k_PG->Write();  hEN_k_PB->Write();  hEN_k_BG->Write();   hEN_k_BB->Write();
  hEN_kA->Write();
  hEN_kA_pho->Write(); hEN_kA_bkg->Write(); hEN_kA_good->Write(); hEN_kA_bad->Write(); 
  hEN_kA_PG->Write();  hEN_kA_PB->Write();  hEN_kA_BG->Write();   hEN_kA_BB->Write();
  hEN_chi2->Write();
  hEN_intgr_300->Write(); hEN_intgr_400->Write(); hEN_intgr_500->Write();
  hEN_ampli_fitampli_290->Write();
  hEN_fitampli_integral_300->Write(); hEN_fitampli_integral_400->Write(); hEN_fitampli_integral_500->Write();
  hEN_fitampli_integral_300_badmatch->Write(); hEN_fitampli_integral_400_badmatch->Write(); hEN_fitampli_integral_500_badmatch->Write();
  hEN_ampli_integral_300->Write(); hEN_ampli_integral_400->Write(); hEN_ampli_integral_500->Write();





  //  ERRORS, AVERAGE
  f_data->mkdir("err_avg");
  f_data->cd("err_avg");
  // A
  TH1D * hEA_A = new TH1D("hEA_A","Fit errors, average: A",1000,-40,40);
  fit_err_avg->Draw("fit_A>>hEA_A","fit_A>-40&&fit_A<40");
  TH1D * hEA_A_pho = new TH1D("hEA_A_pho","Fit errors, average: A; photon events",1000,0,10);
  fit_err_avg->Draw("fit_A>>hEA_A_pho","fit_A>-40&&fit_A<40 && ampli_t0_290 > 25");
  TH1D * hEA_A_bkg = new TH1D("hEA_A_bkg","Fit errors, average: A; background events",1000,-40,40);
  fit_err_avg->Draw("fit_A>>hEA_A_bkg","fit_A>-40&&fit_A<40 && ampli_t0_290 < 25");
  TH1D * hEA_A_good = new TH1D("hEA_A_good","Fit errors, average: A; good fit",1000,-40,40);
  fit_err_avg->Draw("fit_A>>hEA_A_good","fit_A>-40&&fit_A<40 && fit_stat == 0");
  TH1D * hEA_A_bad = new TH1D("hEA_A_bad","Fit errors, average: A; bad fit",1000,-40,40);
  fit_err_avg->Draw("fit_A>>hEA_A_bad","fit_A>-40&&fit_A<40 && fit_stat == 4");
  TH1D * hEA_A_PG = new TH1D("hEA_A_PG","Fit errors, average: A; photon events, good fit",1000,0,10);
  fit_err_avg->Draw("fit_A>>hEA_A_PG","fit_A>-40&&fit_A<40 && ampli_t0_290 > 25 && fit_stat == 0");
  TH1D * hEA_A_PB = new TH1D("hEA_A_PB","Fit errors, average: A; photon events, bad fit",1000,-40,40);
  fit_err_avg->Draw("fit_A>>hEA_A_PB","fit_A>-40&&fit_A<40 && ampli_t0_290 > 25 && fit_stat == 4");
  TH1D * hEA_A_BG = new TH1D("hEA_A_BG","Fit errors, average: A; background events, good fit",1000,-40,40);
  fit_err_avg->Draw("fit_A>>hEA_A_BG","fit_A>-40&&fit_A<40 && ampli_t0_290 < 25 && fit_stat == 0");
  TH1D * hEA_A_BB = new TH1D("hEA_A_BB","Fit errors, average: A; background events, bad fit",1000,-40,40);
  fit_err_avg->Draw("fit_A>>hEA_A_BB","fit_A>-40&&fit_A<40 && ampli_t0_290 < 25 && fit_stat == 4");
  // k
  TH1D * hEA_k = new TH1D("hEA_k","Fit errors, average: k",100,-0.05,0.2);
  fit_err_avg->Draw("fit_k>>hEA_k","fit_k>-0.05&&fit_k<0.2");
  TH1D * hEA_k_pho = new TH1D("hEA_k_pho","Fit errors, average: k; photon events",1000,-0.05,0.2);
  fit_err_avg->Draw("fit_k>>hEA_k_pho","fit_k>-0.05&&fit_k<0.2 && ampli_t0_290 > 25");
  TH1D * hEA_k_bkg = new TH1D("hEA_k_bkg","Fit errors, average: k; background events",1000,-0.05,0.2);
  fit_err_avg->Draw("fit_k>>hEA_k_bkg","fit_k>-0.05&&fit_k<0.2 && ampli_t0_290 < 25");
  TH1D * hEA_k_good = new TH1D("hEA_k_good","Fit errorsi, average: k; good fit",1000,-0.05,0.2);
  fit_err_avg->Draw("fit_k>>hEA_k_good","fit_k>-0.05&&fit_k<0.2 && fit_stat == 0");
  TH1D * hEA_k_bad = new TH1D("hEA_k_bad","Fit errors, average: k; bad fit",1000,-0.05,0.2);
  fit_err_avg->Draw("fit_k>>hEA_k_bad","fit_k>-0.05&&fit_k<0.2 && fit_stat == 4");
  TH1D * hEA_k_PG = new TH1D("hEA_k_PG","Fit errors, average: k; photon events, good fit",1000,-0.05,0.2);
  fit_err_avg->Draw("fit_k>>hEA_k_PG","fit_k>-0.05&&fit_k<0.2 && ampli_t0_290 > 25 && fit_stat == 0");
  TH1D * hEA_k_PB = new TH1D("hEA_k_PB","Fit errors, average: k; photon events, bad fit",1000,-0.05,0.2);
  fit_err_avg->Draw("fit_k>>hEA_k_PB","fit_k>-0.05&&fit_k<0.2 && ampli_t0_290 > 25 && fit_stat == 4");
  TH1D * hEA_k_BG = new TH1D("hEA_k_BG","Fit errors, average: k; background events, good fit",1000,-0.05,0.2);
  fit_err_avg->Draw("fit_k>>hEA_k_BG","fit_k>-0.05&&fit_k<0.2 && ampli_t0_290 < 25 && fit_stat == 0");
  TH1D * hEA_k_BB = new TH1D("hEA_k_BB","Fit errors, average: k; background events, bad fit",1000,-0.05,0.2);
  fit_err_avg->Draw("fit_k>>hEA_k_BB","fit_k>-0.05&&fit_k<0.2 && ampli_t0_290 < 25 && fit_stat == 4");
  // k(A)
  TH2D * hEA_kA = new TH2D("hEA_kA","Fit errors, average: k(A)",1000,-40,40,1000,-0.05,0.2);
  fit_err_avg->Draw("fit_k:fit_A>>hEA_kA","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40","colz");
  TH2D * hEA_kA_pho = new TH2D("hEA_kA_pho", "Fit errors, average: k(A); photon events",1000,-40,40,1000,-0.05,0.2);
  fit_err_avg->Draw("fit_k:fit_A>>hEA_kA_pho","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && ampli_t0_290 > 25","colz");
  TH2D * hEA_kA_bkg = new TH2D("hEA_kA_bkg", "Fit errors, average: k(A); background events",1000,-40,40,1000,-0.05,0.2);
  fit_err_avg->Draw("fit_k:fit_A>>hEA_kA_bkg","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && ampli_t0_290 < 25","colz");
  TH2D * hEA_kA_good = new TH2D("hEA_kA_good", "Fit errors, average: k(A); good fit",1000,-40,40,1000,-0.05,0.2);
  fit_err_avg->Draw("fit_k:fit_A>>hEA_kA_good","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && fit_stat == 0","colz");
  TH2D * hEA_kA_bad = new TH2D("hEA_kA_bad","Fit errors, average: k(A); bad  fit",1000,-40,40,1000,-0.05,0.2);
  fit_err_avg->Draw("fit_k:fit_A>>hEA_kA_bad","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && fit_stat == 4","colz");
  TH2D * hEA_kA_PG = new TH2D("hEA_kA_PG","Fit errors, average: k(A); photon events, good fit",1000,-40,40,1000,-0.05,0.2);
  fit_err_avg->Draw("fit_k:fit_A>>hEA_kA_PG","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && ampli_t0_290 > 25&& fit_stat == 0","colz");
  TH2D * hEA_kA_PB = new TH2D("hEA_kA_PB","Fit errors, average: k(A); photon events, bad  fit",1000,-40,40,1000,-0.05,0.2);
  fit_err_avg->Draw("fit_k:fit_A>>hEA_kA_PB","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && ampli_t0_290 > 25&& fit_stat == 4","colz");
  TH2D * hEA_kA_BG = new TH2D("hEA_kA_BG","Fit errors, average: k(A); background events, good fit",1000,-40,40,1000,-0.05,0.2);
  fit_err_avg->Draw("fit_k:fit_A>>hEA_kA_BG","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && ampli_t0_290 < 25 && fit_stat == 0","colz");
  TH2D * hEA_kA_BB = new TH2D("hEA_kA_BB","Fit errors, average: k(A); background events, bad  fit",1000,-40,40,1000,-0.05,0.2);
  fit_err_avg->Draw("fit_k:fit_A>>hEA_kA_BB","fit_k>-0.05&&fit_k<0.2 && fit_A>-40&&fit_A<40 && ampli_t0_290 < 25 && fit_stat == 4","colz");
  // chi2
  TH1D * hEA_chi2 = new TH1D("hEA_chi2","Fit errors, average: chi2",1000,0,20000);
  fit_err_avg->Draw("fit_chi2>>hEA_chi2","fit_chi2<300000");
  // integrals
  TH1D * hEA_intgr_300 = new TH1D("hEA_intgr_300","Fit errors, average: integral 300 bins",1000,0,100000);
  fit_err_avg->Draw("integral_min_300>>hEA_intgr_300");
  TH1D * hEA_intgr_400 = new TH1D("hEA_intgr_400","Fit errors, average: integral 400 bins",1000,0,100000);
  fit_err_avg->Draw("integral_min_400>>hEA_intgr_400");
  TH1D * hEA_intgr_500 = new TH1D("hEA_intgr_500","Fit errors, average: integral 500 bins",1000,0,100000);
  fit_err_avg->Draw("integral_min_500>>hEA_intgr_500");
  // amplitude at 290 vs. amplitude of fit
  TH2D * hEA_ampli_fitampli_290 = new TH2D("hEA_ampli_fitampli_290","Amplitude at zero of fit in function of amplitude at 290 of signal",1000,-300,1000,1000,0,1300);
  fit_err_avg->Draw("zero_value_fit:ampli_t0_290>>hEA_ampli_fitampli_290","","colz");
  // integral of fit vs. amplitude of fit
  TH2D * hEA_fitampli_integral_300 = new TH2D("hEA_fitampli_integral_300","Amplitude at zero of fit in function of the integral of the fit (300 bins)",500,0,220000,500,0,1600);
  fit_err_avg->Draw("zero_value_fit:integral_min_300>>hEA_fitampli_integral_300","","colz");
  TH2D * hEA_fitampli_integral_400 = new TH2D("hEA_fitampli_integral_400","Amplitude at zero of fit in function of the integral of the fit (400 bins)",500,0,220000,500,0,1600);
  fit_err_avg->Draw("zero_value_fit:integral_min_400>>hEA_fitampli_integral_400","","colz");
  TH2D * hEA_fitampli_integral_500 = new TH2D("hEA_fitampli_integral_500","Amplitude at zero of fit in function of the integral of the fit (500 bins)",500,0,220000,500,0,1600);
  fit_err_avg->Draw("zero_value_fit:integral_min_500>>hEA_fitampli_integral_500","","colz");
  // integral of fit vs. amplitude of fit WHEA ampli_t0_290 doesn't match integral (sum over bins)
  TH2D * hEA_fitampli_integral_300_badmatch = new TH2D("hEA_fitampli_integral_300_badmatch","Amplitude at zero of fit in function of the integral of the fit (300 bins); amplitude of signal doesn't match integral (sum over bins)",500,0,220000,500,0,1600);
  fit_err_avg->Draw("zero_value_fit:integral_min_300>>hEA_fitampli_integral_300_badmatch","ampli_t0_290<50&&integral_bins_400>5000","colz");
  TH2D * hEA_fitampli_integral_400_badmatch = new TH2D("hEA_fitampli_integral_400_badmatch","Amplitude at zero of fit in function of the integral of the fit (400 bins); amplitude of signal doesn't match integral (sum over bins)",500,0,220000,500,0,1600);
  fit_err_avg->Draw("zero_value_fit:integral_min_400>>hEA_fitampli_integral_400_badmatch","ampli_t0_290<50&&integral_bins_400>5000","colz");
  TH2D * hEA_fitampli_integral_500_badmatch = new TH2D("hEA_fitampli_integral_500_badmatch","Amplitude at zero of fit in function of the integral of the fit (500 bins); amplitude of signal doesn't match integral (sum over bins)",500,0,220000,500,0,1600);
  fit_err_avg->Draw("zero_value_fit:integral_min_500>>hEA_fitampli_integral_500_badmatch","ampli_t0_290<50&&integral_bins_400>5000","colz");
  // amplitude at 290 vs. integral of fit
  TH2D * hEA_ampli_integral_300 = new TH2D("hEA_ampli_integral_300","Amplitude at 290 of signal in function of the integral of the fit (300 bins)",500,0,220000,500,0,1600);
  fit_err_avg->Draw("ampli_t0_290:integral_min_300>>hEA_ampli_integral_300","","colz");
  TH2D * hEA_ampli_integral_400 = new TH2D("hEA_ampli_integral_400","Amplitude at 290 of signal in function of the integral of the fit (400 bins)",500,0,220000,500,0,1600);
  fit_err_avg->Draw("ampli_t0_290:integral_min_400>>hEA_ampli_integral_400","","colz");
  TH2D * hEA_ampli_integral_500 = new TH2D("hEA_ampli_integral_500","Amplitude at 290 of signal in function of the integral of the fit (500 bins)",500,0,220000,500,0,1600);
  fit_err_avg->Draw("ampli_t0_290:integral_min_500>>hEA_ampli_integral_500","","colz");
  // compare integrals
  TH2D * hEA_integral_bins_integral_400 = new TH2D("hEA_integral_bins_integral_400","Integral over bins vs. integral of the fitting function (400 bins)",1000,0,220000,1000,0,220000);
  fit_err_avg->Draw("integral_bins_400:integral_min_400>>hEA_integral_bins_integral_400","","colz");
  // compare fitF(0) vs. Integral_bins
  TH2D * hEA_fitampli_integral_bins_400 = new TH2D("hEA_fitampli_integral_bins_400","Amplitude at zero of fit in function of the integral oover 400 bins)",500,0,220000,500,0,1600);
  fit_err_avg->Draw("zero_value_fit:integral_bins_400>>hEA_fitampli_integral_bins_400","","colz");
  // compare minimum amplitude
  TH2D * h_Mamp_ampli = new TH2D("h_Mamp_ampli","Maximum amplitude (amplitude@290)",2000,0,1000,2000,0,1000);
  fit_err_avg->Draw("min:ampli_t0_290>>h_Mamp_ampli","","colz");
  TH2D * h_Mamp_fitampli = new TH2D("h_Mamp_fitampli","Maximum amplitude (amplitude at zero of fit)",2000,0,1000,2000,0,1000);
  fit_err_avg->Draw("min:zero_value_fit>>h_Mamp_fitampli","","colz");
  TH2D * h_Mamp_integral_bins = new TH2D("h_Mamp_integral_bins","Maximum amplitude (integral sum over bins)",2000,0,100000,2000,0,1000);
  fit_err_avg->Draw("min:integral_bins_400>>h_Mamp_integral_bins","","colz");
  TH2D * h_Mamp_integral_400 = new TH2D("h_Mamp_integral_400","Maximum amplitude (integral of fit, 400 bins)",2000,0,100000,2000,0,1000);
  fit_err_avg->Draw("min:integral_min_400>>h_Mamp_integral_400","","colz");
  // "diagonalized" amplitude vs. integral
  TH1D * h_diag_fitampli_intgr_400 = new TH1D("h_diag_fitampli_intgr_400","'diagonalized' amplitude vs. integral of fit",1000,0,3000);
  fit_err_avg->Draw("zero_value_fit+(integral_min_400/100.)>>h_diag_fitampli_intgr_400"); // 45 degrees rotation

  // write on file
  hEA_A->Write(); 
  hEA_A_pho->Write(); hEA_A_bkg->Write(); hEA_A_good->Write(); hEA_A_bad->Write();
  hEA_A_PG->Write();  hEA_A_PB->Write();  hEA_A_BG->Write();   hEA_A_BB->Write(); hEA_k->Write();
  hEA_k_pho->Write(); hEA_k_bkg->Write(); hEA_k_good->Write(); hEA_k_bad->Write();
  hEA_k_PG->Write();  hEA_k_PB->Write();  hEA_k_BG->Write();   hEA_k_BB->Write();
  hEA_kA->Write();
  hEA_kA_pho->Write(); hEA_kA_bkg->Write(); hEA_kA_good->Write(); hEA_kA_bad->Write(); 
  hEA_kA_PG->Write();  hEA_kA_PB->Write();  hEA_kA_BG->Write();   hEA_kA_BB->Write();
  hEA_chi2->Write();
  hEA_intgr_300->Write(); hEA_intgr_400->Write(); hEA_intgr_500->Write();
  hEA_ampli_fitampli_290->Write();
  hEA_fitampli_integral_300->Write(); hEA_fitampli_integral_400->Write(); hEA_fitampli_integral_500->Write();
  hEA_fitampli_integral_300_badmatch->Write(); hEA_fitampli_integral_400_badmatch->Write(); hEA_fitampli_integral_500_badmatch->Write();
  hEA_ampli_integral_300->Write(); hEA_ampli_integral_400->Write(); hEA_ampli_integral_500->Write();
  hEA_integral_bins_integral_400->Write();
  hEA_fitampli_integral_bins_400->Write();
  h_Mamp_ampli->Write(); h_Mamp_fitampli->Write(); h_Mamp_integral_bins->Write(); h_Mamp_integral_400->Write();
  h_diag_fitampli_intgr_400->Write();

  f_data->cd();
  f_data->Close();
  delete f_data;

  return 0;

}























