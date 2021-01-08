//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// PlotTestSpectra
// Erin Conley (erin.conley@duke.edu)
// Description: Superimpose different test spectra contained in the
//              input/test_spectra directory
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include "TH1D.h"
#include "TMatrixDBase.h"
#include "TH2.h"
#include "../headers/ForwardFit.h"

void PlotTestSpectra(){
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // HERE ARE THE VARIABLES THAT CHANGE
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //define the test spectra we would like to superimpose
	std::vector<TString> ts; std::vector<TString> names;
	//example: different cross section models
	ts.emplace_back("MARLEY_DCR_nueOnly_PNXscn_StepEffic"); names.emplace_back("MARLEY (p, n)");
	ts.emplace_back("MARLEY_DCR_nueOnly_SNOwGLoBESXscn_StepEffic"); names.emplace_back("SNOwGLoBES");
	ts.emplace_back("MARLEY_DCR_nueOnly_40TiXscn_StepEffic"); names.emplace_back("40-Ti xscn");
	
    //define supernova distance
    Double_t distanceSN =10.0; //kpc
    
    TString supertitle("Test spectra for different cross section models");
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE INPUTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //mass factor for DUNE statistics
    Double_t massfact = 1.0; //so that we see the full DUNE far detector response

    //convert to string for output file
    TString distanceStr; distanceStr.Form("%.2lf", distanceSN);   
    
    //get histograms
    std::vector<TH1D*> hists;
    for(size_t i = 0; i < ts.size(); ++i){
		TString tmp = "input/test_spectra/pinched_test_smeared_sum_" + ts[i] + ".dat";
		
		TString ts_name = "test_spectrum_" + ts[i];
		
		TH1D *test_spectrum_tmp = fill_spect_hist(tmp, ts_name, distanceSN, massfact);
		//TH1D *test_spectrum_tmp = fill_spect_hist_NoRebin(tmp, ts_name, distanceSN, massfact);
		
		hists.emplace_back(test_spectrum_tmp);
		
	}
    
    //use same style settings that Kate used
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1,0);
    
    TString title = supertitle + ";Observed Energy (MeV);Events per 0.5 MeV";
	
    superimposeTH1D(hists, names, title);
    
}
