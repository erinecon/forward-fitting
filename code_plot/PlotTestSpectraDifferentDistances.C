//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// PlotTestSpectraDifferentDistances
// Erin Conley (erin.conley@duke.edu)
// Description: Superimpose different test spectra contained in the
//              input/test_spectra directory with different SN distances
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include "TH1D.h"
#include "TMatrixDBase.h"
#include "TH2.h"
#include "../headers/ForwardFit.h"

void PlotTestSpectraDifferentDistances(){
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // HERE ARE THE VARIABLES THAT CHANGE
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //define the test spectra we would like to superimpose
	std::vector<Double_t> distances; std::vector<TString> names;
	
	//specific order for tech note plot
	distances.emplace_back(1.0); names.emplace_back("1.0 kpc");
	distances.emplace_back(5.0); names.emplace_back("5.0 kpc");
	distances.emplace_back(10.0); names.emplace_back("10.0 kpc");
	distances.emplace_back(15.0); names.emplace_back("15.0 kpc");
	distances.emplace_back(20.0); names.emplace_back("20.0 kpc");
	
    TString ts_config("MARLEY_DCR_nueOnly_PNXscn_StepEffic");
    
    TString supertitle("Test spectra for different supernova distances");
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE INPUTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //mass factor for DUNE statistics
    Double_t massfact = 1.0; //so that we see the full DUNE far detector response
    
    //get histograms
    std::vector<TH1D*> hists;
    for(size_t i = 0; i < distances.size(); ++i){
		TString tmp = "input/test_spectra/pinched_test_smeared_sum_" + ts_config + ".dat";
		
		TString ts_name = "test_spectrum_" + ts_config + TString::Format("%.2lf",distances[i]);
		
		TH1D *test_spectrum_tmp = fill_spect_hist(tmp, ts_name, distances[i], massfact);
		
		std::cout << ts_name << ": " << test_spectrum_tmp->Integral() << std::endl;
		//std::cout << ts_name << ": " << test_spectrum_tmp->Integral( 1, test_spectrum_tmp->FindBin(20.) ) << std::endl;
		
		hists.emplace_back(test_spectrum_tmp);
		
	}
    
    //use same style settings that Kate used
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1,0);
    
    TString title = supertitle + ";Observed Energy (MeV);Events per 0.5 MeV";
	
    superimposeTH1D(hists, names, title);
    
}
