//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// PlotFluxGridElement
// Erin Conley (erin.conley@duke.edu)
// Description: For given test spectrum parameter assumptions and grid elements
//              described by their (alpha, e0, luminosity) values, plot the test
//              spectrum and spectra corresponding to the grid elements.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include "TH1D.h"
#include "TMatrixDBase.h"
#include "TH2.h"
#include "../headers/ForwardFit.h"
#include "../headers/Asimov.h"

void PlotFluxGridElement(){
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // HERE ARE THE VARIABLES THAT CHANGE
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //define supernova distance
    Double_t distanceSN = 10.0; //kpc
    
    //define the grid we want to use
    TString grid("2019October17");
    
    //define the smearing matrix we want to use
	TString grid_smear("MARLEY_DCR_nueOnly");
    TString ts_smear("MARLEY_DCR_nueOnly");
     
    //define the xscn model(s) we want to use
    TString grid_xscn("PNXscn");
    TString ts_xscn("PNXscn");
     
    //define the efficiency model(s) we want to use
    TString grid_effic("StepEffic");
    TString ts_effic("StepEffic");
    
    //grid elements to plot
    std::vector<double> values_alpha; std::vector<double> values_e0; std::vector<double> values_lum;
    values_alpha.emplace_back(2.5); values_e0.emplace_back(9.5); values_lum.emplace_back(1.0);
    values_alpha.emplace_back(0.1); values_e0.emplace_back(5.1); values_lum.emplace_back(0.3);
    values_alpha.emplace_back(6.0); values_e0.emplace_back(15.0); values_lum.emplace_back(0.6);
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE INPUTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //grid element we care about
    Double_t true_alpha = 2.5;
    Double_t true_e0 = 9.5;
    Double_t true_lum = 0.5;
    
	TString test_spect_filename = "input/test_spectra/pinched_test_smeared_sum_" + ts_smear + "_" + ts_xscn + "_" + ts_effic + ".dat";
        
    TString smeareddir = "input/" + grid + "/smear_" + grid_smear + "_" + grid_xscn + "_" + grid_effic;
    TString pinchedinfo = "input/pinched_info/pinched_info_" + grid + ".dat";
    TString gridinfo = "input/grid_info/grid_info_" + grid + ".dat";
    //define output file
    TString distanceStr; distanceStr.Form("%.2lf", distanceSN);
    
	//Range of physical parameters
    Double_t mingoodalpha = 1.0;
    Double_t maxgoodalpha = 7.0;
    Double_t mingoode0 = 5.0;
    Double_t maxgoode0 = 20.;
    
    Double_t massfact = 1.0; //see the full DUNE far detector response
    Double_t chi2_90contour = 4.61;  //pdg 2018
    
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1,0);
    //set random seed
    gRandom->SetSeed(0);
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE GRID BOUNDS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //define the bounds of the grid
    Double_t first_alpha;
    Double_t last_alpha;
    Double_t step_alpha;
    
    Double_t first_e0;
    Double_t last_e0;
    Double_t step_e0;
    
    Double_t factor = 1e53;
    Double_t first_lum;
    Double_t last_lum;
    Double_t step_lum;
    
    //open grid_info.dat
    ifstream gin;
    gin.open(gridinfo);
    
    Int_t ig = 0;
    while(1){
        gin >> first_alpha >> last_alpha >> step_alpha >> first_e0 >> last_e0 >> step_e0 >> first_lum >> last_lum >> step_lum;
        if(!gin.good()) break;
        
        ++ig;
    }
    
    gin.close();
    
    Double_t minalpha = first_alpha - step_alpha/2.0;
    Double_t maxalpha = last_alpha + step_alpha/2.0;
    Double_t mine0 = first_e0 - step_e0/2.0;
    Double_t maxe0 = last_e0 + step_e0/2.0;
    Double_t minlum = first_lum - step_lum/2.0;
    Double_t maxlum = last_lum + step_lum/2.0;
    
    Double_t minalpha2 = first_alpha - step_alpha*2.0;
    Double_t maxalpha2 = last_alpha + step_alpha*2.0;
    Double_t mine02 = first_e0 - step_e0*2.0;
    Double_t maxe02 = last_e0 + step_e0*2.0;
    Double_t minlum2 = first_lum - step_lum*2.0;
    Double_t maxlum2 = last_lum + step_lum*2.0;
    
    Int_t numalphabins = int(last_alpha - first_alpha)/step_alpha+1;
    Int_t nume0bins = int(last_e0 - first_e0)/step_e0+1;
    Int_t numlumbins = (last_lum - first_lum)/(step_lum) + 1;
    
    std::cout << "The " << grid << " grid follows this definition:" << std::endl;
    std::cout << "    Alpha: [" << first_alpha << ", " << last_alpha << "] with " << step_alpha << " spacing" << std::endl;
    std::cout << "    E0: [" << first_e0 << ", " << last_e0 << "] with " << step_e0 << " spacing" << std::endl;
    std::cout << "    Luminosity: [" << first_lum << ", " << last_lum << "] with " << step_lum << " spacing" << std::endl;
    
   	
	//for the central + corner plot
	Double_t central_alpha = (first_alpha + last_alpha)/2.0;
	Double_t central_e0 = (first_e0 + last_e0)/2.0;
	Double_t central_lum = (first_lum + last_lum)/2.0;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // READ PINCHING PARAMETERS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //open the pinched parameter file
    ifstream pin;
    pin.open(pinchedinfo);
    
    //define arrays, iterators to keep the parameters
    const Int_t maxflux = 350000;
    
    std::vector<Int_t> pnum; pnum.resize(maxflux);
    std::vector<Double_t> anue; anue.resize(maxflux);
    std::vector<Double_t> bnue; bnue.resize(maxflux);
    std::vector<Double_t> cnue; cnue.resize(maxflux);
    std::vector<Double_t> anuebar; anuebar.resize(maxflux);
    std::vector<Double_t> bnuebar; bnuebar.resize(maxflux);
    std::vector<Double_t> cnuebar; cnuebar.resize(maxflux);
    std::vector<Double_t> anux; anux.resize(maxflux);
    std::vector<Double_t> bnux; bnux.resize(maxflux);
    std::vector<Double_t> cnux; cnux.resize(maxflux);
    
    std::cout << "Reading flux parameters: " << std::endl;
    
    Int_t ip = 0;
    while(1){
        
        //Int_t p1;
        //Double_t a1, a2, a3, b1, b2, b3, c1, c2, c3;
        //pin >> p1 >> a1 >> a2 >> a3 >> b1 >> b2 >> b3 >> c1 >> c2 >> c3;
        
        pin >> pnum[ip] >> anue[ip] >> anuebar[ip] >> anux[ip] >> bnue[ip] >> bnuebar[ip] >> bnux[ip] >> cnue[ip] >> cnuebar[ip] >> cnux[ip];
        if(!pin.good()) break;
        
        //std::cout << pnum[ip]<<" "<<anue[ip]<<" "<<anuebar[ip]<<" "<<anux[ip]<<" "<<bnue[ip]<<" "<<bnuebar[ip]<<" "<<bnux[ip]<<" "<<cnue[ip]<<" "<<cnuebar[ip]<<" "<<cnux[ip]<<std::endl;
        
        
        ++ip;
        
    }
    
    pin.close();
    
    std::cout << "This grid contains " << ip << " fluxes." << std::endl;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE HISTOS AND PLOT
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    TH1D* test_spectrum = fill_spect_hist(test_spect_filename,"test_spectrum",distanceSN, massfact);
    
    //now we loop over templates in the grid
	Int_t iFlux; //iterator over the grid
	Int_t numfluxes = ip; //total number of fluxes we care about
	        
	const Int_t maxhist = maxflux;
    
    Double_t chi2;
    Int_t dof;
	
	std::vector<TH1D*> histos;
	std::vector<TGraph*> graphs;
	std::vector<TString> histos_names;
	
	TString tslabel = "Test Spectrum: (" + TString::Format("%.1lf", true_alpha) + ", " + TString::Format("%.1lf", true_e0) + ", " + TString::Format("%.2lf", true_lum) + "e53)";
	histos.emplace_back(test_spectrum);
	histos_names.emplace_back(tslabel);
	
	
	for(size_t i = 0; i < values_alpha.size(); ++i){
		for(iFlux = 0; iFlux < numfluxes; ++iFlux){
			double diff1 = ( values_alpha[i] - anue[iFlux] )*100.0;
			double diff2 = ( values_e0[i] - bnue[iFlux] )*100.0;
			double diff3 = ( values_lum[i] - cnue[iFlux]/factor )*100.0;

			if( (anue[iFlux] == values_alpha[i] && bnue[iFlux] == values_e0[i] && cnue[iFlux] == values_lum[i]*factor) || (diff1 == 0.0 && diff2 ==0.0 && diff3 == 0.0)) {
				
				
				//define filename 
				TString pargridfilename = smeareddir+"/pinched_"+TString::Format("%d",iFlux)+"_smeared_sum.dat";
				TString histname = "pinched_"+grid + "_" + grid_smear + "_" + grid_xscn + "_" + grid_effic + "_" + ts_smear + "_" + ts_xscn + "_" + ts_effic + "_" + distanceStr + "kpc_"+TString::Format("%d",iFlux);
				
				//TH1D *pargridhist = fill_spect_hist_NoRebin(pargridfilename, histname, distanceSN, massfact);
				TH1D *pargridhist = fill_spect_hist(pargridfilename, histname, distanceSN, massfact);
				mychi2(test_spectrum,pargridhist,&chi2,&dof);
				
				TGraph *g = makeTGraphFromFilename(pargridfilename, histname, distanceSN, massfact);
				
				//useful label
				TString histlabel;
				histlabel = "#splitline{Grid Element " + TString::Format("%d", iFlux) + ", Red. #chi^{2} = " + TString::Format("%.2lf", chi2/dof) + ":}{(" + TString::Format("%.1lf", anue[iFlux]) + ", " + TString::Format("%.1lf", bnue[iFlux]) + ", " + TString::Format("%.2lf", values_lum[i]) + "e53)}";
				pargridhist->SetLineWidth(0);
				pargridhist->SetMarkerStyle(20);
				
				//std::cout << chi2/dof << std::endl;
				
				histos.emplace_back(pargridhist);
				graphs.emplace_back(g);
				
				histos_names.emplace_back(histlabel);
				
				
			}//end check for the values we want 
		}//end loop over fluxes
		
	}//end loop over values we want
	
	
	superimposeTH1D(histos, histos_names, "Test Spectrum with Example Grid Element;Observed Energy (MeV);Events per 0.5 MeV");
		
    
}//end PlotFluxGridElement
