//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// PlotBestFitElement
// Erin Conley (erin.conley@duke.edu)
// Description: Superimpose test spectrum and best-fit spectrum for given test
//              spectrum and grid parameter assumptions.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include "TH1D.h"
#include "TMatrixDBase.h"
#include "TH2.h"
#include "../headers/ForwardFit.h"
#include "../headers/Asimov.h"

void PlotBestFitElement(){
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // HERE ARE THE VARIABLES THAT CHANGE
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //define supernova distance
	Double_t distanceSN =10.0; //kpc
    
    //define the smearing matrix we want to use
	TString grid_smear("MARLEY_DCR_nueOnly");
        
    //define the grid we want to use
    TString grid("2019October17");
     
    //define the test spectra we want to use
    TString ts_smear("MARLEY_DCR_nueOnly");
    
    //define the xscn model(s) we want to use
    TString grid_xscn("PNXscn");
    TString gridlabel = "Grid " + grid_xscn;
    TString ts_xscn("PNXscn");
	TString tslabel = "Test Spectrum " + ts_xscn;
     
    //define the efficiency model(s) we want to use
    TString grid_effic("StepEffic");
    TString ts_effic("StepEffic");
    
    //true grid element
    Double_t true_alpha = 2.5;
    Double_t true_e0 = 9.5;
    Double_t true_lum = 0.5;
    
    std::cout << "Producing best-fit plot for a " << distanceSN << "kpc supernova and the following input parameters:" << std::endl;
    std::cout << "	Grid: " << grid << " with smearing " << grid_smear << ", xscn " << grid_xscn << ", and effic " << grid_effic << std::endl;  
    std::cout << "	Test spectrum: " << ts_smear << ", xscn " << ts_xscn << " and effic " << ts_effic << std::endl;    
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE INPUTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //define string to hold grid parameters
    TString grid_pars = "smear_" + grid_smear + "_" + grid_xscn + "_" + grid_effic;
    TString ts_pars = ts_smear + "_" + ts_xscn + "_" + ts_effic;
	TString test_spect_filename = "input/test_spectra/pinched_test_smeared_sum_" + ts_pars + ".dat";
        
    TString smeareddir = "input/" + grid + "/" + grid_pars;
    TString pinchedinfo = "input/pinched_info/pinched_info_" + grid + ".dat";
    TString gridinfo = "input/grid_info/grid_info_" + grid + ".dat";
    //define output file
    TString distanceStr; distanceStr.Form("%.2lf", distanceSN);
    
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

    //Range of physical parameters
    Double_t mingoodalpha = first_alpha;//1.0;
    Double_t maxgoodalpha = last_alpha;//7.0;
    Double_t mingoode0 = first_e0;//0.0;
    Double_t maxgoode0 = last_e0;//20.;
    
    Int_t numalphabins = int(last_alpha - first_alpha)/step_alpha+1;
    Int_t nume0bins = int(last_e0 - first_e0)/step_e0+1;
    Int_t numlumbins = (last_lum - first_lum)/(step_lum) + 1;
    
    std::cout << "The " << grid << " grid follows this definition:" << std::endl;
    std::cout << "    Alpha: [" << first_alpha << ", " << last_alpha << "] with " << step_alpha << " spacing" << std::endl;
    std::cout << "    E0: [" << first_e0 << ", " << last_e0 << "] with " << step_e0 << " spacing" << std::endl;
    std::cout << "    Luminosity: [" << first_lum << ", " << last_lum << "] with " << step_lum << " spacing" << std::endl;

    std::cout << "The parameter-fitting algorithm will consider the following physical range for alpha and E0:" << std::endl;
    std::cout << "    Alpha: [" << mingoodalpha << ", " << maxgoodalpha << "]" << std::endl;
    std::cout << "    E0: [" << mingoode0 << ", " << maxgoode0 << "]" << std::endl;
    
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
    // MAKE HISTOS AND PLOT
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    std::vector<TH1D*> histos; std::vector<TString> histos_names;
    
    TH1D* test_spectrum = fill_spect_hist(test_spect_filename,"test_spectrum",distanceSN, massfact);
    
    //now we loop over templates in the grid
    Int_t iFlux; //iterator over the grid
    Int_t numfluxes = ip; //total number of fluxes we care about
    
    const Int_t maxhist = maxflux;
    
    Double_t chi2min = 100000000000.0;
    Int_t dofbest;
    Int_t ifbest;
    
    Double_t chi2;
    Int_t dof;
    
    for(iFlux = 0; iFlux < numfluxes; ++iFlux){
        
        // Do not go outside the limits
        if(anue[iFlux] < maxgoodalpha && bnue[iFlux] < maxgoode0 && anue[iFlux] > mingoodalpha && bnue[iFlux] > mingoode0 ){
            
            //define filename
            TString pargridfilename = smeareddir+"/pinched_"+TString::Format("%d",iFlux)+"_smeared_sum.dat";
            TString histname = "pinched_"+grid + "_" + grid_smear + "_" + grid_xscn + "_" + grid_effic + "_" + ts_smear + "_" + ts_xscn + "_" + ts_effic + "_" + distanceStr + "kpc_"+TString::Format("%d",iFlux);
            
            TH1D *pargridhist = fill_spect_hist(pargridfilename, histname, distanceSN, massfact);
            mychi2(test_spectrum,pargridhist,&chi2,&dof);
                        
            
            if (chi2<chi2min) {
                chi2min = chi2;
                dofbest = dof;
                ifbest = iFlux;
            }
            
        }
    }

	std::cout << "Best-fit measurement: (" << anue[ifbest] << ", " << bnue[ifbest] << ", " << cnue[ifbest] << ") with reduced chi2 = " << chi2min/dofbest << std::endl;    
    
    
	histos.emplace_back(test_spectrum);
	histos_names.emplace_back(tslabel);
    
    TString bflabel = gridlabel + " best-fit: (" + TString::Format("%.1lf", anue[ifbest]) + ", " + TString::Format("%.1lf", bnue[ifbest]) + ", " + TString::Format("%.2lf", cnue[ifbest]/factor) + "e53)";
	TString bffilename = smeareddir+"/pinched_"+TString::Format("%d",ifbest)+"_smeared_sum.dat";
	TString bfname = "best_fit_" + TString::Format("%.1lf", anue[ifbest]) + "_" + TString::Format("%.1lf", bnue[ifbest]) + "_" + TString::Format("%.2lf", cnue[ifbest]/factor) + "e53";
	TH1D *bfhist = fill_spect_hist(bffilename, bfname, distanceSN, massfact);
	
	histos.emplace_back(bfhist);
	histos_names.emplace_back(bflabel);
		
	superimposeTH1D(histos, histos_names, "Test Spectrum with its Best Fit;Observed Energy (MeV);Events per 5 MeV");
    
}//end PlotBestFitElement
