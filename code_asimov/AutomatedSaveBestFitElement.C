//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// AutomatedSaveBestFitElement.C
// Erin Conley (erin.conley@duke.edu)
// Description: Save best fit spectrum for several combinations of grid/test
//              spectrum parameters, e.g., cross section models
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include "TH1D.h"
#include "TMatrixDBase.h"
#include "TH2.h"
#include "headers/ForwardFit.h"
#include "headers/Asimov.h"

void AutomatedSaveBestFitElement(){
    //define supernova distance
	Double_t distanceSN =10.0; //kpc
    
    //define the smearing matrix we want to use
	TString grid_smear("MARLEY_DCR_nueOnly"); 
        
    //define the grid we want to use
	TString grid("2019October17");
	
    //define the test spectra we want to use
    TString ts_smear("MARLEY_DCR_nueOnly");
     
    //define the xscn models we want to use
    std::vector<TString> grid_xscn;
    grid_xscn.emplace_back("PNXscn");
    grid_xscn.emplace_back("40TiXscn");
    grid_xscn.emplace_back("SNOwGLoBESXscn");
    
    //define the test spectra xscn models we want to use
    std::vector<TString> ts_xscn;
    ts_xscn.emplace_back("PNXscn");
    ts_xscn.emplace_back("40TiXscn");
    ts_xscn.emplace_back("SNOwGLoBESXscn");
     
     //define the efficiency model(s) we want to use
     TString grid_effic("StepEffic");
     TString ts_effic("StepEffic");
     
     for(int i = 0; i < grid_xscn.size(); ++i){
        for(int j = 0; j < ts_xscn.size(); ++j){
			SaveBestFitElement(distanceSN, grid, grid_smear, grid_xscn[i], grid_effic, ts_smear, ts_xscn[j], ts_effic);
        }
    }
    
    
    
    
}//AutomatedSaveBestFitElement
