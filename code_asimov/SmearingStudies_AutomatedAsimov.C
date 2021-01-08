//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// SmearingStudies_AutomatedAsimov.C
// Erin Conley (erin.conley@duke.edu)
// Description: Produces chi2-minimized flux parameter measurements and plots
//              for different combinations of SNOwGLoBES smearing (used to
//              produce the grid and test specta).
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include <string>
#include <iostream>
#include <dirent.h>
#include <vector>
#include <math.h>
#include <algorithm>

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TF1.h"
#include "TMath.h"
#include "TText.h"

#include "../headers/ForwardFit.h"
#include "../headers/Asimov.h"

void SmearingStudies_AutomatedAsimov(){
    //define supernova distance
    Double_t distanceSN =10.0; //kpc
    
    //define the grid we want to use
    TString grid("2019October17");

	//define the xscn model we want to use 
	//same for grid and test spectra for this study
	TString xscn("PNXscn");    
	
	//define the efficiency we want to use
	TString effic("StepEffic");
	
	//define the grid smearing we want to use
	std::vector<TString> grid_smear;
	grid_smear.emplace_back("MARLEY_DCR_NO");
	grid_smear.emplace_back("MARLEY_DCR_NO_nueOnly");
    
    //define the test spectra smearing we want to use 
	std::vector<TString> tss;
	tss.emplace_back("MARLEY_DCR_NO");
	tss.emplace_back("MARLEY_DCR_NO_nueOnly");
	    
	for(int i = 0; i < grid_smear.size(); ++i){
		for(int j = 0; j < tss.size(); ++j){
			
			AsimovMethod(distanceSN, grid, grid_smear[i], xscn, effic, tss[j], xscn, effic);

		}
	}
    
       
}
