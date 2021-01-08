//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ResolutionVsDistanceStudy_AutomatedAsimov
// Erin Conley (erin.conley@duke.edu)
// Description: Produces chi2-minimized flux parameter measurements and plots
//              for different supernova distance and smearing resolution
//              (contained in the SNOwGLoBES smearing used to produce the grid
//              and test spectra).
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

void ResolutionVsDistanceStudy_AutomatedAsimov(){
    //define supernova distance
    std::vector<TString> dists;  std::vector<double> distance;
    dists.emplace_back("1.00kpc"); distance.emplace_back(1.0);
    dists.emplace_back("2.00kpc"); distance.emplace_back(2.0);
    dists.emplace_back("3.00kpc"); distance.emplace_back(3.0);
    dists.emplace_back("4.00kpc"); distance.emplace_back(4.0);
    dists.emplace_back("5.00kpc"); distance.emplace_back(5.0);
    dists.emplace_back("6.00kpc"); distance.emplace_back(6.0);
    dists.emplace_back("7.00kpc"); distance.emplace_back(7.0);
    dists.emplace_back("8.00kpc"); distance.emplace_back(8.0);
    dists.emplace_back("9.00kpc"); distance.emplace_back(9.0);
    dists.emplace_back("10.00kpc"); distance.emplace_back(10.0);
    
    //define the grid we want to use
    TString grid("2019October17");
    
	//define the xscn model we want to use 
	//same for grid and test spectra for this study
	TString xscn("PNXscn");    
	
	//define the efficiency we want to use
	TString effic("StepEffic");
	
	//define the grid smearing we want to use
	std::vector<TString> grid_smear;
	grid_smear.emplace_back("0.00MARLEY_NO_nueOnly");
	grid_smear.emplace_back("0.05MARLEY_NO_nueOnly");
	grid_smear.emplace_back("0.10MARLEY_NO_nueOnly");
	grid_smear.emplace_back("0.15MARLEY_NO_nueOnly");
	grid_smear.emplace_back("0.20MARLEY_NO_nueOnly");
	grid_smear.emplace_back("0.25MARLEY_NO_nueOnly");
	grid_smear.emplace_back("0.30MARLEY_NO_nueOnly");
	    
	for(int i = 0; i < grid_smear.size(); ++i){
		for(int j = 0; j < distance.size(); ++j){
			
			AsimovMethod(distance[j], grid, grid_smear[i], xscn, effic, grid_smear[i], xscn, effic);
			
		}
	}
    
       
}//end ResolutionVsDistanceStudy_AutomatedAsimov
