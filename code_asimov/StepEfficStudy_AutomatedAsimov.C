//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// StepEfficStudy_AutomatedAsimov
// Erin Conley (erin.conley@duke.edu)
// Description: Produces chi2-minimized flux parameter measurements and plots
//              for different combinations of detection threshold values
//              (defined by the SNOwGLoBES post-smearing efficiency used to
//              produce the grid and test spectra).
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

void StepEfficStudy_AutomatedAsimov(){
    //define supernova distance
    Double_t distanceSN =10.0; //kpc
    
    //define the smearing matrix we want to use
    TString smear("MARLEY_DCR_NO_nueOnly");
	
	//define the xscn model we want to use
	TString xscn("PNXscn");
    
    //define the grid we want to use
    TString grid("2019October17");
    
    //define the efficiency model we want to use
    TString effic("StepEffic");
    
    //define the detection threshold for the grid
    std::vector<Double_t> grid_effic;
    grid_effic.emplace_back(3.0);
    grid_effic.emplace_back(3.5);
    grid_effic.emplace_back(4.0);
    grid_effic.emplace_back(4.5);
    grid_effic.emplace_back(5.0);
    grid_effic.emplace_back(5.5);
    grid_effic.emplace_back(6.0);
    grid_effic.emplace_back(6.5);
    grid_effic.emplace_back(7.0);
    
    //define the test spectrum detection threshold
    std::vector<Double_t> ts_effic;
    ts_effic.emplace_back(3.0);
    ts_effic.emplace_back(3.5);
    ts_effic.emplace_back(4.0);
    ts_effic.emplace_back(4.5);
    ts_effic.emplace_back(5.0);
    ts_effic.emplace_back(5.5);
    ts_effic.emplace_back(6.0);
    ts_effic.emplace_back(6.5);
    ts_effic.emplace_back(7.0);
    
	for(int i = 0; i < grid_effic.size(); ++i){
		for(int j = 0; j < ts_effic.size(); ++j){
			
			AsimovMethodStepEfficStudy(distanceSN, grid, smear, xscn, effic, grid_effic[i], ts_effic[j]);
			
		}
	}    
    
       
}//end StepEfficStudy_AutomatedAsimov
