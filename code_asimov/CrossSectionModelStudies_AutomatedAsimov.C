//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// CrossSectionModelStudies_AutomatedAsimov
// Erin Conley (erin.conley@duke.edu)
// Description: Produces chi2-minimized flux parameter measurements and plots
//              for different combinations of nue-Ar40 cross section models
//              (used in SNOwGLoBES to produce the grid and test spectra).
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

void CrossSectionModelStudies_AutomatedAsimov(){
    
    //define supernova distance
    Double_t distanceSN =10.0; //kpc
    
    //define the smearing we want to consider
    TString smear("MARLEY_DCR_NO_nueOnly");
    
    //define the xscn models we want to use
    std::vector<TString> grid_xscn;
    grid_xscn.emplace_back("SNOwGLoBESXscn");
    grid_xscn.emplace_back("PNXscn");
    grid_xscn.emplace_back("40TiXscn");
    
    //define the test spectra xscn models we want to use
    std::vector<TString> ts_xscn;
    ts_xscn.emplace_back("SNOwGLoBESXscn");
    ts_xscn.emplace_back("PNXscn");
    ts_xscn.emplace_back("40TiXscn");
    
    //define the grid we want to use
    TString grid("2019October17");
    
    //define the efficiency we want to use
    TString effic("StepEffic");
    
    for(int i = 0; i < grid_xscn.size(); ++i){
        for(int j = 0; j < ts_xscn.size(); ++j){
            
			AsimovMethod(distanceSN, grid, smear, grid_xscn[i], effic, smear, ts_xscn[j], effic);
			
        }
    }
    
       
}//end of main code
