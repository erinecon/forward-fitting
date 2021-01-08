//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// MassOrderingStudies_AutomatedAsimov
// Erin Conley (erin.conley@duke.edu)
// Description: Produces chi2-minimized flux parameter measurements and plots
//              for different neutrino mass ordering assumptions (contained in
//              the SNOwGLoBES smearing used to produce the grid and test
//              spectra).
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

void MassOrderingStudies_AutomatedAsimov(){
    
    //define supernova distance
    Double_t distanceSN =10.0; //kpc
    
    //define the smearing we want to consider
    std::vector<TString> grid_smear;
    grid_smear.emplace_back("MARLEY_DCR_nueOnly");
    grid_smear.emplace_back("MARLEY_DCR_NO_nueOnly");
    grid_smear.emplace_back("MARLEY_DCR_IO_nueOnly");
    
    //test spectra smearing we want to use
    std::vector<TString> ts_smear;
    ts_smear.emplace_back("MARLEY_DCR_nueOnly");
    ts_smear.emplace_back("MARLEY_DCR_NO_nueOnly");
    ts_smear.emplace_back("MARLEY_DCR_IO_nueOnly");
    
    //define the xscn models we want to use
    TString xscn("PNXscn");
    
    //define the grid we want to use
    TString grid("2019October17");
    
    //define the efficiency we want to use
    TString effic("StepEffic");
    
    for(int i = 0; i < grid_smear.size(); ++i){
        for(int j = 0; j < ts_smear.size(); ++j){
			
			AsimovMethod(distanceSN, grid, grid_smear[i], xscn, effic, ts_smear[j], xscn, effic);
			
        }
    }
    
       
}//end of main code
