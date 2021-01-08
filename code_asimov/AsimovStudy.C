//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// AsimovStudy
// Erin Conley (erin.conley@duke.edu)
// Description: Produces chi2-minimized flux parameter measurements and plots
//              for defined grid and test spectrum assumptions.
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

void AsimovStudy(){
    //define supernova distance
    Double_t distanceSN = 10.0; //kpc
    
    //define the smearing matrix we want to use
    //grid == detector assumptions for reconstruction ability
    TString grid_smear("MARLEY_DCR_nueOnly");
        
    //define the grid we want to use
	TString grid("2019October17");
	
     //define the test spectrum we want to use
     //test spectrum == the "real data"
     TString ts("MARLEY_DCR_nueOnly");
	
     //define the xscn model(s) we want to use
     TString grid_xscn("PNXscn");
     
     TString ts_xscn("PNXscn");
     
     //define the efficiency model(s) we want to use
     TString grid_effic("StepEffic");
     
     TString ts_effic("StepEffic");
     
     //perform the study!
     AsimovMethod(Double_t distanceSN, TString grid, TString grid_smear, TString grid_xscn, TString grid_effic, TString ts, TString ts_xscn, TString ts_effic)
      
    
}
