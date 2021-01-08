#include "TH1D.h"
#include "TMatrixDBase.h"
#include "TH2.h"

#include "../headers/ForwardFit.h"
#include "../headers/Supernovae.h"

void GenerateSupernovaEvents(){
    //define how many events we want to keep
    Int_t numberOfSN = 2000; 
    
    //define supernova distance
    Double_t distanceSN = 10.0; //kpc
    
    //define the grid we want to use
    TString grid("2019October17");
    
    //define the smearing matrix we want to use
	TString grid_smear("MARLEY");
    TString ts("MARLEY");

    //define the xscn model(s) we want to use
    TString grid_xscn("PNXscn");
    TString ts_xscn("PNXscn");
     
    //define the efficiency model(s) we want to use
    TString grid_effic("StepEffic");
    TString ts_effic("StepEffic");
     
    FakeSupernovaMethod(numberOfSN, distanceSN, grid, grid_smear, grid_xscn, grid_effic, ts, ts_xscn, ts_effic);
        
}//end of function
