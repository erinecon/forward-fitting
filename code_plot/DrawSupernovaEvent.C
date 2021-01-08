//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// DrawSupernovaEvent
// Erin Conley (erin.conley@duke.edu)
// Description: Generate and draw energy spectrum with statistical fluctuations
//              (via the "fake supernova" method)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include "TH1D.h"
#include "TMatrixDBase.h"
#include "TH2.h"

#include "../headers/ForwardFit.h"
#include "../headers/Supernovae.h"

void DrawSupernovaEvent(){
    //define supernova distance
    Double_t distanceSN =10.0; //kpc
    
    //define the grid we want to use
    TString grid("2019October17");
    
    //define the smearing we want to use
	TString grid_smear("MARLEY");
    TString ts("MARLEY");
     
     //define the xscn model(s) we want to use
     TString grid_xscn("PNXscn");
     TString ts_xscn("PNXscn");
     
     //define the efficiency model(s) we want to use
     TString grid_effic("StepEffic");
     TString ts_effic("StepEffic");
     
     DrawFakeSupernova(distanceSN, grid, grid_smear, grid_xscn, grid_effic, ts, ts_xscn, ts_effic);
        
}//end of function
