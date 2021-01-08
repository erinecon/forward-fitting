//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// SupernovaDistanceStudy_AutomatedAsimov
// Erin Conley (erin.conley@duke.edu)
// Description: Produces chi2-minimized flux parameter measurements and plots
//              for different supernova distance uncertainty factors.
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

void SupernovaDistanceStudy_AutomatedAsimov(){
	double trueSNdistance = 10.0;
	
	std::vector<double> SNdistance_uncertainties;
	SNdistance_uncertainties.emplace_back(0.10); //10% uncertainty
	SNdistance_uncertainties.emplace_back(0.20); //20% uncertainty
	SNdistance_uncertainties.emplace_back(0.50); //50% uncertainty
    
    //define the grid we want to use
    TString grid("2019October17");

	//define the xscn model we want to use
	TString xscn("PNXscn");    
	
	//define the efficiency we want to use
	TString effic("StepEffic");
	
	//define the grid smearing we want to use
	TString smear("MARLEY_DCR_NO_nueOnly");
	
	//no uncertainty
	AsimovMethod(trueSNdistance, grid, smear, xscn, effic, smear, xscn, effic);
	
	//incorporate the uncertainties into the measurement
	for(size_t i = 0; i < SNdistance_uncertainties.size(); ++i)
		AsimovMethodSNDistance(trueSNdistance, SNdistance_uncertainties[i], grid, smear, xscn, effic, smear, xscn, effic);
	
	
}//end SupernovaDistanceStudy_AutomatedAsimov
