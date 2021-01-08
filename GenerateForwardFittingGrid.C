//GenerateForwardFittingGrid.C
//Erin Conley (erin.conley@duke.edu)
//
//Creates two data files:
//1) A "pinched info" file that contains (alpha, e0, lum) information for every element in the grid
//2) A "grid info" file that contains (alpha, e0, lum) minimum/maximum/step information for the grid
//
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

void GenerateForwardFittingGrid()
{
	//define dat file title
	TString outfile("pinched_info_2019October17.dat");
	TString gridfile("grid_info_2019October17.dat");
	
	//assume all three flavors get the same ranges
	Double_t first_alpha = 0.1;
	Double_t last_alpha = 7.0;
	
	Double_t step_alpha = 0.1;
	
	Double_t first_e0 = 5.0;
	Double_t last_e0 = 20.0;

	Double_t step_e0 = 0.1;
	
	Double_t factor = 1e53;
	
    Double_t factorNuxLum = 0;//1;
	
    Double_t factorNuebarE = 0;//12.0/9.5;
	
	//upper bounds
	Double_t first_lum = 0.2*factor;
	Double_t last_lum = 1.0*factor;
	Double_t step_lum = 0.025*factor;
	
	//open gridfile
	ofstream osgrid;
	osgrid.open(gridfile);
	osgrid << first_alpha << " " << last_alpha << " " << step_alpha << " ";
	osgrid << first_e0 << " " << last_e0 << " " << step_e0 << " ";
	osgrid << first_lum << " " << last_lum << " " << step_lum << std::endl;
	osgrid.close();
	
	//open outfile
	ofstream os;
	os.open(outfile);	 
	Int_t counter = 0;
	
	for(Double_t iAlpha = first_alpha; iAlpha < last_alpha; iAlpha += step_alpha){
	for(Double_t iNueE = first_e0; iNueE < last_e0+step_e0; iNueE += step_e0){
	for(Double_t iLum = first_lum; iLum < last_lum+step_lum; iLum += step_lum){
		
		double nuebarE = (int)( std::round((iNueE*factorNuebarE)*10.0) )/10.0;
		double nuxE = (int)( std::round((nuebarE*(1.3))*10.0) )/10.0;
		double nuxLum = iLum*factorNuxLum;
		
		os << counter << " " << iAlpha << " " << iAlpha << " " << iAlpha << " ";
		os << iNueE << " " << nuebarE << " " << nuxE << " ";
		os << iLum << " " << iLum << " " << nuxLum << std::endl;
		
		
		++counter;
	}}}	
	
	
	os.close();
	
	std::cout << "Number of elements: " << counter << std::endl;
    
}//end of code
