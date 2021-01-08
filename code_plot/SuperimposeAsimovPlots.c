//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// SuperimposeAsimovPlots
// Erin Conley (erin.conley@duke.edu)
// Description: Creates a TCanvas with superimposed contour plots produced by
//              the Asimov method. Specify which plots to superimpose in the
//              "config" vector.
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

void SuperimposeAsimovPlots(){
	
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // HERE ARE THE VARIABLES WE CHANGE
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    TString distance("10.00");
    
    std::vector<TString> config; std::vector<TString> names;
	
	//all config parameters will live in the same TString; names are the legend labels
	//config will look something like this: smearSMEARING_XSCN_EFFIC_spectraSMEARING_XSCN_EFFIC
    //Note that the grid ("smear") and test spectra ("spectra") configuration could be different
    //Example below: different mass ordering assumptions
    config.emplace_back( "smearMARLEY_DCR_nueOnly_PNXscn_StepEffic_spectraMARLEY_DCR_nueOnly_PNXscn_StepEffic"); names.emplace_back("No mass ordering assumptions");
    config.emplace_back( "smearMARLEY_DCR_NO_nueOnly_PNXscn_StepEffic_spectraMARLEY_DCR_NO_nueOnly_PNXscn_StepEffic"); names.emplace_back("Normal mass ordering assumptions");
    config.emplace_back( "smearMARLEY_DCR_IO_nueOnly_PNXscn_StepEffic_spectraMARLEY_DCR_IO_nueOnly_PNXscn_StepEffic"); names.emplace_back("Inverted mass ordering assumptions");
    
	
    std::vector<TString> histNames; //histograms in the ROOT files to find
	histNames.emplace_back("alphavse0_90cut");
	histNames.emplace_back("lumvse0_90cut");
	histNames.emplace_back("lumvsalpha_90cut");
	
	//define the grid we want to use
    TString grid("2019October17");
    
    //define extra directory where the files live, e.g., "/CrossSectionStudy/"
    //ps: if there isn't an extra directory and they live in "out", set the parameter to "/"
    TString dir("/");
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE INPUTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //define file names   
    std::vector<TString> files;
	for(size_t i = 0; i < config.size(); ++i){
		TString tmp = "out" + dir + "chi2plots_" + grid + "_" + config[i] + "_" + distance + "kpc.root";
		files.emplace_back(tmp);
	}
	
	TString gridinfo = "input/grid_info/grid_info_" + grid + ".dat";
    TString pdfout = "out/SuperimposedSensitivityRegions_" + grid + "_" + distance + "kpc_AsmiovMethod.pdf";
    
    //set specific aesthetic things here
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1,0);
    
    //need to first define the truth points
    Double_t alpha_true = 2.5;
    Double_t e0_true = 9.5; //mev
    Double_t lum_true = 5e52; //ergs
	//define contours
	double contours[1];
	contours[0] = 0.001;
 
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE GRID BOUNDS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //define the bounds of the grid
	Double_t first_alpha;
	Double_t last_alpha;
	Double_t step_alpha;
	
	Double_t first_e0;
	Double_t last_e0;
	Double_t step_e0; 
	
	Double_t factor = 1e53;
	Double_t first_lum;
	Double_t last_lum;
	Double_t step_lum;	 
	
	//open grid_info.dat
    ifstream gin;
    gin.open(gridinfo);
    
    Int_t ig = 0;
    while(1){
        gin >> first_alpha >> last_alpha >> step_alpha >> first_e0 >> last_e0 >> step_e0 >> first_lum >> last_lum >> step_lum;
        if(!gin.good()) break;        
        
        ++ig;
    }
    
    gin.close();
    
    Double_t minalpha = first_alpha - step_alpha/2.0;
    Double_t maxalpha = last_alpha + step_alpha/2.0;
    Double_t mine0 = first_e0 - step_e0/2.0;
    Double_t maxe0 = last_e0 + step_e0/2.0;
    Double_t minlum = first_lum - step_lum/2.0;
    Double_t maxlum = last_lum + step_lum/2.0;

	Int_t numalphabins = int(last_alpha - first_alpha)/step_alpha+1;
	Int_t nume0bins = int(last_e0 - first_e0)/step_e0+1;
	Int_t numlumbins = (last_lum - first_lum)/(step_lum) + 1;
    
    //std::cout << "The flux grid follows this definition:" << std::endl;
    //std::cout << "Alpha: [" << first_alpha << ", " << last_alpha << "] with " << step_alpha << " spacing" << std::endl;
    //std::cout << "E0: [" << first_e0 << ", " << last_e0 << "] with " << step_e0 << " spacing" << std::endl;
    //std::cout << "Luminosity: [" << first_lum << ", " << last_lum << "] with " << step_lum << " spacing" << std::endl;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE TRUTH POINTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //define truth point of test spectrum
    Double_t OGalpha[1];
    Double_t OGe0[1];
    Double_t OGlum[1];
    
    OGalpha[0] = alpha_true;
    OGe0[0] = e0_true;
    OGlum[0] = lum_true/factor; //10^{53} erg
    
    //make TGraph object to draw
    TGraph *OGpoint_AlphaVsE0 = new TGraph(1, OGe0, OGalpha);
    OGpoint_AlphaVsE0->SetMarkerStyle(29);
    OGpoint_AlphaVsE0->SetMarkerColor(28);
    OGpoint_AlphaVsE0->SetMarkerSize(2);
    
    TGraph *OGpoint_LumVsE0 = new TGraph(1, OGe0, OGlum);
    OGpoint_LumVsE0->SetMarkerStyle(29);
    OGpoint_LumVsE0->SetMarkerColor(28);
    OGpoint_LumVsE0->SetMarkerSize(2);
    
    TGraph *OGpoint_LumVsAlpha = new TGraph(1, OGalpha, OGlum);
    OGpoint_LumVsAlpha->SetMarkerStyle(29);
    OGpoint_LumVsAlpha->SetMarkerColor(28);
    OGpoint_LumVsAlpha->SetMarkerSize(2);
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GRAB EVENT INFORMATION
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TString a("alpha");
    TString e("e0");
    TString l("lum");
    
    std::vector<std::vector<double>> alpha = getManyParameters(files, a);
    std::vector<std::vector<double>> e0 = getManyParameters(files, e);
    std::vector<std::vector<double>> lum = getManyParameters(files, l);  
    
	std::vector<TGraph*> BestFitPoints_AlphaVsE0;
	std::vector<TGraph*> BestFitPoints_LumVsE0;
	std::vector<TGraph*> BestFitPoints_LumVsAlpha;
	std::vector<TString> BestFitPoints_labels;
	
	for(size_t i = 0; i < alpha.size(); ++i){
		Double_t tmpalpha[1];
		Double_t tmpe0[1];
		Double_t tmplum[1];
		
		tmpalpha[0] = alpha[i][0];
		tmpe0[0] = e0[i][0];
		tmplum[0] = lum[i][0];
		
		TString label = "Best-Fit: #alpha = " + TString::Format("%4.1f",tmpalpha[0]) + ", #LT E_{#nu} #GT = " + TString::Format("%4.1f",tmpe0[0]) + " MeV, #epsilon = " + TString::Format("%4.1f",tmplum[0]) + "e53 ergs";
		
	    TGraph *tmppoint_AlphaVsE0 = new TGraph(1, tmpe0, tmpalpha);
		tmppoint_AlphaVsE0->SetMarkerStyle(29);
		tmppoint_AlphaVsE0->SetMarkerColor(kBlack);
		tmppoint_AlphaVsE0->SetMarkerSize(1.5);
    
		TGraph *tmppoint_LumVsE0 = new TGraph(1, tmpe0, tmplum);
		tmppoint_LumVsE0->SetMarkerStyle(29);
		tmppoint_LumVsE0->SetMarkerColor(kBlack);
		tmppoint_LumVsE0->SetMarkerSize(1.5);
		
		TGraph *tmppoint_LumVsAlpha = new TGraph(1, tmpalpha, tmplum);
		tmppoint_LumVsAlpha->SetMarkerStyle(29);
		tmppoint_LumVsAlpha->SetMarkerColor(kBlack);
		tmppoint_LumVsAlpha->SetMarkerSize(1.5);	
		
		if(i+1==5){ 
			tmppoint_AlphaVsE0->SetMarkerColor(i+2);
			tmppoint_LumVsE0->SetMarkerColor(i+2);
			tmppoint_LumVsAlpha->SetMarkerColor(i+2);
		}
		else if(i+1==3){ 
			tmppoint_AlphaVsE0->SetMarkerColor(8);
			tmppoint_LumVsE0->SetMarkerColor(8);
			tmppoint_LumVsAlpha->SetMarkerColor(8);
		}
		else{ 
			tmppoint_AlphaVsE0->SetMarkerColor(i+1);
			tmppoint_LumVsE0->SetMarkerColor(i+1);
			tmppoint_LumVsAlpha->SetMarkerColor(i+1);
		}
		
		BestFitPoints_AlphaVsE0.emplace_back(tmppoint_AlphaVsE0);
		BestFitPoints_LumVsE0.emplace_back(tmppoint_LumVsE0);
		BestFitPoints_LumVsAlpha.emplace_back(tmppoint_LumVsAlpha);
		BestFitPoints_labels.emplace_back(label);
	}
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GRAB CHI2 PLOTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
	std::vector<std::vector<TH2D*>> hists;
	for(int iHist = 0; iHist < histNames.size(); ++iHist){  
		std::vector<TH2D*> tmp;
		for(int iFile = 0; iFile < files.size(); ++iFile){
			TFile *f = new TFile(files[iFile]);
			TH2D *h = (TH2D*)gDirectory->Get(histNames[iHist]);
			h->SetContour(1, contours);
			h->SetLineWidth(3);
			
			//output contour area 
			//std::cout << "Area for " << files[iFile] << ", " << histNames[iHist] << ": " << findAreaOfContour(h) << std::endl; 
			
			//also try to make the labels/axis titles better
			h->GetXaxis()->CenterTitle();
			h->GetXaxis()->SetTitleOffset(1.2);
			h->GetXaxis()->SetTitleSize(0.04);
			h->GetXaxis()->SetLabelSize(0.05);
			h->GetYaxis()->CenterTitle();
			h->GetYaxis()->SetTitleSize(0.04);
			h->GetYaxis()->SetLabelSize(0.05);
						
			tmp.emplace_back(h);
			//f->Close();
			
		}
		hists.emplace_back(tmp);
	}	
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DRAW PLOTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
	TCanvas *c1 = new TCanvas("c1", "c1", 1000, 700);
    c1->Divide(2, 2);
    
	TLegend *leg = makeLegend(hists[0], names);
	leg->AddEntry(OGpoint_AlphaVsE0, "Truth Point", "p");
	for(size_t i = 0; i < BestFitPoints_AlphaVsE0.size(); ++i)
		leg->AddEntry(BestFitPoints_AlphaVsE0[i], BestFitPoints_labels[i], "p");
	
	c1->cd(2);
	leg->Draw();
	
    for(int i = 0; i < hists.size(); ++i){
		
		int toDraw = i+1;
		if(i >= 1) ++toDraw;
		c1->cd(toDraw);
		
		for(int j = 0; j < hists[i].size(); ++j){
			if(j+1==5) hists[i][j]->SetLineColor(j+2);
			else if(j+1==3) hists[i][j]->SetLineColor(8);
			else hists[i][j]->SetLineColor(j+1);
			
			if(j==0) hists[i][j]->Draw("cont3");
			else hists[i][j]->Draw("cont3 same");
			
		}
		
	}
	
	//also draw the truth points 
	c1->cd(1); 
	OGpoint_AlphaVsE0->Draw("P SAME");
	for(size_t i = 0; i < BestFitPoints_AlphaVsE0.size(); ++i) BestFitPoints_AlphaVsE0[i]->Draw("P SAME");
	c1->cd(3); 
	OGpoint_LumVsE0->Draw("P SAME");
	for(size_t i = 0; i < BestFitPoints_LumVsE0.size(); ++i) BestFitPoints_LumVsE0[i]->Draw("P SAME");
	c1->cd(4); 
	OGpoint_LumVsAlpha->Draw("P SAME");   
	for(size_t i = 0; i < BestFitPoints_LumVsAlpha.size(); ++i) BestFitPoints_LumVsAlpha[i]->Draw("P SAME");
	
	c1->Print(pdfout);
    
}//end of code
