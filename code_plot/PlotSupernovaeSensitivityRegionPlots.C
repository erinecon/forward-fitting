//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// PlotSupernovaeSensitivityRegionPlots
// Erin Conley (erin.conley@duke.edu)
// Description: Creates a TCanvas with superimposed contour plots produced by
//              the "fake supernova" method.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include <string>
#include <iostream>
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
#include "TRint.h"
#include "TROOT.h"

#include "../headers/ForwardFit.h"
#include "../headers/Supernovae.h"

void PlotSupernovaeSensitivityRegionPlots(){
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // HERE ARE THE VARIABLES WE CHANGE
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //define how many events we want to keep
    TString numEvents("1000");	
    
    //define supernova distance
    TString distance("10.00"); //kpc
    
    //define the smearing matrix we want to use
	TString grid_smear("MARLEY_nueOnly");
        
    //define the grid we want to use
    TString grid("2019October17");
     
     //define the test spectra we want to use
     TString ts("MARLEY_nueOnly");
     
     //define the xscn model(s) we want to use
     TString grid_xscn("PNXscn");
     
     TString ts_xscn("PNXscn");
     
     //define the efficiency model(s) we want to use
     TString grid_effic("StepEffic");
     
     TString ts_effic("StepEffic");
    
    //define extra directory where the files live, e.g., "/CrossSectionStudy/"
    //ps: if there isn't an extra directory and they live in "out", set the parameter to "/"
    TString dir("/");
     
     TString options = "smear" + grid_smear + "_" + grid_xscn + "_" + grid_effic + "_spectra" + ts + "_" + ts_xscn + "_" + ts_effic;
	
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE INPUTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//file to open
	TString gridinfo = "input/grid_info/grid_info_" + grid + ".dat";
    TString filename("out"+dir+"supernovae_" + grid + "_" + options + "_" + distance + "kpc_" + numEvents + "events.root");
    
	//set specific aesthetic things here 
	gStyle->SetOptStat(0);
    gStyle->SetPalette(1,0);
    gStyle->SetNumberContours(99);
    
    //need to first define the truth points
    Double_t alpha_true = 2.5;
    Double_t e0_true = 9.5; //mev
    Double_t lum_true = 5e52; //ergs
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GRAB EVENT INFORMATION
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    std::vector<double> eventid;
    std::vector<double> alpha;
    std::vector<double> e0;
    std::vector<double> lum;
    
    TFile *f = new TFile(filename);
    TTree *tr = (TTree*)gDirectory->Get("data;1");
    TTreeReader read(tr);
    TTreeReaderValue<double> eid(read, "eventid");
    TTreeReaderValue<double> a(read, "alpha");
    TTreeReaderValue<double> e(read, "e0");
    TTreeReaderValue<double> l(read, "lum");
    
    while(read.Next()){
        eventid.emplace_back(*eid);
        alpha.emplace_back(*a);
        e0.emplace_back(*e);
        lum.emplace_back(*l);
    }
	
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
    
    first_alpha = *std::min_element(alpha.begin(), alpha.end())-step_alpha*2.0;
    last_alpha = *std::max_element(alpha.begin(), alpha.end())+step_alpha*2.0;
    
    first_e0 = *std::min_element(e0.begin(), e0.end())-step_e0*2.0;
    last_e0 = *std::max_element(e0.begin(), e0.end())+step_e0*2.0;
    
    first_lum = *std::min_element(lum.begin(), lum.end())*factor-step_lum*2.0;
    last_lum = *std::max_element(lum.begin(), lum.end())*factor+step_lum*2.0;
    
    
    Double_t minalpha = first_alpha ;//- step_alpha/2.0;
    Double_t maxalpha = last_alpha ;//+ step_alpha/2.0;
    Double_t mine0 = first_e0 ;//- step_e0/2.0;
    Double_t maxe0 = last_e0 ;//+ step_e0/2.0;
    Double_t minlum = first_lum ;//- step_lum/2.0;
    Double_t maxlum = last_lum ;//+ step_lum/2.0;

	Int_t numalphabins = int(last_alpha - first_alpha)/step_alpha+1;
	Int_t nume0bins = int(last_e0 - first_e0)/step_e0+1;
	Int_t numlumbins = (last_lum - first_lum)/(step_lum) + 1;
    
    //std::cout << "The flux grid follows this definition:" << std::endl;
    std::cout << "Alpha: [" << first_alpha << ", " << last_alpha << "] with " << step_alpha << " spacing" << std::endl;
    std::cout << "E0: [" << first_e0 << ", " << last_e0 << "] with " << step_e0 << " spacing" << std::endl;
    std::cout << "Luminosity: [" << first_lum << ", " << last_lum << "] with " << step_lum << " spacing" << std::endl;
    
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
    OGpoint_AlphaVsE0->SetMarkerColor(kBlack);
    OGpoint_AlphaVsE0->SetMarkerSize(2);
    
    TGraph *OGpoint_LumVsE0 = new TGraph(1, OGe0, OGlum);
    OGpoint_LumVsE0->SetMarkerStyle(29);
    OGpoint_LumVsE0->SetMarkerColor(kBlack);
    OGpoint_LumVsE0->SetMarkerSize(2);
    
    TGraph *OGpoint_LumVsAlpha = new TGraph(1, OGalpha, OGlum);
    OGpoint_LumVsAlpha->SetMarkerStyle(29);
    OGpoint_LumVsAlpha->SetMarkerColor(kBlack);
    OGpoint_LumVsAlpha->SetMarkerSize(2);
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAKE 2D PLOTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	TH2D *hAllowedRegion_AlphaVsE0 = new TH2D("alphavse0", ";;#alpha", nume0bins, mine0, maxe0, numalphabins, minalpha, maxalpha );
	TH2D *hAllowedRegion_LumVsE0 = new TH2D("lumvse0", ";#LT E_{#nu} #GT (MeV);#epsilon (10^{53} erg)", nume0bins, mine0, maxe0, numlumbins, minlum/factor, maxlum/factor );
	TH2D *hAllowedRegion_LumVsAlpha = new TH2D("lumvsalpha", ";#alpha;", numalphabins, minalpha, maxalpha, numlumbins, minlum/factor, maxlum/factor );

    for(int i = 0; i < eventid.size(); ++i){
		hAllowedRegion_AlphaVsE0->Fill(e0[i], alpha[i]);
		hAllowedRegion_LumVsE0->Fill(e0[i], lum[i]);
        hAllowedRegion_LumVsAlpha->Fill(alpha[i], lum[i]);
        
    }
    
    //set some aesthetic stuff here 
	hAllowedRegion_AlphaVsE0->GetYaxis()->CenterTitle();
    hAllowedRegion_AlphaVsE0->GetYaxis()->SetTitleSize(0.05);
    
    hAllowedRegion_LumVsE0->GetXaxis()->CenterTitle();
    hAllowedRegion_LumVsE0->GetXaxis()->SetTitleSize(0.04);
    hAllowedRegion_LumVsE0->GetYaxis()->CenterTitle();
    hAllowedRegion_LumVsE0->GetYaxis()->SetTitleSize(0.04);
    
    hAllowedRegion_LumVsAlpha->GetXaxis()->CenterTitle();
    hAllowedRegion_LumVsAlpha->GetXaxis()->SetTitleSize(0.05);
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FIND CONTOUR LINES
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Double_t contours_AlphaVsE0[1];
    contours_AlphaVsE0[0] = findContourLine(hAllowedRegion_AlphaVsE0, 0.9);
        
    Double_t contours_LumVsE0[1];
    contours_LumVsE0[0] = findContourLine(hAllowedRegion_LumVsE0, 0.9);
    
    Double_t contours_LumVsAlpha[1];
    contours_LumVsAlpha[0] = findContourLine(hAllowedRegion_LumVsAlpha, 0.9);
        
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DRAW 2D PLOTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    TLegend *leg = new TLegend(0.1, 0.1, 0.9, 0.9);
    leg->SetFillColor(10);
    leg->SetTextFont( (Font_t) 62);
    leg->SetTextSize(0.05);
    leg->SetBorderSize(1);
    leg->AddEntry(hAllowedRegion_AlphaVsE0, "90% Contour", "l");
    leg->AddEntry(OGpoint_AlphaVsE0, "Truth: #alpha = 2.5, #LT E_{#nu} #GT = 9.5 MeV, #varepsilon = 5e52 ergs", "p");
    
    TCanvas *c1 = new TCanvas("c1","c1",1000,700);
    c1->Divide(2, 2);
    c1->cd(1);
    hAllowedRegion_AlphaVsE0->DrawCopy("colz");
    hAllowedRegion_AlphaVsE0->SetContour(1, contours_AlphaVsE0);
    hAllowedRegion_AlphaVsE0->SetLineWidth(5);
    hAllowedRegion_AlphaVsE0->SetLineColor(kBlack);
    //hAllowedRegion_AlphaVsE0->Draw("cont3");
    hAllowedRegion_AlphaVsE0->Draw("cont3 same");
    OGpoint_AlphaVsE0->Draw("P SAME");
    c1->cd(2);
    leg->Draw();
    c1->cd(3);
	hAllowedRegion_LumVsE0->DrawCopy("colz");
	hAllowedRegion_LumVsE0->SetContour(1, contours_LumVsE0);
    hAllowedRegion_LumVsE0->SetLineWidth(5);
    hAllowedRegion_LumVsE0->SetLineColor(kBlack);
    //hAllowedRegion_LumVsE0->Draw("CONT3");
    hAllowedRegion_LumVsE0->Draw("CONT3 same");
    OGpoint_LumVsE0->Draw("P SAME");
    c1->cd(4);
    hAllowedRegion_LumVsAlpha->DrawCopy("colz");
    hAllowedRegion_LumVsAlpha->SetContour(1, contours_LumVsAlpha);
    hAllowedRegion_LumVsAlpha->SetLineWidth(5);
    hAllowedRegion_LumVsAlpha->SetLineColor(kBlack);
    //hAllowedRegion_LumVsAlpha->Draw("CONT3");
    hAllowedRegion_LumVsAlpha->Draw("CONT3 same");
    OGpoint_LumVsAlpha->Draw("P SAME");
    
    
}//end of code
