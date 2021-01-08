//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// SuperimposeSupernovaePlots
// Erin Conley (erin.conley@duke.edu)
// Description: Creates a TCanvas with superimposed contour plots produced by
//              the "fake supernova" method. Specify which plots to superimpose
//              in the "options" vector.
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
#include "../headers/Supernovae.h"

void SuperimposeSupernovaePlots(){
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // HERE ARE THE VARIABLES WE CHANGE
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    TString numEvents("1000");
    
    //define the grid we want to use
    TString grid("2019October17");
    
    std::vector<TString> options; std::vector<TString> names;
    options.emplace_back( "smearMARLEY_NO_nueOnly_PNXscn_StepEffic_spectraMARLEY_NO_nueOnly_PNXscn_StepEffic_10.00kpc"); names.emplace_back("10.00 kpc");
    
    options.emplace_back( "smearMARLEY_NO_nueOnly_PNXscn_StepEffic_spectraMARLEY_NO_nueOnly_PNXscn_StepEffic_7.00kpc"); names.emplace_back("7.00 kpc");
    
    options.emplace_back( "smearMARLEY_NO_nueOnly_PNXscn_StepEffic_spectraMARLEY_NO_nueOnly_PNXscn_StepEffic_5.00kpc"); names.emplace_back("5.00 kpc");
    
    //define extra directory where the files live, e.g., "/CrossSectionStudy/"
    //ps: if there isn't an extra directory and they live in "out", set the parameter to "/"
    TString dir("/");
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE INPUTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //define file names   
    std::vector<TString> files;
    for(size_t i = 0; i < options.size(); ++i){
				
		TString tmp = "out"+dir+"supernovae_" + grid + "_" + options[i] + "_" + numEvents + "events.root";
		files.emplace_back(tmp);
		
	}
    
    TString gridinfo = "input/grid_info/grid_info_" + grid + ".dat";
    
    //set specific aesthetic things here
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1,0);
    
    //need to first define the truth points
    Double_t alpha_true = 2.5;
    Double_t e0_true = 9.5; //mev
    Double_t lum_true = 5e52; //ergs

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
    
    last_alpha = 4.5;
    first_e0 = 7.0;
    last_e0 = 12.0;
    first_lum = 3.5e52;
    last_lum = 6.5e52;
    
    Double_t minalpha = first_alpha - 2.5*step_alpha;
    Double_t maxalpha = last_alpha + 2.5*step_alpha;
    Double_t mine0 = first_e0 - 2.5*step_e0;
    Double_t maxe0 = last_e0 + 2.5*step_e0;
    Double_t minlum = first_lum - 2.5*step_lum;
    Double_t maxlum = last_lum + 2.5*step_lum;

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
    OGlum[0] = lum_true/factor; //*10^{53} erg
    
    //make TGraph object to draw
    TGraph *OGpoint_AlphaVsE0 = new TGraph(1, OGe0, OGalpha);
    OGpoint_AlphaVsE0->SetMarkerStyle(29);
    OGpoint_AlphaVsE0->SetMarkerColor(kBlue);
    OGpoint_AlphaVsE0->SetMarkerSize(2);
    
    TGraph *OGpoint_LumVsE0 = new TGraph(1, OGe0, OGlum);
    OGpoint_LumVsE0->SetMarkerStyle(29);
    OGpoint_LumVsE0->SetMarkerColor(kBlue);
    OGpoint_LumVsE0->SetMarkerSize(2);
    
    TGraph *OGpoint_LumVsAlpha = new TGraph(1, OGalpha, OGlum);
    OGpoint_LumVsAlpha->SetMarkerStyle(29);
    OGpoint_LumVsAlpha->SetMarkerColor(kBlue);
    OGpoint_LumVsAlpha->SetMarkerSize(2);
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GRAB EVENT INFORMATION
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TString eid("eventid");
    TString a("alpha");
    TString e("e0");
    TString l("lum");
    
    std::vector<std::vector<double>> eventid = getManyParameters(files, eid);
    std::vector<std::vector<double>> alpha = getManyParameters(files, a);
    std::vector<std::vector<double>> e0 = getManyParameters(files, e);
    std::vector<std::vector<double>> lum = getManyParameters(files, l);
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAKE 2D PLOTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    std::vector<TH2D*> hAllowedRegion_AlphaVsE0;
    std::vector<TH2D*> hAllowedRegion_LumVsE0;
    std::vector<TH2D*> hAllowedRegion_LumVsAlpha;
    
    for(int iFile = 0; iFile < files.size(); ++iFile){
        TString name_AlphaVsE0 = "alphavse0_" + options[iFile];
        TString name_LumVsE0 = "lumvse0_" + options[iFile];
        TString name_LumVsAlpha = "lumvsalpha_" + options[iFile];
        
        TH2D *tmp_AlphaVsE0 = new TH2D(name_AlphaVsE0, ";#LT E_{#nu} #GT (MeV);#alpha", nume0bins, mine0, maxe0, numalphabins, minalpha, maxalpha );
        TH2D *tmp_LumVsE0 = new TH2D(name_LumVsE0, ";#LT E_{#nu} #GT (MeV);#varepsilon (10^{53} erg)", nume0bins, mine0, maxe0, numlumbins, minlum/factor, maxlum/factor );
        TH2D *tmp_LumVsAlpha = new TH2D(name_LumVsAlpha, ";#alpha;#varepsilon (10^{53} erg)", numalphabins, minalpha, maxalpha, numlumbins, minlum/factor, maxlum/factor );
        
        for(int i = 0; i < eventid[iFile].size(); ++i){
            tmp_AlphaVsE0->Fill(e0[iFile][i], alpha[iFile][i]);
            tmp_LumVsE0->Fill(e0[iFile][i], lum[iFile][i]);
            tmp_LumVsAlpha->Fill(alpha[iFile][i], lum[iFile][i]);
        }
        
        //set some aesthetic stuff here
        tmp_AlphaVsE0->GetXaxis()->CenterTitle();
        tmp_AlphaVsE0->GetXaxis()->SetTitleSize(0.06);
        tmp_AlphaVsE0->GetXaxis()->SetLabelSize(0.06);
        tmp_AlphaVsE0->GetYaxis()->CenterTitle();
        tmp_AlphaVsE0->GetYaxis()->SetTitleSize(0.06);
        tmp_AlphaVsE0->GetYaxis()->SetTitleOffset(1.0);
        tmp_AlphaVsE0->GetYaxis()->SetLabelSize(0.06);
        
        tmp_LumVsE0->GetXaxis()->CenterTitle();
        tmp_LumVsE0->GetXaxis()->SetTitleSize(0.06);
        tmp_LumVsE0->GetXaxis()->SetLabelSize(0.06);
        tmp_LumVsE0->GetYaxis()->CenterTitle();
        tmp_LumVsE0->GetYaxis()->SetTitleSize(0.06);
        tmp_LumVsE0->GetYaxis()->SetLabelSize(0.06);
    
        tmp_LumVsAlpha->GetXaxis()->CenterTitle();
        tmp_LumVsAlpha->GetXaxis()->SetTitleSize(0.06);
        tmp_LumVsAlpha->GetXaxis()->SetLabelSize(0.06);
        tmp_LumVsAlpha->GetYaxis()->CenterTitle();
        tmp_LumVsAlpha->GetYaxis()->SetTitleSize(0.06);
        tmp_LumVsAlpha->GetYaxis()->SetLabelSize(0.06);
                
        //smooth them here b/c we only care about contour lines
        tmp_AlphaVsE0->Smooth();
        tmp_LumVsE0->Smooth();
        tmp_LumVsAlpha->Smooth();
        
        //keep them
        hAllowedRegion_AlphaVsE0.emplace_back(tmp_AlphaVsE0);
        hAllowedRegion_LumVsE0.emplace_back(tmp_LumVsE0);
        hAllowedRegion_LumVsAlpha.emplace_back(tmp_LumVsAlpha);
    }
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DRAW 2D PLOTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    TCanvas *c1 = new TCanvas("c1", "c1", 1000, 700);
    c1->Divide(2, 2);//, 0.02, 0.02); //nx, ny, xmargin, ymargin
    
    c1->cd(1);
    //TPad *subpad1 = (TPad*)c1->GetPad(1);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.14);
    for(size_t i = 0; i < hAllowedRegion_AlphaVsE0.size(); ++i){
        
        hAllowedRegion_AlphaVsE0[i] = setTH2Dcontour(hAllowedRegion_AlphaVsE0[i], 0.9);
                
        hAllowedRegion_AlphaVsE0[i]->SetLineWidth(5);
        if(i+1==5) hAllowedRegion_AlphaVsE0[i]->SetLineColor(i+2);
        else if(i+1==3) hAllowedRegion_AlphaVsE0[i]->SetLineColor(8);
        else hAllowedRegion_AlphaVsE0[i]->SetLineColor(i+1);
        
        hAllowedRegion_AlphaVsE0[i]->GetXaxis()->SetRangeUser(6, 13);
        hAllowedRegion_AlphaVsE0[i]->GetYaxis()->SetRangeUser(1, 5);
        
        if(i == 0) hAllowedRegion_AlphaVsE0[i]->Draw("cont3");
        else hAllowedRegion_AlphaVsE0[i]->Draw("cont3 same");
    }
    OGpoint_AlphaVsE0->Draw("P SAME");
    
    TLegend *leg = makeLegend(0.1, 0.1, 0.9, 0.9);
    leg->SetBorderSize(0);
    for(size_t i = 0; i < hAllowedRegion_AlphaVsE0.size(); ++i) leg->AddEntry(hAllowedRegion_AlphaVsE0[i], names[i], "pe");
    
    c1->cd(2);
    leg->SetTextSize(0.06);
    leg->Draw();
    
    c1->cd(3);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.01);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.16);
    
    for(size_t i = 0; i < hAllowedRegion_LumVsE0.size(); ++i){
        
        hAllowedRegion_LumVsE0[i] = setTH2Dcontour(hAllowedRegion_LumVsE0[i], 0.9);
        hAllowedRegion_LumVsE0[i]->SetLineWidth(5);
        if(i+1==5) hAllowedRegion_LumVsE0[i]->SetLineColor(i+2);
        else if(i+1==3) hAllowedRegion_LumVsE0[i]->SetLineColor(8);
        else hAllowedRegion_LumVsE0[i]->SetLineColor(i+1);
        
        hAllowedRegion_LumVsE0[i]->GetXaxis()->SetRangeUser(6, 13);
        
        if(i == 0) hAllowedRegion_LumVsE0[i]->Draw("cont3");
        else hAllowedRegion_LumVsE0[i]->Draw("cont3 same");
    }
    OGpoint_LumVsE0->Draw("P SAME");
    
    c1->cd(4);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.01);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.16);
    
    for(size_t i = 0; i < hAllowedRegion_LumVsAlpha.size(); ++i){
        
        hAllowedRegion_LumVsAlpha[i] = setTH2Dcontour(hAllowedRegion_LumVsAlpha[i], 0.9);
        hAllowedRegion_LumVsAlpha[i]->SetLineWidth(5);
        if(i+1==5) hAllowedRegion_LumVsAlpha[i]->SetLineColor(i+2);
        else if(i+1==3) hAllowedRegion_LumVsAlpha[i]->SetLineColor(8);
        else hAllowedRegion_LumVsAlpha[i]->SetLineColor(i+1);
        
        if(i == 0) hAllowedRegion_LumVsAlpha[i]->Draw("cont3");
        else hAllowedRegion_LumVsAlpha[i]->Draw("cont3 same");
    }
    OGpoint_LumVsAlpha->Draw("P SAME");
    
    TCanvas *c2 = new TCanvas("c2", "Compare to Nikrant", 800,700);
    gPad->SetTopMargin(0.015);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.14);
    gPad->SetLeftMargin(0.19);
    
    hAllowedRegion_LumVsE0[0]->GetYaxis()->SetRangeUser(0.3, 1.0);
    hAllowedRegion_LumVsE0[0]->Draw("cont3");
    for(int i = 1; i < hAllowedRegion_LumVsE0.size(); ++i){
        
        hAllowedRegion_LumVsE0[i]->Draw("cont3 same");
        
    }
    TLegend *leg2 = makeLegend(0.45,0.65,0.8,0.95);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.035);
    leg2->SetHeader("Fake supernovae study","C");
    for(size_t i = 0; i < hAllowedRegion_LumVsE0.size(); ++i){
        leg2->AddEntry(hAllowedRegion_LumVsE0[i], names[i], "l");
    }
    leg2->AddEntry(OGpoint_LumVsE0, "Truth: #LT E_{#nu} #GT = 9.5 MeV, #varepsilon = 5e52 erg", "p");
    leg2->Draw("same");
    
    
    OGpoint_LumVsE0->Draw("P SAME");
    
    
}//end of code
