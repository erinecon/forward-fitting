//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// InterpolateSensitivityRegions
// Erin Conley (erin.conley@duke.edu)
// Description: Test interpolation techniques for TH2D contours made by the
//              Asimov method
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

void InterpolateSensitivityRegions(){
	
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // HERE ARE THE VARIABLES WE CHANGE
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    TString distance("10.00");
    
    //define extra directory where the files live, e.g., "/CrossSectionStudy/"
    //ps: if there isn't an extra directory, set the parameter to "/"
    TString dir("/");
    
    TString grid_config("2019October17_smearMARLEY_nueOnly_PNXscn_StepEffic");
    TString ts_config("spectraMARLEY_nueOnly_PNXscn_StepEffic");
    TString name("MARLEY + (p, n) xscn");
    
    std::vector<TString> histNames;
	histNames.emplace_back("alphavse0_90cut");
	histNames.emplace_back("lumvse0_90cut");
	histNames.emplace_back("lumvsalpha_90cut");
    
    std::vector<TString> histNamesForInterpolation;
    histNamesForInterpolation.emplace_back("alphavse0");
    histNamesForInterpolation.emplace_back("lumvse0");
    histNamesForInterpolation.emplace_back("lumvsalpha");
	
	//define the grid we want to use
    TString grid("2019October17");
    
    //name of the legend
    TString legendTitle("10 kpc supernova");
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE INPUTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    TString filename = "out" + dir + "chi2plots_" + grid_config + "_" + ts_config + "_" + distance + "kpc.root";
	
	TString gridinfo = "input/grid_info/grid_info_" + grid + ".dat";
    
    //set specific aesthetic things here
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1,0);
    gStyle->SetNumberContours(99);
    
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
    // GRAB EVENT INFORMATION
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TString a("alpha");
    TString e("e0");
    TString l("lum");
    
    std::vector<double> alpha = getParameter(filename, a);
    std::vector<double> e0 = getParameter(filename, e);
    std::vector<double> lum = getParameter(filename, l);
    
    Double_t tmpalpha[1];
    Double_t tmpe0[1];
    Double_t tmplum[1];
    
    tmpalpha[0] = alpha[0];
    tmpe0[0] = e0[0];
    tmplum[0] = lum[0];
    
    TGraph* BestFitPoints_AlphaVsE0 = new TGraph(1, tmpe0, tmpalpha);
    BestFitPoints_AlphaVsE0->SetMarkerStyle(29);
    BestFitPoints_AlphaVsE0->SetMarkerColor(kBlack);
    BestFitPoints_AlphaVsE0->SetMarkerSize(2.0);
    
    TGraph* BestFitPoints_LumVsE0 = new TGraph(1, tmpe0, tmplum);
    BestFitPoints_LumVsE0->SetMarkerStyle(29);
    BestFitPoints_LumVsE0->SetMarkerColor(kBlack);
    BestFitPoints_LumVsE0->SetMarkerSize(2.0);
    
    TGraph* BestFitPoints_LumVsAlpha = new TGraph(1, tmpalpha, tmplum);
    BestFitPoints_LumVsAlpha->SetMarkerStyle(29);
    BestFitPoints_LumVsAlpha->SetMarkerColor(kBlack);
    BestFitPoints_LumVsAlpha->SetMarkerSize(2.0);
    
    TString BestFitPoints_label = "Best-Fit: #alpha = " + TString::Format("%4.1f",tmpalpha[0]) + ", #LT E_{#nu} #GT = " + TString::Format("%4.1f",tmpe0[0]) + ", #varepsilon = " + TString::Format("%4.1f",tmplum[0]) + "e53";
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GRAB CHI2 PLOTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
	std::vector<TH2D*> hists;
    std::vector<TH2D*> hists_colormap;
    std::vector<TH2D*> hists_interpolated_colormap;
    std::vector<TH2D*> hists_interpolated;
    
    
	for(int iHist = 0; iHist < histNames.size(); ++iHist){
		
        TFile *f = new TFile(filename);
		TH2D *h = (TH2D*)gDirectory->Get(histNames[iHist]);
        
        //clone histogram to keep
        TH2D *h_colormap = (TH2D*)h->Clone();
        
		h->SetContour(1, contours);
        h->SetLineWidth(2);
			
        //also try to make the labels/axis titles better
        h->GetXaxis()->CenterTitle();
        h->GetXaxis()->SetTitleOffset(1.0);
        h->GetXaxis()->SetTitleSize(0.05);
		h->GetXaxis()->SetLabelSize(0.07);
        h->GetYaxis()->CenterTitle();
        h->GetYaxis()->SetTitleOffset(1.0);
        h->GetYaxis()->SetTitleSize(0.05);
        h->GetYaxis()->SetLabelSize(0.07);
        
        TH2D *hToInterpolate = (TH2D*)gDirectory->Get(histNamesForInterpolation[iHist]);
        
        TH2D *h_interpolated = InterpolateChi2Hist(hToInterpolate, 1000, 4.61);
        
        //clone histogram to keep
        TH2D *h_inter_colormap = (TH2D*)h_interpolated->Clone();
        
        h_interpolated->SetContour(1, contours);
        h_interpolated->SetLineWidth(2);
		
		hists.emplace_back(h);
        hists_colormap.emplace_back(h_colormap);
        hists_interpolated_colormap.emplace_back(h_inter_colormap);
        hists_interpolated.emplace_back(h_interpolated);
	}
    
    TCanvas *ColorMap = new TCanvas("ColorMap", "ColorMap", 800,700);
    hists_colormap[1]->SetTitle("Raw colormap with superimposed contour");
    hists_colormap[1]->Draw("colz");
    hists[1]->Draw("cont3 same");
    
    ColorMap = new TCanvas("ColorMapInterpolated", "ColorMapInterpolated", 800,700);
    hists_interpolated_colormap[1]->SetTitle("Interpolated colormap with superimposed contour;#LT E_{#nu} #GT (MeV);#varepsilon (10^{53} erg)");
    hists_interpolated_colormap[1]->Draw("colz");
    hists_interpolated[1]->Draw("cont3 same");
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DRAW PLOTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
	TCanvas *c1 = new TCanvas("c1", "c1", 1000, 700);
    c1->Divide(2, 2);
    c1->cd(2);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.01);
    gPad->SetRightMargin(0.015);
	TLegend *leg = makeLegend(0.1,0.1,0.99,0.99);
    leg->SetTextSize(0.07);
    leg->SetHeader(legendTitle,"C");
    leg->AddEntry(hists[0], name);
    leg->AddEntry(hists_interpolated[0], "TH2D Interpolate");
    
	leg->Draw();
    
    
    std::vector<TH2D*> axisHistograms;
    for(int i = 0; i < hists.size(); ++i){
        
        int toDraw = i+1;
        if(i >= 1) ++toDraw;
        //hists[i][0]
        
        TString histAxisName = "axis_" + TString::Format("%d",toDraw);
        
        double xmin= hists[i]->GetXaxis()->GetXmin();
        double xmax= hists[i]->GetXaxis()->GetXmax();
        double ymin= hists[i]->GetYaxis()->GetXmin();
        //if(toDraw == 1 || toDraw == 3) xmax = 15.0;
        double ymax= hists[i]->GetYaxis()->GetXmax();
        
        
        if(toDraw == 1){ //alpha vs e0
            xmax-=0.06;
            ymin+=0.14;
            ymax-=0.14;
        }
        else if(toDraw == 3){ //lum vs e0
            xmax-=0.06;
            
        }
        else{ //lum vs alpha
            xmin+=0.14;
            xmax-=0.14;
            
        }
        
        TH2D *histAxis = new TH2D(histAxisName, " ;;", 2, xmin, xmax, 2, ymin, ymax);
        histAxis->GetXaxis()->SetTitle( hists[i]->GetXaxis()->GetTitle() );
        histAxis->GetYaxis()->SetTitle( hists[i]->GetYaxis()->GetTitle() );
        
        histAxis->GetXaxis()->CenterTitle();
        histAxis->GetXaxis()->SetTitleOffset(1.1);
        histAxis->GetXaxis()->SetTitleSize(0.06);
        histAxis->GetXaxis()->SetLabelSize(0.07);
        histAxis->GetYaxis()->CenterTitle();
        histAxis->GetYaxis()->SetTitleOffset(1.0);
        if(toDraw == 3 || toDraw == 4){
            histAxis->GetYaxis()->SetTitleOffset(1.1);
            histAxis->GetYaxis()->SetTitle("#varepsilon (10^{53} erg)");
        }
        histAxis->GetYaxis()->SetTitleSize(0.06);
        histAxis->GetYaxis()->SetLabelSize(0.07);
        
        axisHistograms.emplace_back(histAxis);
        
    }
    
	
    for(int i = 0; i < hists.size(); ++i){
		
		int toDraw = i+1;
		if(i >= 1) ++toDraw;
		c1->cd(toDraw);
        gPad->SetBottomMargin(0.18);
        gPad->SetTopMargin(0.01);
        gPad->SetRightMargin(0.04);
        gPad->SetLeftMargin(0.13);
		
        hists[i]->SetLineColor(kBlack);
        hists_interpolated[i]->SetLineColor(kBlue);
        
        axisHistograms[i]->GetXaxis()->SetLabelSize(0.06);
        axisHistograms[i]->GetYaxis()->SetLabelSize(0.06);
        axisHistograms[i]->GetXaxis()->SetTitleSize(0.06);
        axisHistograms[i]->GetYaxis()->SetTitleSize(0.06);
        axisHistograms[i]->GetXaxis()->SetTitleOffset(1.3);
        axisHistograms[i]->GetYaxis()->SetTitleOffset(1.1);
        axisHistograms[i]->Draw();
        
        hists[i]->Draw("cont3 same");
        hists_interpolated[i]->Draw("cont3 same");
		
	}
	
	//also draw the truth points
	c1->cd(1); 
	OGpoint_AlphaVsE0->Draw("P SAME");
    BestFitPoints_AlphaVsE0->Draw("P SAME");
    
	c1->cd(3); 
	OGpoint_LumVsE0->Draw("P SAME");
	BestFitPoints_LumVsE0->Draw("P SAME");
	c1->cd(4); 
	OGpoint_LumVsAlpha->Draw("P SAME");
	BestFitPoints_LumVsAlpha->Draw("P SAME");
	
    
    TCanvas *c2 = new TCanvas("c2", "Compare to Nikrant", 800,700);
    gPad->SetTopMargin(0.015);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.14);
    gPad->SetLeftMargin(0.14);
    
    TLegend *leg2 = makeLegend(0.45,0.65,0.85,0.95);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.035);
    leg2->SetHeader(legendTitle,"C");
    leg2->AddEntry(hists[0], name);
    OGpoint_LumVsE0->SetMarkerColor(kBlack);
    leg2->AddEntry(hists_interpolated[0], "TH2D Interpolate");
    leg2->AddEntry(OGpoint_LumVsE0, "Truth: #LT E_{#nu} #GT = 9.5 MeV, #varepsilon = 5e52 erg", "p");
    for(int i = 0; i < hists.size(); ++i){
        
        int toDraw = i+1;
        if(i >= 1) ++toDraw;
        
        if(toDraw != 3) continue; //only draw if lum vs e0
        
        hists[i]->SetLineColor(kBlack);
        hists_interpolated[i]->SetLineColor(kBlue);
        
        
        axisHistograms[i]->GetXaxis()->SetLabelSize(0.06);
        axisHistograms[i]->GetYaxis()->SetLabelSize(0.06);
        axisHistograms[i]->GetXaxis()->SetTitleSize(0.06);
        axisHistograms[i]->GetYaxis()->SetTitleSize(0.06);
        axisHistograms[i]->GetXaxis()->SetTitleOffset(1.1);
        axisHistograms[i]->GetYaxis()->SetTitleOffset(1.1);
        axisHistograms[i]->Draw();
        
        
        hists[i]->Draw("cont3 same");
        hists_interpolated[i]->Draw("cont3 same");
        
    }
    leg2->Draw("same");
    OGpoint_LumVsE0->Draw("P SAME");
    
    
}//end of code
