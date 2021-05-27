//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// SuperimposeDifferentDistances
// Erin Conley (erin.conley@duke.edu)
// Description: Creates a TCanvas with superimposed contour plots produced by
//              the Asimov method and studying different supernova distances.
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

#include "../headers/ForwardFit.h"
#include "../headers/Asimov.h"

void SuperimposeDifferentDistances(){
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // HERE ARE THE VARIABLES WE CHANGE
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //define extra directory where the files live, e.g., "/CrossSectionStudy/"
    //if there isn't an extra directory, set the parameter to "/"
    TString dir("/ResolutionVsDistanceStudies/");
    
    //define the smearing we want to use
    TString smear("MARLEY_DCR_nueOnly");
	
	//define the cross section we want to use
	TString xscn("PNXscn");
	
	//define the efficiency we want to use
	TString effic("StepEffic");
	
	//define the grid we want to use
    TString grid("2019October17");
    
	std::vector<TString> distances;
	distances.emplace_back("10.00kpc");
    distances.emplace_back("7.00kpc");
    distances.emplace_back("4.00kpc");
	
	std::vector<TString> graphNames;
	graphNames.emplace_back("Chi2VsAlpha;1");
	graphNames.emplace_back("Chi2VsE0;1");
	graphNames.emplace_back("Chi2VsLum;1");
	
	std::vector<TString> histNames;
	histNames.emplace_back("alphavse0_90cut");
	histNames.emplace_back("lumvse0_90cut");
	histNames.emplace_back("lumvsalpha_90cut");
    
    std::vector<TString> histNamesForInterpolation;
    histNamesForInterpolation.emplace_back("alphavse0");
    histNamesForInterpolation.emplace_back("lumvse0");
    histNamesForInterpolation.emplace_back("lumvsalpha");
    
    //name of the legend
    TString legendTitle("90% C.L.");
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE INPUTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//file to open
	std::vector<TString> files;
	for(size_t i = 0; i < distances.size(); ++i){
		
		TString tmp = "out" + dir + "chi2plots_" + grid + "_smear" + smear + "_" + xscn + "_" + effic + "_spectra" + smear + "_" + xscn + "_" + effic + "_" + distances[i] + ".root";
		files.emplace_back(tmp);
	}
	
	TString gridinfo = "input/grid_info/grid_info_"+grid+".dat";
    
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
    // GRAB CHI2 PLOTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    std::vector<std::vector<TGraph*>> graphs;
    std::vector<std::vector<TString>> graphTitles;
	for(int iGraph = 0; iGraph < graphNames.size(); ++iGraph){  
		std::vector<TGraph*> tmp;
		std::vector<TString> tmpname;
		for(int iFile = 0; iFile < files.size(); ++iFile){
			TFile *f = new TFile(files[iFile]);
			TGraph *g = (TGraph*)gDirectory->Get(graphNames[iGraph]);
			TString title = g->GetTitle();
			tmp.emplace_back(g);
			tmpname.emplace_back(title);
			f->Close();
		}
		graphs.emplace_back(tmp);
		graphTitles.emplace_back(tmpname);
	}
	
	std::vector<std::vector<TH2D*>> hists;
	for(int iHist = 0; iHist < graphNames.size(); ++iHist){  
		std::vector<TH2D*> tmp;
		for(int iFile = 0; iFile < files.size(); ++iFile){
			TFile *f = new TFile(files[iFile]);
			//TH2D *h = (TH2D*)gDirectory->Get(histNames[iHist]);
            
            TH2D *hToInterpolate = (TH2D*)gDirectory->Get(histNamesForInterpolation[iHist]);
            TH2D *h = InterpolateChi2Hist(hToInterpolate, 1000, 4.61);
            
			h->SetContour(1, contours);
			h->SetLineWidth(3);
						
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
    c1->cd(2);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.01);
    gPad->SetRightMargin(0.015);
    TLegend *leg = makeLegend(0.1,0.1,0.99,0.99);
    leg->SetTextSize(0.07);
    leg->SetHeader(legendTitle,"C");
    for(size_t i = 0; i < hists[0].size(); ++i) leg->AddEntry(hists[0][i], distances[i], "f");
    leg->Draw();
    
    std::vector<TH2D*> axisHistograms;
    for(int i = 0; i < hists.size(); ++i){
        
        int toDraw = i+1;
        if(i >= 1) ++toDraw;
        //hists[i][0]
        
        TString histAxisName = "axis_" + TString::Format("%d",toDraw);
        
        double xmin;//= hists[i][0]->GetXaxis()->GetXmin();//*0.5;
        double xmax;//= hists[i][0]->GetXaxis()->GetXmax();//*1.2;
        double ymin;//= hists[i][0]->GetYaxis()->GetXmin();//0.0;
        //if(toDraw == 1 || toDraw == 3) xmax = 15.0;
        double ymax;//= hists[i][0]->GetYaxis()->GetXmax();//*1.2;
        
        getMinOrMaxFromHists(hists[i], xmin, xmax, ymin, ymax);
        
        TString xTitle("");
        TString yTitle("");
        
        if(toDraw == 1){ //alpha vs e0
            xTitle = "#LT E_{#nu} #GT (MeV)";
            yTitle = "#alpha";
        }
        else if(toDraw == 3){ //lum vs e0
            xTitle = "#LT E_{#nu} #GT (MeV)";
            yTitle = "#varepsilon (10^{53} erg)";
        }
        else{ //lum vs alpha
            xTitle = "#alpha";
            yTitle = "#varepsilon (10^{53} erg)";
        }
        
        TH2D *histAxis = new TH2D(histAxisName, " ;;", 2, xmin, xmax, 2, ymin, ymax);
        histAxis->GetXaxis()->SetTitle( xTitle );
        histAxis->GetYaxis()->SetTitle( yTitle );
        
        histAxis->GetXaxis()->CenterTitle();
        histAxis->GetXaxis()->SetTitleOffset(1.1);
        histAxis->GetXaxis()->SetTitleSize(0.06);
        histAxis->GetXaxis()->SetLabelSize(0.07);
        histAxis->GetYaxis()->CenterTitle();
        histAxis->GetYaxis()->SetTitleOffset(1.0);
        if(toDraw == 3 || toDraw == 4){
            histAxis->GetYaxis()->SetTitle("#varepsilon (10^{53} erg)");
            histAxis->GetYaxis()->SetTitleOffset(1.1);
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
        
        for(int j = 0; j < hists[i].size(); ++j){
            
            if(j+1==5){
                hists[i][j]->SetLineColor(kYellow+2);
                hists[i][j]->SetFillColorAlpha(kYellow+2, 0.5);
            }
            else if(j+1==3){
                hists[i][j]->SetLineColor(kGreen+2);
                hists[i][j]->SetFillColorAlpha(kGreen+2, 0.5);
            }
            else{
                hists[i][j]->SetLineColor(j+1);
                hists[i][j]->SetFillColorAlpha(j+1, 0.5);
            }
            
            if(j == 0){
                axisHistograms[i]->GetXaxis()->SetLabelSize(0.06);
                axisHistograms[i]->GetYaxis()->SetLabelSize(0.06);
                axisHistograms[i]->GetXaxis()->SetTitleSize(0.06);
                axisHistograms[i]->GetYaxis()->SetTitleSize(0.06);
                axisHistograms[i]->GetXaxis()->SetTitleOffset(1.3);
                axisHistograms[i]->GetYaxis()->SetTitleOffset(1.1);
                axisHistograms[i]->Draw();
            }
            
            
            hists[i][j]->Draw("cont0 same");
        }
        
    }
    c1->cd(1);
    OGpoint_AlphaVsE0->Draw("P SAME");
    c1->cd(3);
    OGpoint_LumVsE0->Draw("P SAME");
    c1->cd(4);
    OGpoint_LumVsAlpha->Draw("P SAME");
    
    TCanvas *c2 = new TCanvas("c2", "Compare to Nikrant", 800,700);
    gPad->SetTopMargin(0.015);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.14);
    gPad->SetLeftMargin(0.14);
    
    TLegend *leg2 = makeLegend(0.1,0.7,0.3,0.85);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.04);
    for(size_t i = 0; i < hists[0].size(); ++i) leg2->AddEntry(hists[0][i], distances[i], "f");
    leg2->SetHeader(legendTitle,"C");
    leg2->AddEntry(OGpoint_LumVsE0, "Truth: #LT E_{#nu} #GT = 9.5 MeV, #varepsilon = 5e52 erg", "p");
    for(int i = 0; i < hists.size(); ++i){
        
        int toDraw = i+1;
        if(i >= 1) ++toDraw;
        
        if(toDraw != 3) continue; //only draw if lum vs e0
        
        for(int j = 0; j < hists[i].size(); ++j){
            
            hists[i][j]->SetLineWidth(5);
            
            if(j+1==5){
                hists[i][j]->SetLineColor(kYellow+2);
                hists[i][j]->SetFillColorAlpha(kYellow+2, 0.5);
            }
            else if(j+1==3){
                hists[i][j]->SetLineColor(kGreen+2);
                hists[i][j]->SetFillColorAlpha(kGreen+2, 0.5);
            }
            else{
                hists[i][j]->SetLineColor(j+1);
                hists[i][j]->SetFillColorAlpha(j+1, 0.5);
            }
            
            
            if(j == 0){
                axisHistograms[i]->GetXaxis()->SetLimits(5.0, 20.0);
                axisHistograms[i]->GetYaxis()->SetLimits(0.2, 1.0);
                axisHistograms[i]->GetXaxis()->SetLabelSize(0.06);
                axisHistograms[i]->GetYaxis()->SetLabelSize(0.06);
                axisHistograms[i]->GetXaxis()->SetTitleSize(0.06);
                axisHistograms[i]->GetYaxis()->SetTitleSize(0.06);
                axisHistograms[i]->GetXaxis()->SetTitleOffset(1.1);
                axisHistograms[i]->GetYaxis()->SetTitleOffset(1.1);
                axisHistograms[i]->Draw();
            }
            
            
            //hists[i][j]->Draw("cont3 same");
            hists[i][j]->Draw("cont0 same");
            
            
        }
        
    }
    OGpoint_LumVsE0->Draw("P same");
    leg2->Draw("same");

}//end of code
