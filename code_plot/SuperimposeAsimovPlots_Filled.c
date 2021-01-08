//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// SuperimposeAsimovPlots_Filled
// Erin Conley (erin.conley@duke.edu)
// Description: Creates a TCanvas with superimposed contour plots produced by
//              the Asimov method. Specify which plots to superimpose in the
//              "grid_config" and "ts_config" vectors. Draws "filled" (or
//              shaded) contours instead of lines
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

void SuperimposeAsimovPlots_Filled(){
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // HERE ARE THE VARIABLES WE CHANGE
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TString distance("10.00");
    
    //name of the legend
    //TString legendTitle("");
    TString legendTitle("10 kpc supernova, 90% C.L.");
    
    //define extra directory where the files live, e.g., "/CrossSectionStudy/"
    //ps: if there isn't an extra directory, set the parameter to "/"
    TString dir("/");
    
    std::vector<TString> grid_config;
    std::vector<TString> ts_config;
    std::vector<TString> names;
   
    //example configurations: mass ordering assumptions
    grid_config.emplace_back("2019October17_smearMARLEY_nueOnly_PNXscn_StepEffic");
    ts_config.emplace_back("spectraMARLEY_nueOnly_PNXscn_StepEffic");
    names.emplace_back("No mass ordering assumptions");
    
    grid_config.emplace_back("2019October17_smearMARLEY_NO_nueOnly_PNXscn_StepEffic");
    ts_config.emplace_back("spectraMARLEY_NO_nueOnly_PNXscn_StepEffic");
    names.emplace_back("Normal mass ordering assumptions");
    
    grid_config.emplace_back("2019October17_smearMARLEY_IO_nueOnly_PNXscn_StepEffic");
    ts_config.emplace_back("spectraMARLEY_IO_nueOnly_PNXscn_StepEffic");
    names.emplace_back("Inverted mass ordering assumptions");
    
    std::vector<TString> histNames;
	histNames.emplace_back("alphavse0_90cut");
	histNames.emplace_back("lumvse0_90cut");
	histNames.emplace_back("lumvsalpha_90cut");
    
    std::vector<TString> histNamesForInterpolation;
    histNamesForInterpolation.emplace_back("alphavse0");
    histNamesForInterpolation.emplace_back("lumvse0");
    histNamesForInterpolation.emplace_back("lumvsalpha");
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE INPUTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //define file names
    
    std::vector<TString> files;
	for(size_t i = 0; i < names.size(); ++i){
		
        TString tmp = "out" + dir + "chi2plots_" + grid_config[i] + "_" + ts_config[i] + "_" + distance + "kpc.root";
        
		files.emplace_back(tmp);
	}
	
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
    // DEFINE TRUTH POINTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //define truth point of test spectrum
    Double_t OGalpha[1];
    Double_t OGe0[1];
    Double_t OGlum[1];
    
    Double_t factor = 1e53;
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
		
		//TString label = "Best-Fit: #alpha = " + TString::Format("%4.1f",tmpalpha[0]) + ", #LT E_{#nu} #GT = " + TString::Format("%4.1f",tmpe0[0]) + " MeV, #varepsilon = " + TString::Format("%4.1f",tmplum[0]) + "e53 ergs";
        TString label = "Best-Fit: #alpha = " + TString::Format("%4.1f",tmpalpha[0]) + ", #LT E_{#nu} #GT = " + TString::Format("%4.1f",tmpe0[0]) + ", #varepsilon = " + TString::Format("%4.1f",tmplum[0]) + "e53";
        
		
	    TGraph *tmppoint_AlphaVsE0 = new TGraph(1, tmpe0, tmpalpha);
		tmppoint_AlphaVsE0->SetMarkerStyle(29);
		tmppoint_AlphaVsE0->SetMarkerColor(kBlack);
		tmppoint_AlphaVsE0->SetMarkerSize(2.0);
    
		TGraph *tmppoint_LumVsE0 = new TGraph(1, tmpe0, tmplum);
		tmppoint_LumVsE0->SetMarkerStyle(29);
		tmppoint_LumVsE0->SetMarkerColor(kBlack);
		tmppoint_LumVsE0->SetMarkerSize(2.0);
		
		TGraph *tmppoint_LumVsAlpha = new TGraph(1, tmpalpha, tmplum);
		tmppoint_LumVsAlpha->SetMarkerStyle(29);
		tmppoint_LumVsAlpha->SetMarkerColor(kBlack);
		tmppoint_LumVsAlpha->SetMarkerSize(2.0);
		
		if(i+1==5){ 
			tmppoint_AlphaVsE0->SetMarkerColor(kYellow+2);
			tmppoint_LumVsE0->SetMarkerColor(kYellow+2);
			tmppoint_LumVsAlpha->SetMarkerColor(kYellow+2);
		}
		else if(i+1==3){ 
			tmppoint_AlphaVsE0->SetMarkerColor(kGreen+2);
			tmppoint_LumVsE0->SetMarkerColor(kGreen+2);
			tmppoint_LumVsAlpha->SetMarkerColor(kGreen+2);
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
			TH2D *h1 = (TH2D*)gDirectory->Get(histNames[iHist]);
            
            
            TH2D *hToInterpolate = (TH2D*)gDirectory->Get(histNamesForInterpolation[iHist]);
            TH2D *h = InterpolateChi2Hist(hToInterpolate, 1000, 4.61);
            
            //std::cout << std::abs(findAreaOfContour(h) - findAreaOfContour(h1))/findAreaOfContour(h) << std::endl;
            
            if( std::abs(findAreaOfContour(h) - findAreaOfContour(h1))/findAreaOfContour(h) > 0.7 || (iHist == 0 && iFile == 0 ) ){
                
                h->Reset();
                h = (TH2D*)gDirectory->Get(histNames[iHist]);
                
                
            }
            
            
			h->SetContour(1, contours);
			//h->SetLineWidth(6);
			
			//output contour area 
			//std::cout << "Area for " << files[iFile] << ", " << histNames[iHist] << ": " << findAreaOfContour(h) << std::endl; 
			
			//also try to make the labels/axis titles better
			h->GetXaxis()->CenterTitle();
			h->GetXaxis()->SetTitleOffset(1.0);
			h->GetXaxis()->SetTitleSize(0.05);
			h->GetXaxis()->SetLabelSize(0.07);
			h->GetYaxis()->CenterTitle();
            h->GetYaxis()->SetTitleOffset(1.0);
			h->GetYaxis()->SetTitleSize(0.05);
			h->GetYaxis()->SetLabelSize(0.07);
						
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
    for(size_t i = 0; i < hists[0].size(); ++i) leg->AddEntry(hists[0][i], names[i], "f");
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
	
	//also draw the truth points
    
	c1->cd(1); 
	//OGpoint_AlphaVsE0->Draw("P SAME");
    for(size_t i = 0; i < BestFitPoints_AlphaVsE0.size(); ++i) BestFitPoints_AlphaVsE0[i]->Draw("P SAME");
    
	c1->cd(3); 
	//OGpoint_LumVsE0->Draw("P SAME");
	for(size_t i = 0; i < BestFitPoints_LumVsE0.size(); ++i) BestFitPoints_LumVsE0[i]->Draw("P SAME");
	c1->cd(4); 
	//OGpoint_LumVsAlpha->Draw("P SAME");
	for(size_t i = 0; i < BestFitPoints_LumVsAlpha.size(); ++i)            BestFitPoints_LumVsAlpha[i]->Draw("P SAME");
	
    
    TCanvas *c2 = new TCanvas("c2", "Compare to Nikrant", 800,700);
    gPad->SetTopMargin(0.015);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.14);
    gPad->SetLeftMargin(0.14);
    
    for(int i = 0; i < hists.size(); ++i){
        
        int toDraw = i+1;
        if(i >= 1) ++toDraw;
        
        if(toDraw != 3) continue; //only draw if lum vs e0
        
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
                axisHistograms[i]->GetXaxis()->SetTitleOffset(1.1);
                axisHistograms[i]->GetYaxis()->SetTitleOffset(1.1);
                axisHistograms[i]->Draw();
            }
            
            
            //hists[i][j]->Draw("cont3 same");
            hists[i][j]->Draw("cont0 same");
            
            
        }
        
    }
    TLegend *leg2 = makeLegend(0.45,0.65,0.8,0.95);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.04);
    leg2->SetHeader(legendTitle,"C");
    for(size_t i = 0; i < hists[2].size(); ++i){
        leg2->AddEntry(hists[2][i], names[i], "f");
    }
    leg2->AddEntry(BestFitPoints_LumVsE0[0], "Truth: #LT E_{#nu} #GT = 9.5 MeV, #varepsilon = 5e52 erg", "p");
    leg2->Draw("same");
    
    for(size_t i = 0; i < BestFitPoints_LumVsE0.size(); ++i){
        //BestFitPoints_LumVsE0[i]->SetMarkerSize(2.5);
        //BestFitPoints_LumVsE0[i]->SetMarkerColor(kBlack);
        BestFitPoints_LumVsE0[i]->Draw("P SAME");
    }
    
}//end of code
