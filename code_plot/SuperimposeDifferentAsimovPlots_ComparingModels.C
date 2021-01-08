//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// SuperimposeDifferentAsimovPlots_ComparingModels
// Erin Conley (erin.conley@duke.edu)
// Description: Superimpose Asimov study results with flux parameter values
//              corresponding to different supernova model databases (Nakazato
//              and Huedepohl)
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

void SuperimposeDifferentAsimovPlots_ComparingModels(){
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // HERE ARE THE VARIABLES WE CHANGE
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TString distance("10.00");
    
    //define sub-directory to store flux parameter values for the SN database(s)
    //store sub-directory in input directory
    //if there isn't an extra directory, set the parameter to "/"
    TString dir("/DifferentModels/");
    
    std::vector<TString> grid_config;
    std::vector<TString> ts_config;
    std::vector<TString> names_config;
    
    grid_config.emplace_back( "2019October17_smearMARLEY_DCR_nueOnly_PNXscn_StepEffic");
    ts_config.emplace_back("spectraMARLEY_DCR_nueOnly_PNXscn_StepEffic");
    names_config.emplace_back("#splitline{MARLEY smearing + (p, n)}{xscn + 5 MeV detection thresh.}");
    
    std::vector<TString> otherModelNames; std::vector<TString> otherModelLegNames;
    otherModelNames.emplace_back("nakazato_fitParameters_nueOnly_table.dat"); otherModelLegNames.emplace_back("Nakazato");
    //otherModelNames.emplace_back("Huedepohl_fitParameters_nueOnly.dat"); otherModelLegNames.emplace_back("Huedepohl");
    otherModelNames.emplace_back("Huedepohl_fitParameters_nueOnly_BlackHole.dat"); otherModelLegNames.emplace_back("Huedepohl, Black Hole");
    otherModelNames.emplace_back("Huedepohl_fitParameters_nueOnly_Cooling.dat"); otherModelLegNames.emplace_back("Huedepohl, Cooling");
    
    std::vector<TString> histNames;
	histNames.emplace_back("alphavse0_90cut");
	histNames.emplace_back("lumvse0_90cut");
	histNames.emplace_back("lumvsalpha_90cut");
    
    std::vector<TString> histNamesForInterpolation;
    histNamesForInterpolation.emplace_back("alphavse0");
    histNamesForInterpolation.emplace_back("lumvse0");
    histNamesForInterpolation.emplace_back("lumvsalpha");
    
    //name of the legend
    //TString legendTitle("");
    TString legendTitle("10kpc supernova, 90% C.L.");
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE INPUTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    std::vector<TString> files;
    for(size_t i = 0; i < grid_config.size(); ++i){
        
        TString tmp = "out" + "chi2plots_" + grid_config[i] + "_" + ts_config[i] + "_" + distance + "kpc.root";
        
        files.emplace_back(tmp);
    }
    
    std::vector<TString> file_OtherModels;
    for(int i = 0; i < otherModelNames.size(); ++i){
        TString tmpfile = "input" + dir + otherModelNames[i];
        file_OtherModels.emplace_back(tmpfile);
    }
    
    
    //set specific aesthetic things here
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1,0);
    
    //need to first define the truth points
    Double_t alpha_true = 2.5;
    Double_t e0_true = 9.5; //mev
    Double_t lum_true = 5e52; //ergs
    //define factor
    Double_t factor = 1e53;
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
    // GRAB OTHER SN DATABASE INFORMATION
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    std::vector<std::vector<TString>> code_OtherModels;
    std::vector<std::vector<Double_t>> alpha_OtherModels;
    std::vector<std::vector<Double_t>> alphaerr_OtherModels;
    std::vector<std::vector<Double_t>> meane_OtherModels;
    std::vector<std::vector<Double_t>> meaneerr_OtherModels;
    std::vector<std::vector<Double_t>> toten_OtherModels;
    std::vector<std::vector<Double_t>> totenerr_OtherModels;
    
    
    for(int i = 0; i < file_OtherModels.size(); ++i){
    
        //open grid_info.dat
        ifstream nin;
        nin.open(file_OtherModels[i]);
        
        std::vector<TString> tmp_code;
        std::vector<Double_t> tmp_alpha;
        std::vector<Double_t> tmp_alphaerr;
        std::vector<Double_t> tmp_meane;
        std::vector<Double_t> tmp_meaneerr;
        std::vector<Double_t> tmp_toten;
        std::vector<Double_t> tmp_totenerr;
        
        while(1){
            TString code;
            Double_t alpha, alphaerr, meane, meaneerr, toten, totenerr;
            nin >> code >> alpha >> alphaerr >> meane >> meaneerr >> toten >> totenerr;
            if(!nin.good()) break;
            
            tmp_code.emplace_back(code);
            tmp_alpha.emplace_back(alpha);
            tmp_alphaerr.emplace_back(alphaerr);
            tmp_meane.emplace_back(meane);
            tmp_meaneerr.emplace_back(meaneerr);
            tmp_toten.emplace_back(toten);
            tmp_totenerr.emplace_back(totenerr);
            
        }
        
        nin.close();
        
        code_OtherModels.emplace_back(tmp_code);
        alpha_OtherModels.emplace_back(tmp_alpha);
        alphaerr_OtherModels.emplace_back(tmp_alphaerr);
        meane_OtherModels.emplace_back(tmp_meane);
        meaneerr_OtherModels.emplace_back(tmp_meaneerr);
        toten_OtherModels.emplace_back(tmp_toten);
        totenerr_OtherModels.emplace_back(tmp_totenerr);
    }
    
    std::vector<std::vector<TGraph*>> points_AlphaVsE0;
    std::vector<std::vector<TGraph*>> points_LumVsE0;
    std::vector<std::vector<TGraph*>> points_LumVsAlpha;
    
    std::vector<std::vector<TGraphErrors*>> pointserr_AlphaVsE0;
    std::vector<std::vector<TGraphErrors*>> pointserr_LumVsE0;
    std::vector<std::vector<TGraphErrors*>> pointserr_LumVsAlpha;
    
    for(int iFile = 0; iFile < file_OtherModels.size(); ++iFile){
        
        std::vector<TGraph*> tmp_points_AlphaVsE0;
        std::vector<TGraph*> tmp_points_LumVsE0;
        std::vector<TGraph*> tmp_points_LumVsAlpha;
        
        std::vector<TGraphErrors*> tmp_pointserr_AlphaVsE0;
        std::vector<TGraphErrors*> tmp_pointserr_LumVsE0;
        std::vector<TGraphErrors*> tmp_pointserr_LumVsAlpha;
        
        for(int i = 0; i < code_OtherModels[iFile].size(); ++i){
            Double_t alphaval[1];
            Double_t meaneval[1];
            Double_t totenval[1];
            
            alphaval[0] = alpha_OtherModels[iFile][i];
            meaneval[0] = meane_OtherModels[iFile][i];
            totenval[0] = toten_OtherModels[iFile][i]/factor; //10^{53} erg
            
            Double_t alphavalerr[1];
            Double_t meanevalerr[1];
            Double_t totenvalerr[1];
            
            alphavalerr[0] = alphaerr_OtherModels[iFile][i];
            meanevalerr[0] = meaneerr_OtherModels[iFile][i];
            totenvalerr[0] = totenerr_OtherModels[iFile][i]/factor; //10^{53} erg
            
            //std::cout << alphaval[0] << " " << meaneval[0] << " " << totenval[0] << std::endl;
            
            //make TGraph object to draw
            TGraph *point_AlphaVsE0 = new TGraph(1, meaneval, alphaval);
            point_AlphaVsE0->SetMarkerStyle(20);
            point_AlphaVsE0->SetMarkerSize(1);
            
            TGraph *point_LumVsE0 = new TGraph(1, meaneval, totenval);
            point_LumVsE0->SetMarkerStyle(20);
            point_LumVsE0->SetMarkerSize(1);
            
            TGraph *point_LumVsAlpha = new TGraph(1, alphaval, totenval);
            point_LumVsAlpha->SetMarkerStyle(20);
            point_LumVsAlpha->SetMarkerColor(kRed);
            point_LumVsAlpha->SetMarkerSize(1);
            
            TGraphErrors *point_err_AlphaVsE0 = new TGraphErrors(1, meaneval, alphaval, meanevalerr, alphavalerr);
            point_err_AlphaVsE0->SetMarkerStyle(20);
            point_err_AlphaVsE0->SetMarkerColor(kRed); point_err_AlphaVsE0->SetLineColor(kRed);
            point_err_AlphaVsE0->SetMarkerSize(0.5);
            
            TGraphErrors *point_err_LumVsE0 = new TGraphErrors(1, meaneval, totenval, meanevalerr, totenvalerr);
            point_err_LumVsE0->SetMarkerStyle(20);
            point_err_LumVsE0->SetMarkerColor(kRed); point_err_LumVsE0->SetLineColor(kRed);
            point_err_LumVsE0->SetMarkerSize(0.5);
            
            TGraphErrors *point_err_LumVsAlpha = new TGraphErrors(1, alphaval, totenval, alphavalerr, totenvalerr);
            point_err_LumVsAlpha->SetMarkerStyle(20);
            point_err_LumVsAlpha->SetMarkerColor(kRed); point_err_LumVsAlpha->SetLineColor(kRed);
            point_err_LumVsAlpha->SetMarkerSize(0.5);
            
            if(iFile+2==5){
                point_AlphaVsE0->SetMarkerColor(kYellow+2);
                point_LumVsE0->SetMarkerColor(kYellow+2);
                point_LumVsAlpha->SetMarkerColor(kYellow+2);
                point_err_AlphaVsE0->SetMarkerColor(kYellow+2);
                point_err_LumVsE0->SetMarkerColor(kYellow+2);
                point_err_LumVsAlpha->SetMarkerColor(kYellow+2);
                
            }
            else if(iFile+2==3){
                point_AlphaVsE0->SetMarkerColor(kGreen+2);
                point_LumVsE0->SetMarkerColor(kGreen+2);
                point_LumVsAlpha->SetMarkerColor(kGreen+2);
                point_err_AlphaVsE0->SetMarkerColor(kGreen+2);
                point_err_LumVsE0->SetMarkerColor(kGreen+2);
                point_err_LumVsAlpha->SetMarkerColor(kGreen+2);
            }
            else{
                point_AlphaVsE0->SetMarkerColor(iFile+2);
                point_LumVsE0->SetMarkerColor(iFile+2);
                point_LumVsAlpha->SetMarkerColor(iFile+2);
                point_err_AlphaVsE0->SetMarkerColor(iFile+2);
                point_err_LumVsE0->SetMarkerColor(iFile+2);
                point_err_LumVsAlpha->SetMarkerColor(iFile+2);
            }
            
            tmp_points_AlphaVsE0.emplace_back(point_AlphaVsE0);
            tmp_points_LumVsE0.emplace_back(point_LumVsE0);
            tmp_points_LumVsAlpha.emplace_back(point_LumVsAlpha);
            
            tmp_pointserr_AlphaVsE0.emplace_back(point_err_AlphaVsE0);
            tmp_pointserr_LumVsE0.emplace_back(point_err_LumVsE0);
            tmp_pointserr_LumVsAlpha.emplace_back(point_err_LumVsAlpha);
            
        }
        
        points_AlphaVsE0.emplace_back(tmp_points_AlphaVsE0);
        points_LumVsE0.emplace_back(tmp_points_LumVsE0);
        points_LumVsAlpha.emplace_back(tmp_points_LumVsAlpha);
        
        pointserr_AlphaVsE0.emplace_back(tmp_pointserr_AlphaVsE0);
        pointserr_LumVsE0.emplace_back(tmp_pointserr_LumVsE0);
        pointserr_LumVsAlpha.emplace_back(tmp_pointserr_LumVsAlpha);
    }
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GRAB EVENT INFORMATION
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TString a("alpha");
    TString e("e0");
    TString l("lum");
    
    std::vector<double> alpha = getParameter(files[0], a);
    std::vector<double> e0 = getParameter(files[0], e);
    std::vector<double> lum = getParameter(files[0], l);
	
    Double_t tmpalpha[1];
    Double_t tmpe0[1];
    Double_t tmplum[1];
		
    tmpalpha[0] = alpha[0];
    tmpe0[0] = e0[0];
    tmplum[0] = lum[0];
		
    TString BestFitPoints_labels = "Best-Fit: #alpha = " + TString::Format("%4.1f",tmpalpha[0]) + ", #LT E_{#nu} #GT = " + TString::Format("%4.1f",tmpe0[0]) + ", #varepsilon = " + TString::Format("%4.1f",tmplum[0]) + "e53";
    
    TGraph *BestFitPoints_AlphaVsE0 = new TGraph(1, tmpe0, tmpalpha);
    BestFitPoints_AlphaVsE0->SetMarkerStyle(29);
    BestFitPoints_AlphaVsE0->SetMarkerColor(kBlack);
    BestFitPoints_AlphaVsE0->SetMarkerSize(2.0);
    
    TGraph *BestFitPoints_LumVsE0 = new TGraph(1, tmpe0, tmplum);
    BestFitPoints_LumVsE0->SetMarkerStyle(29);
    BestFitPoints_LumVsE0->SetMarkerColor(kBlack);
    BestFitPoints_LumVsE0->SetMarkerSize(2.0);
		
    TGraph *BestFitPoints_LumVsAlpha = new TGraph(1, tmpalpha, tmplum);
    BestFitPoints_LumVsAlpha->SetMarkerStyle(29);
    BestFitPoints_LumVsAlpha->SetMarkerColor(kBlack);
    BestFitPoints_LumVsAlpha->SetMarkerSize(2.0);
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GRAB CHI2 PLOTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
	std::vector<std::vector<TH2D*>> hists;
	for(int iHist = 0; iHist < histNames.size(); ++iHist){
        std::vector<TH2D*> tmp;
        for(int iFile = 0; iFile < files.size(); ++iFile){
            TFile *f = new TFile(files[iFile]);
            //TH2D *h = (TH2D*)gDirectory->Get(histNames[iHist]);
            
            TH2D *hToInterpolate = (TH2D*)gDirectory->Get(histNamesForInterpolation[iHist]);
            TH2D *h = InterpolateChi2Hist(hToInterpolate, 200, 4.61);
            
            tmp.emplace_back(h);
            
        }//end loop over files
        
        
        //okay going to try and merge relevant hists
        std::vector<TH2D*> mergedhists;
 
        TH2D *hmerged = mergeHistograms(tmp, histNamesForInterpolation[iHist]);
        hmerged->SetContour(1, contours);
        hmerged->SetLineWidth(3);
        //also try to make the labels/axis titles better
        hmerged->GetXaxis()->CenterTitle();
        hmerged->GetXaxis()->SetTitleOffset(1.0);
        hmerged->GetXaxis()->SetTitleSize(0.05);
        hmerged->GetXaxis()->SetLabelSize(0.07);
        hmerged->GetYaxis()->CenterTitle();
        hmerged->GetYaxis()->SetTitleOffset(1.0);
        hmerged->GetYaxis()->SetTitleSize(0.05);
        hmerged->GetYaxis()->SetLabelSize(0.07);
        hmerged->SetLineColor(kBlack);
            
        mergedhists.emplace_back(hmerged);
           
        hists.emplace_back(mergedhists);
	}
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DRAW PLOTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	TCanvas *c1 = new TCanvas("c1", "c1", 1000, 700);
    c1->Divide(2, 2);
    
    std::vector<TH2D*> axisHistograms;
    for(int i = 0; i < hists.size(); ++i){
        
        int toDraw = i+1;
        if(i >= 1) ++toDraw;
        //hists[i][0]
        
        TString histAxisName = "axis_" + TString::Format("%d",toDraw);
        
        double xmin;// = hists[i]->GetXaxis()->GetXmin();
        double xmax;// = hists[i]->GetXaxis()->GetXmax();
        double ymin;// = hists[i]->GetYaxis()->GetXmin();
        double ymax;// = hists[i]->GetYaxis()->GetXmax();
        
        getMinOrMaxFromHists(hists[i], xmin, xmax, ymin, ymax);
        
        TString xTitle("");
        TString yTitle("");
        
        if(toDraw == 1){ //alpha vs e0
            xmax-=0.06;
            ymin+=0.15;
            ymax-=0.15;
            xTitle = "#LT E_{#nu} #GT (MeV)";
            yTitle = "#alpha";
        }
        else if(toDraw == 3){ //lum vs e0
            xmax-=0.06;
            
            //ymin = 0.01;
            ymax -= 0.01;
            
            xTitle = "#LT E_{#nu} #GT (MeV)";
            yTitle = "#varepsilon (10^{53} erg)";
        }
        else{ //lum vs alpha
            xmin+=0.15;
            xmax-=0.15;
            
            //ymin = 0.01;
            //ymax = 1.4;
            ymax -= 0.01;
            
            xTitle = "#alpha";
            yTitle = "#varepsilon (10^{53} erg)";
        }
        
        TH2D *histAxis = new TH2D(histAxisName, " ;;", 2, xmin, xmax, 2, ymin, ymax);
        histAxis->GetXaxis()->SetTitle( xTitle );
        histAxis->GetYaxis()->SetTitle( yTitle );
        
        histAxis->GetXaxis()->CenterTitle();
        histAxis->GetXaxis()->SetTitleOffset(1.1);
        histAxis->GetYaxis()->CenterTitle();
        histAxis->GetYaxis()->SetTitleOffset(1.0);
        if(toDraw == 3 || toDraw == 4) histAxis->GetYaxis()->SetTitleOffset(1.1);
        
        axisHistograms.emplace_back(histAxis);
        
    }
    
	
    for(int i = 0; i < hists.size(); ++i){
		
		int toDraw = i+1;
		if(i >= 1) ++toDraw;
		c1->cd(toDraw);
        gPad->SetBottomMargin(0.18);
        gPad->SetTopMargin(0.03);
        gPad->SetRightMargin(0.03);
        gPad->SetLeftMargin(0.14);
        
        for(int j = 0; j < hists[i].size(); ++j){
            //if(j+1==5) hists[i][j]->SetLineColor(j+2);
            //else if(j+1==3) hists[i][j]->SetLineColor(kGreen+2);
            //else hists[i][j]->SetLineColor(j+1);
            
            hists[i][j]->SetLineColor(kBlack);
             
            if(j == 0){
                axisHistograms[i]->GetXaxis()->SetLabelSize(0.065);
                axisHistograms[i]->GetYaxis()->SetLabelSize(0.065);
                axisHistograms[i]->GetXaxis()->SetTitleSize(0.065);
                axisHistograms[i]->GetYaxis()->SetTitleSize(0.065);
                axisHistograms[i]->GetXaxis()->SetTitleOffset(1.3);
                axisHistograms[i]->GetYaxis()->SetTitleOffset(1.1);
                axisHistograms[i]->Draw();
            }
            
            hists[i][j]->Draw("cont3 same");
            gPad->Update();
            
        }
		
	}
	
	//also draw the truth points
    
	c1->cd(1);
    for(int i = 0; i < points_AlphaVsE0.size(); ++i)
        for(int j = 0; j < points_AlphaVsE0[i].size(); ++j) points_AlphaVsE0[i][j]->Draw("P SAME");
    BestFitPoints_AlphaVsE0->Draw("P SAME");
	c1->cd(3);
    for(int i = 0; i < points_LumVsE0.size(); ++i)
        for(int j = 0; j < points_LumVsE0[i].size(); ++j) points_LumVsE0[i][j]->Draw("P SAME");
    BestFitPoints_LumVsE0->Draw("P SAME");
	c1->cd(4);
    for(int i = 0; i < points_LumVsAlpha.size(); ++i)
        for(int j = 0; j < points_LumVsAlpha[i].size(); ++j) points_LumVsAlpha[i][j]->Draw("P SAME");
    BestFitPoints_LumVsAlpha->Draw("P SAME");
    
    c1->cd(2);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.01);
    gPad->SetRightMargin(0.015);
    TLegend *leg = makeLegend(0.1,0.1,0.99,0.99);
    leg->SetTextSize(0.07);
    leg->SetHeader(legendTitle,"C");
    leg->AddEntry(BestFitPoints_AlphaVsE0, "Truth: #LT E_{#nu} #GT = 9.5 MeV, #varepsilon = 5e52 erg", "p");
    for(int i = 0; i < names_config.size(); ++i) leg->AddEntry(hists[0][i], names_config[i], "l");
    //leg->AddEntry(hists[0][0], name);
    for(int i = 0; i < otherModelLegNames.size(); ++i) leg->AddEntry(points_AlphaVsE0[i][0], otherModelLegNames[i], "p");
    leg->Draw();
    
    
    
}//end of code
