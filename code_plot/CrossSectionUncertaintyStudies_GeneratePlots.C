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

void CrossSectionUncertaintyStudies_GeneratePlots(){
    
    // ********************************
    // HERE ARE THE VARIABLES WE CHANGE
    // ********************************

    //define the smearing we used
    TString smear("0.10MARLEY_NO_nueOnly");
    
    //number of entries
    const Int_t numEntries = 12;
    
    std::vector<TString> grid_xscn; std::vector<double> grid_axis;
    grid_xscn.emplace_back("Bhattacharya2009Xscn-50Percent"); grid_axis.emplace_back(0);
    grid_xscn.emplace_back("Bhattacharya2009Xscn-20Percent"); grid_axis.emplace_back(1);
    grid_xscn.emplace_back("Bhattacharya2009Xscn-15Percent"); grid_axis.emplace_back(2);
    grid_xscn.emplace_back("Bhattacharya2009Xscn-10Percent"); grid_axis.emplace_back(3);
    grid_xscn.emplace_back("Bhattacharya2009Xscn-5Percent"); grid_axis.emplace_back(4);
    grid_xscn.emplace_back("Bhattacharya2009Xscn"); grid_axis.emplace_back(5);
    grid_xscn.emplace_back("Bhattacharya2009Xscn+5Percent"); grid_axis.emplace_back(6);
    grid_xscn.emplace_back("Bhattacharya2009Xscn+10Percent"); grid_axis.emplace_back(7);
    grid_xscn.emplace_back("Bhattacharya2009Xscn+15Percent"); grid_axis.emplace_back(8);
    grid_xscn.emplace_back("Bhattacharya2009Xscn+20Percent"); grid_axis.emplace_back(9);
    grid_xscn.emplace_back("Bhattacharya2009Xscn+50Percent"); grid_axis.emplace_back(10);
    grid_xscn.emplace_back("Bhattacharya2009Xscn+100Percent"); grid_axis.emplace_back(11);
    
    std::vector<TString> ts_xscn; std::vector<double> spectra_axis;
    ts_xscn.emplace_back("Bhattacharya2009Xscn-50Percent"); spectra_axis.emplace_back(0);
    ts_xscn.emplace_back("Bhattacharya2009Xscn-20Percent"); spectra_axis.emplace_back(1);
    ts_xscn.emplace_back("Bhattacharya2009Xscn-15Percent"); spectra_axis.emplace_back(2);
    ts_xscn.emplace_back("Bhattacharya2009Xscn-10Percent"); spectra_axis.emplace_back(3);
    ts_xscn.emplace_back("Bhattacharya2009Xscn-5Percent"); spectra_axis.emplace_back(4);
    ts_xscn.emplace_back("Bhattacharya2009Xscn"); spectra_axis.emplace_back(5);
    ts_xscn.emplace_back("Bhattacharya2009Xscn+5Percent"); spectra_axis.emplace_back(6);
    ts_xscn.emplace_back("Bhattacharya2009Xscn+10Percent"); spectra_axis.emplace_back(7);
    ts_xscn.emplace_back("Bhattacharya2009Xscn+15Percent"); spectra_axis.emplace_back(8);
    ts_xscn.emplace_back("Bhattacharya2009Xscn+20Percent"); spectra_axis.emplace_back(9);
    ts_xscn.emplace_back("Bhattacharya2009Xscn+50Percent"); spectra_axis.emplace_back(10);
    ts_xscn.emplace_back("Bhattacharya2009Xscn+100Percent"); spectra_axis.emplace_back(11);
    
    const char *labels[numEntries] = {"-50%", "-20%", "-15%", "-10%", "-5%", "0%", "+5%", "+10%", "+15%", "+20%", "+50%", "+100%"};
    std::vector<double> labelVals = {-50.0, -20.0, -15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 50.0, 100.0};
    
    //define distance
    TString distanceSN("10.00kpc");
    
	std::vector<TString> histNames;
	histNames.emplace_back("alphavse0_90cut");
	histNames.emplace_back("lumvse0_90cut");
	histNames.emplace_back("lumvsalpha_90cut");
    
    std::vector<TString> fancyHistNames;
    fancyHistNames.emplace_back("#alpha vs. #LT E_{#nu} #GT");
    fancyHistNames.emplace_back("#varepsilon vs. #LT E_{#nu} #GT");
    fancyHistNames.emplace_back("#varepsilon vs. #alpha");
	
	std::vector<TString> graphNames;
	graphNames.emplace_back("Chi2VsAlpha");
	graphNames.emplace_back("Chi2VsE0");
	graphNames.emplace_back("Chi2VsLum");
	
	std::vector<TString> parNames;
	parNames.emplace_back("#alpha");
	parNames.emplace_back("#LT E_{#nu} #GT");
	parNames.emplace_back("#varepsilon");
	
    //define the grid we want to use
    TString grid("2019October17");
    //TString grid("2018October25");
    
    //define the efficiency we want to use
    TString effic("StepEffic");
    
	//finally, define steps for distance and grid_axis (for binning)
	Double_t step_res = 5.0;
	Double_t step_ts = 5.0;
	
    // ********************************
    // DEFINE INPUTS
    // ********************************
    //define file names   
    std::vector<TString> files;
    std::vector<TString> filenames;
    std::vector<std::vector<double>> pairs;
    
    TString filename_for_scaling;
    
    for(size_t i = 0; i < grid_xscn.size(); ++i){
		
		for(size_t j = 0; j < ts_xscn.size(); ++j){
				//chi2plots_" + grid + "_smear" + dc + "_spectra" + ts + "_" + xscn + "_" + distanceStr + "kpc.root"
			TString tmp = "out/chi2plots_" + grid + "_smear" + smear + "_" + grid_xscn[i] + "_" + effic + "_spectra" + smear + "_" + ts_xscn[j] + "_" + effic + "_" + distanceSN + ".root";
            
			files.emplace_back(tmp);
            
            if(i==5 && j == 5) filename_for_scaling = tmp;
			
			TString tmpname = grid + "_smear" + smear + "_" + grid_xscn[i] + "_" + effic + "_spectra" + smear + "_" + ts_xscn[j] + "_" + effic;
			filenames.emplace_back(tmpname);
			
			//std::cout << tmpname << std::endl;
			//std::cout << tmp << " " << grid_axis[i] << " " << spectra_axis[j] << std::endl;
			
			std::vector<double> tmppair;
			tmppair.emplace_back(grid_axis[i]);
			tmppair.emplace_back(spectra_axis[j]);
			pairs.emplace_back(tmppair);
						
		}	
		
	}
    
	//set specific aesthetic things here
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1,0);
	gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(99);
    //gStyle->SetTitleFontSize(0.07);
    gStyle->SetTitleSize(0.05);
	
	//define contours
	double contours[1];
	contours[0] = 0.001;
	
	//define binning, axis information for 2d plots
	Double_t first_ts = *std::min_element(spectra_axis.begin(), spectra_axis.end());
	Double_t last_ts = *std::max_element(spectra_axis.begin(), spectra_axis.end());
	Int_t numtsbins = int( last_ts - first_ts )/step_ts+1;
	Double_t first_res = *std::min_element(grid_axis.begin(), grid_axis.end());
	Double_t last_res = *std::max_element(grid_axis.begin(), grid_axis.end());
	Int_t numresbins = int( last_res - first_res )/step_res+1;
	
	Double_t tsmin = first_ts - step_ts/2.0;
	Double_t tsmax = last_ts + step_ts/2.0;
	Double_t resmin = first_res - step_res/2.0;
	Double_t resmax = last_res + step_res/2.0;
	
	//truth points
	Double_t alpha_true = 2.5;
	Double_t e0_true = 9.5; //mev
	Double_t lum_true = 0.5; //10^{53} ergs
    
    // ********************************
    // GRID INFORMATION
    // ********************************
    TString gridinfo = "input/grid_info/grid_info_" + grid + ".dat";
    
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
    
    //define biases
    Double_t minbias_alpha = (first_alpha-alpha_true)/alpha_true;
    Double_t maxbias_alpha = (last_alpha-alpha_true)/alpha_true;
    
    Double_t minbias_e0 = (first_e0-e0_true)/e0_true;
    Double_t maxbias_e0 = (last_e0-e0_true)/e0_true;
    
    Double_t minbias_lum = (first_lum/factor-lum_true)/lum_true;
    Double_t maxbias_lum = (last_lum/factor-lum_true)/lum_true;
    
    std::cout << "Color-scale boundaries used: " << std::endl;
    std::cout <<"\talpha: (" << minbias_alpha << ", " << maxbias_alpha << ")" << std::endl;
    std::cout <<"\tE0: (" << minbias_e0 << ", " << maxbias_e0 << ")" << std::endl;
    std::cout <<"\tlum: (" << minbias_lum << ", " << maxbias_lum << ")" << std::endl;
    
    // ********************************
    // GRAB EVENT INFORMATION
    // ********************************
    TString a("alpha");
    TString e("e0");
    TString l("lum");
    TString ch("chi2");
    TString d("dof");
    TString numgrid_ave("numgrid_alphavse0");
    TString numgrid_lve("numgrid_lumvse0");
    TString numgrid_lva("numgrid_lumvsalpha");
    TString numgood_ave("numgood_alphavse0");
    TString numgood_lve("numgood_lumvse0");
    TString numgood_lva("numgood_lumvsalpha");
    
    std::vector<std::vector<double>> alpha = getManyParameters(files, a);
    std::vector<std::vector<double>> e0 = getManyParameters(files, e);
    std::vector<std::vector<double>> lum = getManyParameters(files, l);
    std::vector<std::vector<double>> chi2 = getManyParameters(files, ch);
    std::vector<std::vector<double>> dof = getManyParameters(files, d);
    std::vector<std::vector<double>> numgrid_alphavse0 = getManyParameters(files, numgrid_ave);
    std::vector<std::vector<double>> numgrid_lumvse0 = getManyParameters(files, numgrid_lve);
    std::vector<std::vector<double>> numgrid_lumvsalpha = getManyParameters(files, numgrid_lva);
    std::vector<std::vector<double>> numgood_alphavse0 = getManyParameters(files, numgood_ave);
    std::vector<std::vector<double>> numgood_lumvse0 = getManyParameters(files, numgood_lve);
    std::vector<std::vector<double>> numgood_lumvsalpha = getManyParameters(files, numgood_lva);
    
    // ********************************
    // SCALING FACTORS
    // ********************************
    std::vector<double> scalingFactors;
    std::vector<double> tmp_ngrid_ave = getParameter(filename_for_scaling, numgrid_ave);
    std::vector<double> tmp_ngrid_lve = getParameter(filename_for_scaling, numgrid_lve);
    std::vector<double> tmp_ngrid_lva = getParameter(filename_for_scaling, numgrid_lva);
    std::vector<double> tmp_ngood_ave = getParameter(filename_for_scaling, numgood_ave);
    std::vector<double> tmp_ngood_lve = getParameter(filename_for_scaling, numgood_lve);
    std::vector<double> tmp_ngood_lva = getParameter(filename_for_scaling, numgood_lva);
    
    scalingFactors.emplace_back(tmp_ngood_ave[0]/tmp_ngrid_ave[0]);
    scalingFactors.emplace_back(tmp_ngood_lve[0]/tmp_ngrid_lve[0]);
    scalingFactors.emplace_back(tmp_ngood_lva[0]/tmp_ngrid_lva[0]);
    
    // ********************************
    // MAKE 2D HISTS
    // ********************************
    
	std::vector<TH2D*> hists_absdiff;
    std::vector<TH2D*> hists_fracelements;
    
    //so we can loop through pairs assuming that
    //pairs.size() == alpha.size() == number of files
    
    //loop over the parameters
    for(Int_t iPar = 0; iPar < parNames.size(); ++iPar){
        
        TString titlediff = "Fractional difference from truth for " + parNames[iPar];
        TString namediff = "fracdiff_" + parNames[iPar];
        
        //TH2D *hdiff = new TH2D(namediff, titlediff, numresbins, resmin, resmax, numtsbins, tsmin, tsmax);
        TH2D *hdiff = new TH2D(namediff, titlediff, numEntries, 0, numEntries-1, numEntries, 0, numEntries-1);
        //TH2D *hdiff = new TH2D(namediff, titlediff, 11, 0, 10, 11, 0, 10);
        
        TString titlenum = "Fractional area: " + fancyHistNames[iPar];
        TString namenum = "fracNumElements_" + fancyHistNames[iPar];
        //TH2D *hnum = new TH2D(namenum, titlenum, numEntries, 0, numEntries-1, numEntries, 0, numEntries-1);
        TH2D *hnum = new TH2D(namenum, titlenum, 11, 0, 10, 11, 0, 10);
        
        Double_t smallestDiff = 0.;
        
        
        //now loop through the pairs
        for(Int_t i = 0; i < pairs.size(); ++i){
            //pairs[i][0] = grid_axis
            //pairs[i][1] = spectra
            //now we fill the hist:
            
            double diffToFill(0.0);
            double fracToFill(0.0);
            
            if(iPar == 0){
                diffToFill = (alpha[i][0] - alpha_true)/alpha_true;
                //std::cout << grid_xscn[(Int_t)pairs[i][0]] << " " << ts_xscn[(Int_t)pairs[i][1]] << " " << alpha[i][0] << std::endl;
                
                //alphavse0
                fracToFill = numgood_alphavse0[i][0]/numgrid_alphavse0[i][0];
                
            }
            else if(iPar == 1){
                diffToFill = (e0[i][0] - e0_true)/e0_true;
                
                //lumvse0
                fracToFill = numgood_lumvse0[i][0]/numgrid_lumvse0[i][0];
            }
            else{
                diffToFill = (lum[i][0] - lum_true)/lum_true;
                //double fracDiff = (lum[i][0] - lum_true)/lum_true;
                //if(fracDiff < 0) diffToFill = -log(1+std::abs(fracDiff) );
                //else diffToFill =log(1+std::abs(fracDiff) );
                
                //std::cout << grid_xscn[(Int_t)pairs[i][0]] << " " << ts_xscn[(Int_t)pairs[i][1]] << " " << fracDiff << std::endl;
                
                //lumvsalpha
                fracToFill = numgood_lumvsalpha[i][0]/numgrid_lumvsalpha[i][0];
                
            }
            
            //std::cout << grid_xscn[(Int_t)pairs[i][0]] << " " << ts_xscn[(Int_t)pairs[i][1]] << " " << diffToFill << std::endl;
            
            //if(diffToFill == 0.0) diffToFill = 1e-20;
            
            Int_t index1 = (Int_t)pairs[i][0];
            Int_t index2 = (Int_t)pairs[i][1];
            hdiff->Fill(labels[index1], labels[index2], diffToFill);
            
            fracToFill /= scalingFactors[iPar];
            hnum->Fill(labels[index1], labels[index2], fracToFill);
            
            //std::cout << labels[index1] << " " << labels[index2] << " " << fancyHistNames[iPar] << " " << fracToFill << std::endl;
            
            //hdiff->Fill(pairs[i][0], pairs[i][1], diffToFill);
            
            //std::cout << parNames[iPar] << " " << grid_xscn[index1] << " " << ts_xscn[index2] << " " << diffToFill << std::endl;
            
            //std::cout << parNames[iPar] << ": (" << labels[index1] << ", " << labels[index2] << ") " << diffToFill << std::endl;
            
            //std::cout << fancyHistNames[iPar] << " " << grid_xscn[index1] << " " << ts_xscn[index2] << " " << fracToFill << std::endl;
            
            if(diffToFill < smallestDiff) smallestDiff = diffToFill;
            
        }
        
        hdiff->SetMinimum(smallestDiff-1e-20);
        hdiff->GetXaxis()->SetLabelSize(0.07);
        hdiff->GetXaxis()->SetLabelOffset(0.01);
        hdiff->GetYaxis()->SetLabelSize(0.07);
        hdiff->GetZaxis()->SetLabelSize(0.06);
        hdiff->LabelsDeflate("X");
        hdiff->LabelsDeflate("Y");
        hdiff->GetYaxis()->LabelsOption("v");
        hdiff->GetXaxis()->LabelsOption("v");
        hdiff->GetXaxis()->SetTickLength(0.);
        hdiff->GetYaxis()->SetTickLength(0.);
        hdiff->SetTitleSize(0.09);
        
        
        hnum->SetMinimum(-1e-20);
        hnum->GetXaxis()->SetLabelSize(0.09);
        hnum->GetXaxis()->SetLabelOffset(0.01);
        hnum->GetYaxis()->SetLabelSize(0.09);
        hnum->GetZaxis()->SetLabelSize(0.065);
        hnum->LabelsDeflate("X");
        hnum->LabelsDeflate("Y");
        hnum->GetYaxis()->LabelsOption("v");
        hnum->GetXaxis()->LabelsOption("v");
        hnum->GetXaxis()->SetTickLength(0.);
        hnum->GetYaxis()->SetTickLength(0.);
        hnum->SetTitleSize(0.09);
        
        
        //hdiff->GetYaxis()->RotateTitle(1);
        //hdiff->GetYaxis()->LabelsOption("v");
        
        //now we save the hist by doing
        hists_absdiff.emplace_back(hdiff);
        hists_fracelements.emplace_back(hnum);
        
    }
	
	
    
    // ********************************
    // DRAW PLOTS
    // ********************************
    TCanvas *c3 = new TCanvas("c3","c3",1100,700);
    c3->Divide(2, 2);
    c3->cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.13);
    gPad->SetBottomMargin(0.19);
    gStyle->SetTitleFontSize(0.08);
    hists_absdiff[0]->GetZaxis()->SetRangeUser(minbias_alpha, maxbias_alpha);
    SetColorScale(hists_absdiff[0]);
    hists_absdiff[0]->Draw("colz");
    gPad->Update();
    c3->cd(3);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.13);
    gPad->SetBottomMargin(0.19);
    gStyle->SetTitleFontSize(0.08);
    hists_absdiff[1]->GetZaxis()->SetRangeUser(minbias_e0, maxbias_e0);
    SetColorScale(hists_absdiff[1]);
    hists_absdiff[1]->Draw("colz");
    gPad->Update();
    c3->cd(4);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.13);
    gPad->SetBottomMargin(0.19);
    gStyle->SetTitleFontSize(0.08);
    hists_absdiff[2]->GetZaxis()->SetRangeUser(minbias_lum, maxbias_lum);
    SetColorScale(hists_absdiff[2]);
    hists_absdiff[2]->Draw("colz");
    gPad->Modified();
    gPad->Update();
    
    /*
    c3 = new TCanvas("c4","c4",400,1000);
    c3->Divide(1, 3);
    c3->cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.13);
    gPad->SetBottomMargin(0.19);
    hists_fracelements[0]->Draw("colz");
    gPad->Update();
    c3->cd(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.13);
    gPad->SetBottomMargin(0.19);
    hists_fracelements[1]->Draw("colz");
    c3->cd(3);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.13);
    gPad->SetBottomMargin(0.19);
    //gPad->SetLogz();
    //hists_absdiff[2]->GetZaxis()->SetRangeUser(-1.0, 4.0);
    //hists_absdiff[2]->GetZaxis()->SetNdivisions(6);
    
    hists_fracelements[2]->Draw("colz");
    gPad->Modified();
    */
    
}//end of code
