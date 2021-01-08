//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// AsimovSmearingStudies_GeneratePlots
// Erin Conley (erin.conley@duke.edu)
// Description: For given test spectrum parameter assumptions and scaling
//              factors, creates FOM plots for different combinations of
//              SNOwGLoBES smearing (used to produce the grid and test specta).
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

void AsimovSmearingStudies_GeneratePlots(){
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // HERE ARE THE VARIABLES WE CHANGE
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
	//define the grid we want to use
    TString grid("2019October17");
    
    //define the xscn model we want to use
    TString xscn("PNXscn");
    
    //define the efficiency we want to use
    TString effic("StepEffic");
    
    //updated plotting with char labels
    //example: energy scaling
    //change config and labels to suit the study
    std::vector<TString> grid_smear_config; std::vector<double> grid_smear_num;
    grid_smear_config.emplace_back("MARLEY_DCR_NO_nueOnly-15Percent"); grid_smear_num.emplace_back(0);
    grid_smear_config.emplace_back("MARLEY_DCR_NO_nueOnly-10Percent"); grid_smear_num.emplace_back(1);
    grid_smear_config.emplace_back("MARLEY_DCR_NO_nueOnly-5Percent"); grid_smear_num.emplace_back(2);
    grid_smear_config.emplace_back("MARLEY_DCR_NO_nueOnly-1Percent"); grid_smear_num.emplace_back(3);
    grid_smear_config.emplace_back("MARLEY_DCR_NO_nueOnly"); grid_smear_num.emplace_back(4);
    grid_smear_config.emplace_back("MARLEY_DCR_NO_nueOnly+1Percent"); grid_smear_num.emplace_back(5);
    grid_smear_config.emplace_back("MARLEY_DCR_NO_nueOnly+5Percent"); grid_smear_num.emplace_back(6);
    grid_smear_config.emplace_back("MARLEY_DCR_NO_nueOnly+10Percent"); grid_smear_num.emplace_back(7);
    grid_smear_config.emplace_back("MARLEY_DCR_NO_nueOnly+15Percent"); grid_smear_num.emplace_back(8);
    
    //define the test spectrum we want to use
    std::vector<TString> test_smear_config; std::vector<double> test_smear_num;
    test_smear_config.emplace_back("MARLEY_DCR_NO_nueOnly-15Percent"); test_smear_num.emplace_back(0);
    test_smear_config.emplace_back("MARLEY_DCR_NO_nueOnly-10Percent"); test_smear_num.emplace_back(1);
    test_smear_config.emplace_back("MARLEY_DCR_NO_nueOnly-5Percent"); test_smear_num.emplace_back(2);
    test_smear_config.emplace_back("MARLEY_DCR_NO_nueOnly-1Percent"); test_smear_num.emplace_back(3);
    test_smear_config.emplace_back("MARLEY_DCR_NO_nueOnly"); test_smear_num.emplace_back(4);
    test_smear_config.emplace_back("MARLEY_DCR_NO_nueOnly+1Percent"); test_smear_num.emplace_back(5);
    test_smear_config.emplace_back("MARLEY_DCR_NO_nueOnly+5Percent"); test_smear_num.emplace_back(6);
    test_smear_config.emplace_back("MARLEY_DCR_NO_nueOnly+10Percent"); test_smear_num.emplace_back(7);
    test_smear_config.emplace_back("MARLEY_DCR_NO_nueOnly+15Percent"); test_smear_num.emplace_back(8);
    
    const char *labels[9] = {"-15%", "-10%", "-5%", "-1%", "0%", "1%", "5%", "10%", "15%"};
    std::vector<double> labelVals = {-15.0, -10.0, -5.0, -1.0, 0.0, 1.0, 5.0, 10.0, 15.0};
    
    //define distance
    TString distanceSN("10.00kpc"); Double_t distance(10.0);
    
	std::vector<TString> histNames;
	histNames.emplace_back("alphavse0_90cut");
	histNames.emplace_back("lumvse0_90cut");
	histNames.emplace_back("lumvsalpha_90cut");  
	
	std::vector<TString> graphNames;
	graphNames.emplace_back("Chi2VsAlpha");
	graphNames.emplace_back("Chi2VsE0");
	graphNames.emplace_back("Chi2VsLum");
	
	std::vector<TString> parNames;
	parNames.emplace_back("#alpha");
	parNames.emplace_back("#LT E_{#nu} #GT");
	parNames.emplace_back("#varepsilon");

	//define extra directory where the files live, e.g., "/CrossSectionStudy/"
    //ps: if there isn't an extra directory and they live in "out", set the parameter to "/"
    TString dir("/");	

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE INPUTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //define file names   
    std::vector<TString> files;
    std::vector<TString> filenames;
    std::vector<Double_t> refareas_alphavse0;
    std::vector<Double_t> refareas_lumvse0;
    std::vector<Double_t> refareas_lumvsalpha;
    std::vector<std::vector<double>> pairs;
    for(size_t i = 0; i < grid_smear_config.size(); ++i){
		
		for(size_t j = 0; j < test_smear_config.size(); ++j){
			
			TString tmp = "out" + dir +"chi2plots_" + grid + "_smear" + grid_smear_config[i] + "_" + xscn + "_" + effic + "_spectra" + test_smear_config[j] + "_" + xscn + "_" + effic + "_" + distanceSN + ".root";
            
			files.emplace_back(tmp);
			
			TString tmpname = grid + "_smear" + grid_smear_config[i] + "_" + xscn + "_" + effic + "_spectra" + test_smear_config[j] + "_" + xscn + "_" + effic;
			filenames.emplace_back(tmpname);
            
            Double_t ref_alphavse0; Double_t ref_lumvse0; Double_t ref_lumvsalpha;
            
            //getReferenceContourAreas(distance, grid, xscn, effic, xscn, effic, &ref_alphavse0, &ref_lumvse0, &ref_lumvsalpha);
            
            refareas_alphavse0.emplace_back(ref_alphavse0);
            refareas_lumvse0.emplace_back(ref_lumvse0);
            refareas_lumvsalpha.emplace_back(ref_lumvsalpha);
			
			//std::cout << tmpname << std::endl;
			//std::cout << tmp << " " << grid_smear_num[i] << " " << test_smear_num[j] << std::endl;
			
			std::vector<double> tmppair;
			tmppair.emplace_back(grid_smear_num[i]);
			tmppair.emplace_back(test_smear_num[j]);
			pairs.emplace_back(tmppair);
						
		}	
		
	}
    
    std::vector<std::vector<Double_t>> ReferenceAreas;
    ReferenceAreas.emplace_back(refareas_alphavse0);
    ReferenceAreas.emplace_back(refareas_lumvse0);
    ReferenceAreas.emplace_back(refareas_lumvsalpha);
	
	//define contours
	double contours[1];
	contours[0] = 0.001;
	
	//define binning, axis information for 2d plots
	Double_t *binsSmear = new Double_t[grid_smear_num.size()];
	for(size_t i = 0; i < grid_smear_num.size(); ++i) binsSmear[i] = grid_smear_num[i];
	Int_t numsmearbins = grid_smear_num.size()-1;
	
	Double_t *binsTS = new Double_t[test_smear_num.size()];
	for(size_t i = 0; i < test_smear_num.size(); ++i) binsTS[i] = test_smear_num[i];
	Int_t numtsbins = test_smear_num.size()-1;
	
	//truth points
	Double_t alpha_true = 2.5;
	Double_t e0_true = 9.5; //mev
	Double_t lum_true = 0.5; //10^{53} ergs
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GRAB EVENT INFORMATION
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::vector<std::vector<TH2D*>> hists;
	std::vector<std::vector<Double_t>> areas;
	std::vector<Double_t> areas2;
	
	for(int iHist = 0; iHist < histNames.size(); ++iHist){  
		std::vector<TH2D*> tmp;
		std::vector<Double_t> tmp2;
		for(int iFile = 0; iFile < files.size(); ++iFile){
			TFile *f = new TFile(files[iFile]);
			TH2D *h = (TH2D*)gDirectory->Get(histNames[iHist]);
			
			Double_t a = findAreaOfContour(h);
			
			if(a == -99.0) a = 0.0; 
						
			tmp.emplace_back(h);
			tmp2.emplace_back(a);
			
			areas2.emplace_back(a);
			
		}
		
		//std::cout << "-------------------" << std::endl;
		hists.emplace_back(tmp);
		areas.emplace_back(tmp2);
	}
	
	
	std::vector<std::vector<TGraph*>> graphs;
	
	for(int iGraph = 0; iGraph < graphNames.size(); ++iGraph){  
		std::vector<TGraph*> tmp;
		for(int iFile = 0; iFile < files.size(); ++iFile){
			TFile *f = new TFile(files[iFile]);
			TGraph *g = (TGraph*)gDirectory->Get(graphNames[iGraph]);
					
			tmp.emplace_back(g);
			
			//std::cout << files[iFile] << " " << graphNames[iHist] << std::endl;
			
		}
		
		//std::cout << "-------------------" << std::endl;
		graphs.emplace_back(tmp);
		
	}
	
    
    TString a("alpha");
    TString e("e0");
    TString l("lum");
    TString ch("chi2");
    TString d("dof");
    
    std::vector<std::vector<double>> alpha = getManyParameters(files, a);
    std::vector<std::vector<double>> e0 = getManyParameters(files, e);
    std::vector<std::vector<double>> lum = getManyParameters(files, l);
    std::vector<std::vector<double>> chi2 = getManyParameters(files, ch);
    std::vector<std::vector<double>> dof = getManyParameters(files, d);
	
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAKE 2D HISTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
	std::vector<TH2D*> hists_area;
    std::vector<TH2D*> hists_area_scaled;
	
	//so we can loop through pairs assuming that
	//pairs.size() == areas[i].size()
	
	//loop over the contour plots
	for(Int_t iPlot = 0; iPlot < areas.size(); ++iPlot){
		
		//histNames[iPlot] will give the appropriate histName
		//use this to define name, title of 2D hist 
        TString title = "Contour Areas for " + histNames[iPlot] ;//+ " ;Grid Energy Spectra Scaling (Percent);Test Energy Spectra Scaling (Percent)";
		TString name = "areas_" + histNames[iPlot];
		
		//TH2D *h = new TH2D(name, title, numresbins, resmin, resmax, numtsbins, tsmin, tsmax);
		//TH2D *h = new TH2D(name, title, numsmearbins, binsSmear, numtsbins, binsTS);
		TH2D *h = new TH2D(name, title, 9, 0, 8, 9, 0, 8);
        
        TString title1 = "Normalized Contour Areas for " + histNames[iPlot] ;//+ " ;Grid Energy Spectra Scaling (Percent);Test Energy Spectra Scaling (Percent)";
        TString name1 = "areas_scaled_" + histNames[iPlot];
        
        //TH2D *h1 = new TH2D(name1, title1, numdistbins, distmin, distmax, numresbins, resmin, resmax);
        TH2D *h1 = new TH2D(name1, title1, 9, 0, 8, 9, 0, 8);
        
		//now loop through the pairs
		for(Int_t i = 0; i < pairs.size(); ++i){
			//pairs[i][0] = grid_smear_num
			//pairs[i][1] = test_smear_num
			//now we fill the hist:
			//h->Fill(pairs[i][0], pairs[i][1], areas[iPlot][i]);
			
			Int_t index1 = (Int_t)pairs[i][0];
			Int_t index2 = (Int_t)pairs[i][1];
			h->Fill(labels[index1], labels[index2], areas[iPlot][i]);
			h1->Fill(labels[index1], labels[index2], areas[iPlot][i]/ReferenceAreas[iPlot][i]);
			
		}
		
		h->SetMinimum(0.0);
		h->GetXaxis()->SetLabelSize(0.06);
		h->GetYaxis()->SetLabelSize(0.06);
		
		h->GetXaxis()->SetTitleSize(0.04);
		h->GetYaxis()->SetTitleSize(0.04);
		
		h->GetXaxis()->SetTickLength(0.);
		h->GetYaxis()->SetTickLength(0.);
        
        h->GetXaxis()->SetLabelSize(0.08);
        h->GetYaxis()->SetLabelSize(0.08);
        h->GetZaxis()->SetLabelSize(0.05);
        
        //h1->SetMinimum(0.0);
        h1->GetXaxis()->SetLabelSize(0.08);
        h1->GetYaxis()->SetLabelSize(0.08);
        h1->GetZaxis()->SetLabelSize(0.05);
				
		//now we save the hist by doing
		hists_area.emplace_back(h);
		hists_area_scaled.emplace_back(h1);
	}
	
    
	std::vector<TH2D*> hists_absdiff;
	
	//so we can loop through pairs assuming that
	//pairs.size() == alpha.size() == number of files
	
	//loop over the parameters
	for(Int_t iPar = 0; iPar < parNames.size(); ++iPar){
		
		TString titlediff = "Fractional difference from truth for " + parNames[iPar] ;//+ " ;Grid Energy Spectra Scaling (Percent);Test Energy Spectra Scaling (Percent)";
		TString namediff = "fracdiff_" + parNames[iPar];
		
		//TH2D *hdiff = new TH2D(namediff, titlediff, numresbins, resmin, resmax, numtsbins, tsmin, tsmax);
		//TH2D *hdiff = new TH2D(namediff, titlediff, numsmearbins, binsSmear, numtsbins, binsTS);
		TH2D *hdiff = new TH2D(namediff, titlediff, 9, 0, 8, 9, 0, 8);
		
		Double_t smallestDiff = 0.;
		
		//now loop through the pairs
		for(Int_t i = 0; i < pairs.size(); ++i){
			//pairs[i][0] = grid_smear_num
			//pairs[i][1] = test_smear_num
			//now we fill the hist:
			
			double parameterToFill(0.0);
			double diffToFill(0.0);
			
			if(iPar == 0) diffToFill = (alpha[i][0] - alpha_true)/alpha_true;
			else if(iPar == 1) diffToFill = (e0[i][0] - e0_true)/e0_true;
			else diffToFill = (lum[i][0] - lum_true)/lum_true;
			
			//if(diffToFill == 0.0) diffToFill = 1e-20;
			
			//hdiff->Fill(pairs[i][0], pairs[i][1], diffToFill);
			
			Int_t index1 = (Int_t)pairs[i][0];
			Int_t index2 = (Int_t)pairs[i][1];
			hdiff->Fill(labels[index1], labels[index2], diffToFill);
			
			if(diffToFill < smallestDiff) smallestDiff = diffToFill;
			
            //if(abs(diffToFill) <= 0.15) std::cout << parNames[iPar] << ": (" << labels[index1] << ", " << labels[index2] << ") " << diffToFill << std::endl;
            
            //std::cout << parNames[iPar] << ": (" << labels[index1] << ", " << labels[index2] << ") " << diffToFill << std::endl;
            
		}
		
				
		hdiff->SetMinimum(smallestDiff-1e-20);
		hdiff->GetXaxis()->SetLabelSize(0.08);
		hdiff->GetYaxis()->SetLabelSize(0.08);
        hdiff->GetZaxis()->SetLabelSize(0.055);
		
		hdiff->GetXaxis()->SetTitleSize(0.04);
		hdiff->GetYaxis()->SetTitleSize(0.04);
		
		hdiff->GetXaxis()->SetTickLength(0.);
		hdiff->GetYaxis()->SetTickLength(0.);
		
		//now we save the hist by doing
		hists_absdiff.emplace_back(hdiff);
		
	}
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DRAW 2D PLOTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	TCanvas *c1 = new TCanvas("c1","c1",1100,700);
	c1->Divide(2, 2);
	c1->cd(1);
	hists_area[0]->SetTitle("Contour Areas for #alpha vs. #LT E_{#nu} #GT");
	hists_area[0]->Draw("colz");
	c1->cd(3);
	hists_area[1]->SetTitle("Contour Areas for #epsilon vs. #LT E_{#nu} #GT");
	hists_area[1]->Draw("colz");
	c1->cd(4);
	hists_area[2]->SetTitle("Contour Areas for #epsilon vs. #alpha");
	hists_area[2]->Draw("colz");
	
	//print to keep
	//c1->Print(pngout_2d);
    c1 = new TCanvas("c2","c2",1100,700);
    c1->Divide(2, 2);
    c1->cd(1);
    hists_area_scaled[0]->SetTitle("Normalized Contour Areas for #alpha vs. #LT E_{#nu} #GT");
    hists_area_scaled[0]->Draw("colz");
    c1->cd(3);
    hists_area_scaled[1]->SetTitle("Normalized Contour Areas for #epsilon vs. #LT E_{#nu} #GT");
    hists_area_scaled[1]->Draw("colz");
    c1->cd(4);
    hists_area_scaled[2]->SetTitle("Normalized Contour Areas for #epsilon vs. #alpha");
    hists_area_scaled[2]->Draw("colz");
	
	
    //set specific aesthetic things here
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1,0);
    gStyle->SetPalette(kBird);
    gStyle->SetTitleFontSize(0.07);
    gStyle->SetNumberContours(99);
    
	TCanvas *c3 = new TCanvas("c3","c3",1100,700);
	c3->Divide(2, 2);
	c3->cd(1);
    gPad->SetRightMargin(0.14);
	//hists_absdiff[0]->SetTitle("Contour Areas for #alpha vs. #LT E_{#nu} #GT");
	hists_absdiff[0]->Draw("colz");
	c3->cd(3);
    gPad->SetRightMargin(0.12);
	//hists_absdiff[1]->SetTitle("Contour Areas for #epsilon vs. #LT E_{#nu} #GT");
	hists_absdiff[1]->Draw("colz");
	c3->cd(4);
    gPad->SetRightMargin(0.12);
	//hists_absdiff[2]->SetTitle("Contour Areas for #epsilon vs. #alpha");
	hists_absdiff[2]->Draw("colz");
     
    
    
}//end of code
