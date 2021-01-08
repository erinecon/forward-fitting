//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// AsimovResolutionVsDistanceStudy_GeneratePlots
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

void AsimovResolutionVsDistanceStudy_GeneratePlots(){
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // HERE ARE THE VARIABLES WE CHANGE
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    std::vector<TString> dists;  std::vector<double> distance;
    dists.emplace_back("1.00kpc"); distance.emplace_back(1.0);
    dists.emplace_back("2.00kpc"); distance.emplace_back(2.0);
    dists.emplace_back("3.00kpc"); distance.emplace_back(3.0);
    dists.emplace_back("4.00kpc"); distance.emplace_back(4.0);
    dists.emplace_back("5.00kpc"); distance.emplace_back(5.0);
    dists.emplace_back("6.00kpc"); distance.emplace_back(6.0);
    dists.emplace_back("7.00kpc"); distance.emplace_back(7.0);
    dists.emplace_back("8.00kpc"); distance.emplace_back(8.0);
    dists.emplace_back("9.00kpc"); distance.emplace_back(9.0);
    dists.emplace_back("10.00kpc"); distance.emplace_back(10.0);
    
    std::vector<TString> dcs; std::vector<TString> names; std::vector<double> resolution;
    dcs.emplace_back("0.00MARLEY_NO_nueOnly"); names.emplace_back("Unit Smearing Matrix"); resolution.emplace_back(0.0);
    dcs.emplace_back("0.05MARLEY_NO_nueOnly"); names.emplace_back("MARLEY Smearing, 5% Resolution"); resolution.emplace_back(5.0);
    dcs.emplace_back("0.10MARLEY_NO_nueOnly"); names.emplace_back("Gaussian Smearing, 10% Resolution"); resolution.emplace_back(10.0);
    dcs.emplace_back("0.15MARLEY_NO_nueOnly"); names.emplace_back("Gaussian Smearing, 15% Resolution"); resolution.emplace_back(15.0);
    dcs.emplace_back("0.20MARLEY_NO_nueOnly"); names.emplace_back("Gaussian Smearing, 20% Resolution"); resolution.emplace_back(20.0);
    dcs.emplace_back("0.25MARLEY_NO_nueOnly"); names.emplace_back("Gaussian Smearing, 25% Resolution"); resolution.emplace_back(25.0);
    dcs.emplace_back("0.30MARLEY_NO_nueOnly"); names.emplace_back("Gaussian Smearing, 30% Resolution"); resolution.emplace_back(30.0);
	
	std::vector<TString> histNames;
	histNames.emplace_back("alphavse0_90cut");
	histNames.emplace_back("lumvse0_90cut");
	histNames.emplace_back("lumvsalpha_90cut");
    
    std::vector<TString> fancyHistNames;
    fancyHistNames.emplace_back("#alpha vs. #LT E_{#nu} #GT");
    fancyHistNames.emplace_back("#varepsilon vs. #LT E_{#nu} #GT");
    fancyHistNames.emplace_back("#varepsilon vs. #alpha");
    
    std::vector<TString> parNames;
    parNames.emplace_back("#alpha");
    parNames.emplace_back("#LT E_{#nu} #GT");
    parNames.emplace_back("#varepsilon");
	
    //define the grid we want to use
    TString grid("2019October17");
    
    //define the cross section model we want to use
    TString xscn("PNXscn");
    
    //define the efficiency we want to use
    TString effic("StepEffic");
	
	//finally, define steps for distance and resolution (for binning)
	Double_t step_dist = 1.0;
	Double_t step_res = 5.0;
    
    //define extra directory where the files live, e.g., "/CrossSectionStudy/"
    //ps: if there isn't an extra directory and they live in "out", set the parameter to "/"
    TString dir("/ResolutionVsDistanceStudies/");
	
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE INPUTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //define file names   
    std::vector<TString> files;
    std::vector<TString> filenames;
    std::vector<Double_t> refareas_alphavse0;
    std::vector<Double_t> refareas_lumvse0;
    std::vector<Double_t> refareas_lumvsalpha;
    //also define constant reference areas to scale everyone by the same factor
    Double_t referenceArea_alphavse0;
    Double_t referenceArea_lumvse0;
    Double_t referenceArea_lumvsalpha;
    getReferenceContourAreas(10.0, grid, xscn, effic, xscn, effic, &referenceArea_alphavse0, &referenceArea_lumvse0, &referenceArea_lumvsalpha);
    std::vector<Double_t> ConstantReferenceAreas;
    ConstantReferenceAreas.emplace_back(referenceArea_alphavse0);
    ConstantReferenceAreas.emplace_back(referenceArea_lumvse0);
    ConstantReferenceAreas.emplace_back(referenceArea_lumvsalpha);
    
    std::vector<std::vector<double>> pairs;
    for(size_t i = 0; i < dcs.size(); ++i){
		
		for(size_t j = 0; j < dists.size(); ++j){
			
			TString tmp = "out"+dir+"chi2plots_" + grid + "_smear" + dcs[i] + "_" + xscn + "_" + effic + "_spectra" + dcs[i] + "_" + xscn + "_" + effic + "_" + dists[j] + ".root";
			files.emplace_back(tmp);
			
			TString tmpname = grid + "_" + dcs[i] + "_" + xscn + "_" + effic + "_" + dists[j];
			filenames.emplace_back(tmpname);
            
            Double_t ref_alphavse0; Double_t ref_lumvse0; Double_t ref_lumvsalpha;
            
            getReferenceContourAreas(distance[j], grid, xscn, effic, xscn, effic, &ref_alphavse0, &ref_lumvse0, &ref_lumvsalpha);
            
            //std::cout << ref_alphavse0 << " " << ref_lumvse0 << " " << ref_lumvsalpha << std::endl;
            
            refareas_alphavse0.emplace_back(ref_alphavse0);
            refareas_lumvse0.emplace_back(ref_lumvse0);
            refareas_lumvsalpha.emplace_back(ref_lumvsalpha);
			
			std::vector<double> tmppair;
			tmppair.emplace_back(resolution[i]);
			tmppair.emplace_back(distance[j]);
			pairs.emplace_back(tmppair);
						
		}	
		
	}
    
    std::vector<std::vector<Double_t>> ReferenceAreas;
    ReferenceAreas.emplace_back(refareas_alphavse0);
    ReferenceAreas.emplace_back(refareas_lumvse0);
    ReferenceAreas.emplace_back(refareas_lumvsalpha);
    
	//set specific aesthetic things here
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1,0);
	gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(99);
    gStyle->SetTitleFontSize(0.06);
	
	//define binning, axis information for 2d plots
	Double_t first_dist = *std::min_element(distance.begin(), distance.end());
	Double_t last_dist = *std::max_element(distance.begin(), distance.end());
	Int_t numdistbins = int( last_dist - first_dist )/step_dist+1;
	Double_t first_res = *std::min_element(resolution.begin(), resolution.end());
	Double_t last_res = *std::max_element(resolution.begin(), resolution.end());
	Int_t numresbins = int( last_res - first_res )/step_res+1;
	
	Double_t distmin = first_dist - step_dist/2.0;
	Double_t distmax = last_dist + step_dist/2.0;
	Double_t resmin = first_res - step_res/2.0;
	Double_t resmax = last_res + step_res/2.0;
    
    //std::cout << "dist " << numdistbins << " " << distmin << " " << distmax << std::endl;
    //std::cout << "res " << numresbins << " " << resmin << " " << resmax << std::endl;
	
	//define contours
	double contours[1];
	contours[0] = 0.001;
    
    //define binning, axis information for 2d plots
    Double_t *binsSmear = new Double_t[resolution.size()];
    for(size_t i = 0; i < resolution.size(); ++i) binsSmear[i] = resolution[i];
    Int_t numsmearbins = resolution.size()-1;
    
    Double_t *binsDistance = new Double_t[distance.size()];
    for(size_t i = 0; i < distance.size(); ++i) binsDistance[i] = distance[i];
    Int_t numdistancebins = distance.size()-1;
	
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
			
            //std::cout << pairs[iFile][1] << " " << pairs[iFile][0] << " ";
            
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
    std::vector<TH2D*> hists_area_scaled_constant;
	//we want to plot:
	//x axis: distance from SN
	//y axis: resolution (percent of smearing)
	//z axis: contour area
	
	//so we can loop through pairs assuming that
	//pairs.size() == areas[i].size()
	
	//loop over the contour plots
	for(Int_t iPlot = 0; iPlot < areas.size(); ++iPlot){
		
		//histNames[iPlot] will give the appropriate histName
		//use this to define name, title of 2D hist 
        TString title = "Contour Areas for " + fancyHistNames[iPlot] ;//+ " ;Grid Energy Spectra Scaling (Percent);Test Energy Spectra Scaling (Percent)";
        TString name = "areas_" + histNames[iPlot];
		
		TH2D *h = new TH2D(name, title, numdistbins, distmin, distmax, numresbins, resmin, resmax);
        //TH2D *h = new TH2D(name, title, numdistancebins, binsDistance, numsmearbins, binsSmear);
        
        TString title1 = "Normalized Contour Areas for " + fancyHistNames[iPlot] ;//+ " ;Grid Energy Spectra Scaling (Percent);Test Energy Spectra Scaling (Percent)";
        TString name1 = "areas_scaled_" + histNames[iPlot];
        
        TH2D *h1 = new TH2D(name1, title1, numdistbins, distmin, distmax, numresbins, resmin, resmax);
        //TH2D *h = new TH2D(name, title, numdistancebins, binsDistance, numsmearbins, binsSmear);
        
        TString title2 = "Normalized Contour Areas for " + fancyHistNames[iPlot];//+ " ;Grid Energy Spectra Scaling (Percent);Test Energy Spectra Scaling (Percent)";
        TString name2 = "areas_scaled_constant_" + histNames[iPlot];
        
        TH2D *h2 = new TH2D(name2, title2, numdistbins, distmin, distmax, numresbins, resmin, resmax);
        
        //now loop through the pairs
		for(Int_t i = 0; i < pairs.size(); ++i){
			//pairs[i][0] = resolution
			//pairs[i][1] = distance
			//now we fill the hist by doing
			h->Fill(pairs[i][1], pairs[i][0], areas[iPlot][i]);
            h1->Fill(pairs[i][1], pairs[i][0], areas[iPlot][i]/ReferenceAreas[iPlot][i]);
            h2->Fill(pairs[i][1], pairs[i][0], areas[iPlot][i]/ConstantReferenceAreas[iPlot]);
			
		}
		
		h->SetMinimum(0.0);
        h->GetXaxis()->SetLabelSize(0.06);
        h->GetYaxis()->SetLabelSize(0.06);
        h->GetZaxis()->SetLabelSize(0.06);
        
        //h1->SetMinimum(0.0);
        h1->GetXaxis()->SetLabelSize(0.08);
        h1->GetYaxis()->SetLabelSize(0.08);
        h1->GetYaxis()->SetLabelOffset(0.01);
        h1->GetZaxis()->SetLabelSize(0.07);
        
        //h2->SetMinimum(0.0);
        h2->GetXaxis()->SetLabelSize(0.08);
        h2->GetYaxis()->SetLabelSize(0.08);
        h2->GetYaxis()->SetLabelOffset(0.01);
        h2->GetZaxis()->SetLabelSize(0.07);
		
		//now we save the hist by doing
		hists_area.emplace_back(h);
        hists_area_scaled.emplace_back(h1);
        hists_area_scaled_constant.emplace_back(h2);
		
	}
	
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DRAW 2D PLOTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	TCanvas *c1 = new TCanvas("c1","c1",1100,700);
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetTitleSize(0.06);
	c1->Divide(2, 2);
	c1->cd(1);
	hists_area[0]->Draw("colz");
	c1->cd(3);
	hists_area[1]->Draw("colz");
	c1->cd(4);
	hists_area[2]->Draw("colz");
    
    TCanvas *c2 = new TCanvas("c2","c2",400,1000);
    c2->Divide(1, 3);
    c2->cd(1);
    gPad->SetRightMargin(0.14);
    hists_area_scaled[0]->SetTitle("Normalized Contour Areas for #alpha vs. #LT E_{#nu} #GT");
    hists_area_scaled[0]->Draw("colz");
    c2->cd(2);
    gPad->SetRightMargin(0.14);
    hists_area_scaled[1]->SetTitle("Normalized Contour Areas for #varepsilon vs. #LT E_{#nu} #GT");
    hists_area_scaled[1]->Draw("colz");
    c2->cd(3);
    gPad->SetRightMargin(0.14);
    hists_area_scaled[2]->SetTitle("Normalized Contour Areas for #varepsilon vs. #alpha");
    hists_area_scaled[2]->Draw("colz");
    
    
    
    TCanvas *c3 = new TCanvas("c3","c3",400,1000);
    c3->Divide(1, 3);
    c3->cd(1);
    gPad->SetRightMargin(0.14);
    hists_area_scaled[0]->SetTitle("Normalized Contour Areas for #alpha vs. #LT E_{#nu} #GT");
    hists_area_scaled_constant[0]->Draw("colz");
    c3->cd(2);
    gPad->SetRightMargin(0.14);
    hists_area_scaled[1]->SetTitle("Normalized Contour Areas for #varepsilon vs. #LT E_{#nu} #GT");
    hists_area_scaled_constant[1]->Draw("colz");
    c3->cd(3);
    gPad->SetRightMargin(0.14);
    hists_area_scaled[2]->SetTitle("Normalized Contour Areas for #varepsilon vs. #alpha");
    hists_area_scaled_constant[2]->Draw("colz");
    
	

}//end of code
