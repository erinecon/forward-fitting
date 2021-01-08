//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// AsimovSmearingStudies_DetermineParameterKnowledge
// Erin Conley (erin.conley@duke.edu)
// Description: For given test spectrum parameter assumptions and scaling
//              factors, creates FOM plots for different combinations of
//              SNOwGLoBES smearing (used to produce the grid and test specta).
//              Helps to determine the flux parameter knowledge necessary to
//              constrain measurement bias
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

void AsimovSmearingStudies_DetermineParameterKnowledge(){
    
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
    //example: energy scaling study
    std::vector<TString> dcs; std::vector<double> smearing;
    dcs.emplace_back("MARLEY_DCR_NO_nueOnly-15Percent"); smearing.emplace_back(0);
    dcs.emplace_back("MARLEY_DCR_NO_nueOnly-10Percent"); smearing.emplace_back(1);
    dcs.emplace_back("MARLEY_DCR_NO_nueOnly-5Percent"); smearing.emplace_back(2);
    dcs.emplace_back("MARLEY_DCR_NO_nueOnly-1Percent"); smearing.emplace_back(3);
    dcs.emplace_back("MARLEY_DCR_NO_nueOnly"); smearing.emplace_back(4);
    dcs.emplace_back("MARLEY_DCR_NO_nueOnly+1Percent"); smearing.emplace_back(5);
    dcs.emplace_back("MARLEY_DCR_NO_nueOnly+5Percent"); smearing.emplace_back(6);
    dcs.emplace_back("MARLEY_DCR_NO_nueOnly+10Percent"); smearing.emplace_back(7);
    dcs.emplace_back("MARLEY_DCR_NO_nueOnly+15Percent"); smearing.emplace_back(8);
    
    //define the test spectrum we want to use
    std::vector<TString> tss; std::vector<double> spectra;
    tss.emplace_back("MARLEY_DCR_NO_nueOnly-15Percent"); spectra.emplace_back(0);
    tss.emplace_back("MARLEY_DCR_NO_nueOnly-10Percent"); spectra.emplace_back(1);
    tss.emplace_back("MARLEY_DCR_NO_nueOnly-5Percent"); spectra.emplace_back(2);
    tss.emplace_back("MARLEY_DCR_NO_nueOnly-1Percent"); spectra.emplace_back(3);
    tss.emplace_back("MARLEY_DCR_NO_nueOnly"); spectra.emplace_back(4);
    tss.emplace_back("MARLEY_DCR_NO_nueOnly+1Percent"); spectra.emplace_back(5);
    tss.emplace_back("MARLEY_DCR_NO_nueOnly+5Percent"); spectra.emplace_back(6);
    tss.emplace_back("MARLEY_DCR_NO_nueOnly+10Percent"); spectra.emplace_back(7);
    tss.emplace_back("MARLEY_DCR_NO_nueOnly+15Percent"); spectra.emplace_back(8);
    
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
    
    //define extra directory where the files live, e.g., "/EnergyScalingStudy/"
    //ps: if there isn't an extra directory and they live in "out", set the parameter to "/"
    TString dir("/");
	
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE INPUTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //define file names   
    std::vector<TString> files;
    std::vector<TString> filenames;
    std::vector<std::vector<double>> pairs;
    for(size_t i = 0; i < dcs.size(); ++i){
		
		for(size_t j = 0; j < tss.size(); ++j){
			
			TString tmp = "out"+dir+"EnergyScalingStudy/chi2plots_" + grid + "_smear" + dcs[i] + "_" + xscn + "_" + effic + "_spectra" + tss[j] + "_" + xscn + "_" + effic + "_" + distanceSN + ".root";
            
			files.emplace_back(tmp);
			
			TString tmpname = grid + "_smear" + dcs[i] + "_" + xscn + "_" + effic + "_spectra" + tss[j] + "_" + xscn + "_" + effic;
			filenames.emplace_back(tmpname);
			
			//std::cout << tmpname << std::endl;
			//std::cout << tmp << " " << smearing[i] << " " << spectra[j] << std::endl;
			
			std::vector<double> tmppair;
			tmppair.emplace_back(smearing[i]);
			tmppair.emplace_back(spectra[j]);
			pairs.emplace_back(tmppair);
						
		}	
		
	}
	
	//define binning, axis information for 2d plots
	Double_t *binsSmear = new Double_t[smearing.size()];
	for(size_t i = 0; i < smearing.size(); ++i) binsSmear[i] = smearing[i];
	Int_t numsmearbins = smearing.size()-1;
	
	Double_t *binsTS = new Double_t[spectra.size()];
	for(size_t i = 0; i < spectra.size(); ++i) binsTS[i] = spectra[i];
	Int_t numtsbins = spectra.size()-1;
	
	//truth points
	Double_t alpha_true = 2.5;
	Double_t e0_true = 9.5; //mev
	Double_t lum_true = 0.5; //10^{53} ergs
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GRAB EVENT INFORMATION
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
	
	//so we can loop through pairs assuming that
	//pairs.size() == areas[i].size()
	
	std::vector<TH2D*> hists_absdiff;
    std::vector<TH2D*> hists_absdiff_cut;
    
	//so we can loop through pairs assuming that
	//pairs.size() == alpha.size() == number of files
	
	//loop over the parameters
	for(Int_t iPar = 0; iPar < parNames.size(); ++iPar){
		
		TString titlediff = "Fractional difference from truth for " + parNames[iPar] ;//+ " ;Grid Energy Spectra Scaling (Percent);Test Energy Spectra Scaling (Percent)";
		TString namediff = "fracdiff_" + parNames[iPar];
        
        TString titlediff_cut = "Fractional difference from truth for " + parNames[iPar] + " known to >~ 10%";
        TString namediff_cut = "fracdiff_cut_" + parNames[iPar];
		
		TH2D *hdiff = new TH2D(namediff, titlediff, 9, 0, 8, 9, 0, 8);
        TH2D *hdiff_cut = new TH2D(namediff_cut, titlediff_cut, 9, 0, 8, 9, 0, 8);
		
		Double_t smallestDiff = 0.;
		
		//now loop through the pairs
		for(Int_t i = 0; i < pairs.size(); ++i){
			//pairs[i][0] = smearing
			//pairs[i][1] = spectra
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
            
            hdiff_cut->Fill(labels[index1], labels[index2], diffToFill);
            
		}
		
				
		hdiff->SetMinimum(smallestDiff-1e-20);
		hdiff->GetXaxis()->SetLabelSize(0.08);
		hdiff->GetYaxis()->SetLabelSize(0.08);
        hdiff->GetZaxis()->SetLabelSize(0.055);
		
		hdiff->GetXaxis()->SetTitleSize(0.04);
		hdiff->GetYaxis()->SetTitleSize(0.04);
		
		hdiff->GetXaxis()->SetTickLength(0.);
		hdiff->GetYaxis()->SetTickLength(0.);
        
        hdiff_cut->SetMarkerStyle(29);
        hdiff_cut->SetMarkerSize(2.0);
        hdiff_cut->SetMarkerColor(kRed);
		
		//now we save the hist by doing
		hists_absdiff.emplace_back(hdiff);
        hists_absdiff_cut.emplace_back(hdiff_cut);
		
	}
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DRAW 2D PLOTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //set specific aesthetic things here
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1,0);
    gStyle->SetPalette(kBird);
    gStyle->SetTitleFontSize(0.07);
    gStyle->SetPaintTextFormat("4.3f");
    gStyle->SetNumberContours(99);
    
	TCanvas *c3 = new TCanvas("c3","c3",1100,700);
	c3->Divide(2, 2);
	c3->cd(1);
    gPad->SetRightMargin(0.11);
	//hists_absdiff[0]->SetTitle("Contour Areas for #alpha vs. #LT E_{#nu} #GT");
	hists_absdiff[0]->Draw("colz");
    hists_absdiff_cut[0]->Draw("SAME TEXT");
	c3->cd(3);
    gPad->SetRightMargin(0.11);
	//hists_absdiff[1]->SetTitle("Contour Areas for #epsilon vs. #LT E_{#nu} #GT");
	hists_absdiff[1]->Draw("colz");
    hists_absdiff_cut[1]->Draw("SAME TEXT");
	c3->cd(4);
    gPad->SetRightMargin(0.11);
	//hists_absdiff[2]->SetTitle("Contour Areas for #epsilon vs. #alpha");
	hists_absdiff[2]->Draw("colz");
    hists_absdiff_cut[2]->Draw("SAME TEXT");
    
    
}//end of code
