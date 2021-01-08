//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ScaleTestSpectrum
// Erin Conley (erin.conley@duke.edu)
// Description: For given test spectrum parameter assumptions and scaling
//              factors, plots a spectrum and its scaled counterparts. 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include "TH1D.h"
#include "TMatrixDBase.h"
#include "TH2.h"
#include "../headers/ForwardFit.h"

void ScaleTestSpectrum(){
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // HERE ARE THE VARIABLES THAT CHANGE
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //define the test spectrum we would like to shift
    TString ts("MARLEY_DCR_nueOnly");
    //for the legend
    TString fancy_ts("MARLEY w/ drift correction, (p, n) xscn, step efficiency");
    
    //define the cross section model we want to consider
    TString xscn("PNXscn");
    //define the efficiency model we want to consider
    TString effic("StepEffic");
    
    //define the shifts we want + names to superimpose
    std::vector<Double_t> scales; std::vector<TString> names;
    //define the scaling that we want to apply to the grid
	scales.emplace_back(-15.0); names.emplace_back("-15%");
    scales.emplace_back(-10.0); names.emplace_back("-10%");
    scales.emplace_back(-5.0); names.emplace_back("-5%");
    scales.emplace_back(-1.0); names.emplace_back("-1%");
    scales.emplace_back(1.0); names.emplace_back("+1%");
    scales.emplace_back(5.0); names.emplace_back("+5%");
    scales.emplace_back(10.0); names.emplace_back("+10%");
    scales.emplace_back(15.0); names.emplace_back("+15%");
    
    //for the legend
    std::vector<TString> legendLabels;
    legendLabels.emplace_back("Unshifted"); Color_t colUnshifted = (Color_t)kBlack;
    legendLabels.emplace_back("Shifted #pm 15%"); Color_t col15 = (Color_t)kBlue;
    legendLabels.emplace_back("Shifted #pm 10%"); Color_t col10 = (Color_t)kRed;
    legendLabels.emplace_back("Shifted #pm 5%"); Color_t col05 = (Color_t)(kGreen+2);
    legendLabels.emplace_back("Shifted #pm 1%"); Color_t col01 = (Color_t)kViolet;
    
    //define title for superimposed plot
    TString supertitle = "Scaled Test Spectra: " + fancy_ts + ";Observed Energy (MeV);Number of Events";
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE INPUTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //mass factor for DUNE statistics
    Double_t massfact = 1.0;//2.4; //so that we see the full DUNE far detector response
    Double_t distanceSN = 10.0; //input to make the histogram (not the dat file)
    
    //I want to make TGraph objects in order to (hopefully) make the superimposed plot
    //aesthetically pleasing and more readable.
    TString test_spectrum_file = "input/test_spectra/pinched_test_smeared_sum_" + ts + "_" + xscn + "_" + effic + ".dat";
    TString tmp_name = "test_spectrum_" + ts + "_" + xscn + "_" + effic;
    ifstream in;
    in.open(test_spectrum_file);
    
    std::vector<double> energyVector; std::vector<double> nueVector;
    while (1) {
		Double_t en, nue;
        in >> en>>nue;
        
        // Energy in file is in eV
        en *= 1000.;  // in MeV for the plots
                        
        if (!in.good()) break;
		
		energyVector.emplace_back(en);
		nueVector.emplace_back(nue);
		        
    }
    TGraph *g = makeTGraphFromVectors(energyVector, nueVector, tmp_name);
    g->SetLineColor(colUnshifted); g->SetLineWidth(3);
    
    std::vector<TGraph*> AllTGraphs;
    AllTGraphs.emplace_back(g);
    //the energy axis stays the same for everyone
    for(size_t iScale = 0; iScale < scales.size(); ++iScale){
		std::vector<double> nueShifted;
		TString nueShiftedTitle = "test_spectrum_shifted_" + ts + "_" + xscn + "_" + effic;
		for(int i = 0; i < g->GetN(); ++i){
			Double_t en = g->GetX()[i];
					
			Double_t enShift = en*(1.0+ (scales[iScale]/100.0) ); 
			Double_t evShift = g->Eval(enShift);
			
			//watch out for negative event rates!
			if(evShift < 0.0) evShift = 0.0;
			//also watch out for when the scaling is zero 
			//don't change the event rate!
			if(scales[iScale] == 0.0) evShift = g->GetY()[i];
			
			// File has actually events per bin
			Double_t binfact=1;
			evShift /=  binfact;
			
			//now scale for the distance we want
			evShift *= std::pow(10.0/distanceSN, 2);
			//and scale for the mass we want
			evShift *= massfact;
			
			//keep shifted event value
			nueShifted.emplace_back(evShift);
		}
		TGraph *g_tmp = makeTGraphFromVectors(energyVector, nueShifted, nueShiftedTitle);
		
		if(std::abs(scales[iScale]) == 15.0) g_tmp->SetLineColor(col15);
		else if(std::abs(scales[iScale]) == 10.0) g_tmp->SetLineColor(col10);
		else if(std::abs(scales[iScale]) == 5.0) g_tmp->SetLineColor(col05);
		else if(std::abs(scales[iScale]) == 1.0) g_tmp->SetLineColor(col01);
		
		g_tmp->SetLineWidth(3);
		
		AllTGraphs.emplace_back(g_tmp);
	}
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // PLOT
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	TCanvas *c = new TCanvas("c", "c", 1000, 700);
	TLegend *leg = makeLegend(0.5, 0.6, 0.9, 0.9);
	leg->SetTextSize(0.03);
	for(size_t iLeg = 0; iLeg < legendLabels.size(); ++iLeg) leg->AddEntry(AllTGraphs[iLeg], legendLabels[iLeg], "l");
	
	//now make a multi graph object
	auto mg = new TMultiGraph();
	for(size_t iGraph = 0; iGraph < AllTGraphs.size(); ++iGraph) mg->Add(AllTGraphs[iGraph]);
	
	mg->Draw("AL");
	mg->SetTitle(supertitle);
	mg->GetXaxis()->SetRangeUser(0,60);
	mg->GetXaxis()->SetTitleSize(0.04);
	mg->GetXaxis()->SetLabelSize(0.04);
	mg->GetYaxis()->SetTitleSize(0.04);
	mg->GetYaxis()->SetTitleOffset(1.);
	mg->GetYaxis()->SetLabelSize(0.04);
	mg->SetMinimum(0.0);
	gPad->Modified();
	leg->Draw("SAME");
    
}//end ScaleTestSpectrum
