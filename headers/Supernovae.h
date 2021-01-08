//
//  Supernovae.h
//  Erin Conley (erin.conley@duke.edu)
//

#ifndef Supernovae_h
#define Supernovae_h

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// findContourLine: Given TH2D histogram and percent contained, finds the line
//                  needed to draw as contour
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Double_t findContourLine(TH2D *hist, Double_t percent){
    Double_t line = 0.0;
    
    Int_t numberColumns = hist->GetNbinsX();
    Int_t numberRows = hist->GetNbinsY();
    Double_t totalBinContent = 0.0;
    std::vector<Double_t> binContents;
    for(int col = 1; col <= numberColumns; ++col){
        for(int row = 1; row <= numberRows; ++row){
            totalBinContent += hist->GetBinContent(col, row);
            binContents.emplace_back(hist->GetBinContent(col, row));
        }
        
    }
    
    std::sort(binContents.begin(), binContents.end());
    //loop from back to front
    Double_t count = 0.0;
    for (size_t i = 0; i < binContents.size(); ++i) {
        //now define the end
        size_t fromEnd = binContents.size()-1-i;
        
        count += binContents[fromEnd];
        if(count >= totalBinContent*percent){
            line = binContents[fromEnd];
            break;
        }
        
    }
    
    return line;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// setTH2Dcontour: Finds and sets contour on TH2D histogram for a given percent
//                 contained
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TH2D* setTH2Dcontour(TH2D *hist, Double_t percent){
    
    Double_t contours[1];
    contours[0] = findContourLine(hist, percent);
    hist->SetContour(1, contours);
    
    return hist;
    
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// findAreaOfContour: For given TH2D histogram, finds contour corresponding to
//                    percent contained and calculates area
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Double_t findAreaOfContour(TH2D *hist, Double_t percent){
	Double_t area = -99.0;
	TCanvas *c = new TCanvas("c","c",1000,700);
	Double_t contours[1];
	contours[0] = findContourLine(hist, percent);
	hist->SetContour(1, contours);
	hist->Draw("CONT LIST");
	c->Update();
	TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	TList* contLevel = NULL;
	TGraph* curv = NULL;
	TGraph* gc = NULL;
		
	Int_t nGraphs = 0;
	Int_t TotalConts = 0;
		
	if(conts == NULL){
		std::cout << "*** No Contours Were Extracted!" << std::endl;
		TotalConts = 0;
	}
	else TotalConts = conts->GetSize();
		
	//std::cout << "TotalConts = " << TotalConts << std::endl;
		
	for(int iCont = 0; iCont < TotalConts; ++iCont){
		contLevel = (TList*)conts->At(iCont);
		
		if(contLevel->GetSize() > 0){
			curv = (TGraph*)contLevel->First();
			area = curv->Integral();	
		}
	}
	delete c;
	
	return area;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// make_supernova_hist: Returns TH1D histogram with statistical fluctuations
//                      (sampled from spectrum)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TH1D* make_supernova_hist(TH1D *spectrum, TString histname) {
    
    // returns a histogram
    gRandom->SetSeed(0);
    
    TH1D* hist = new TH1D(histname,"Spectrum",200, 0.5, 100.);
    
    //define total number of events to make; increase to eliminate bias 
   // Double_t totEvents = spectrum->Integral()*2.0;
    Double_t totEvents = spectrum->Integral();
        
    for(Double_t i = 0.0; i < totEvents; ++i){
        Double_t ent = spectrum->GetRandom();
        hist->Fill(ent);
    }
    
    for(int i = 1; i <= hist->GetNbinsX(); ++i){
        hist->SetBinError(i, std::sqrt(hist->GetBinContent(i)) );
    }
    
    
    //5 = merges five bins in one; so output will have 40 bins
    //"" = histogram is modified
    //0 = means we are not defining new bin low-edges
    hist->Rebin(5,"",0);
    
    return hist;
    
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// FakeSupernovaMethod: Stores flux parameter measurements for supernova event
//                      rates with statistical fluctuations. Generates spectra
//                      numberOfSN times.  The supernova is defined by a test
//                      spectrum with SNOwGLoBES smearing ts, nue-Ar40 cross
//                      section ts_xscn, and post-smearing efficiency model
//                      ts_effic. The detector assumptions are defined by a grid
//                      described by its parameter bounds (grid), SNOwGLoBES
//                      smearing (grid_smear), nue-Ar40 cross section model
//                      (grid_xscn), and efficiency model (grid_xscn).
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void FakeSupernovaMethod(Int_t numberOfSN, Double_t distanceSN, TString grid, TString grid_smear, TString grid_xscn, TString grid_effic, TString ts, TString ts_xscn, TString ts_effic){
	
	std::cout << "Producing fake supernovae fits for " << numberOfSN << " " << distanceSN << "kpc supernovae and the following input parameters:" << std::endl;
    std::cout << "	Grid: " << grid << " with smearing " << grid_smear << ", xscn " << grid_xscn << ", and effic " << grid_effic << std::endl;  
    std::cout << "	Test spectra: " << ts << " with xscn " << ts_xscn << " and effic " << ts_effic << std::endl;
	
	///////////////////////////////////////////////////////////////////////
    ////////////////////////////DEFINE INPUTS//////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    TString test_spect_filename = "input/test_spectra/pinched_test_smeared_sum_" + ts + "_" + ts_xscn + "_" + ts_effic + ".dat";
    
    TString smeareddir = "input/" + grid + "/smear_" + grid_smear + "_" + grid_xscn + "_" + grid_effic;    
    TString pinchedinfo = "input/pinched_info/pinched_info_" + grid + ".dat";
    TString gridinfo = "input/grid_info/grid_info_" + grid + ".dat";
	
    Double_t massfact = 1.0;//2.4; //so that we see the full DUNE far detector response
		
	//define string for output file
    //convert to string for output file
    TString numStr; numStr.Form("%d", numberOfSN);
    //convert to string for output file
    TString distanceStr; distanceStr.Form("%.2lf", distanceSN);	
    TString outfile = "out/supernovae_" + grid + "_smear" + grid_smear + "_" + grid_xscn + "_" + grid_effic + "_spectra" + ts + "_" + ts_xscn + "_" + ts_effic + "_" + distanceStr + "kpc_" + numStr + "events.root";
    
    //styling options
    gRandom->SetSeed(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1,0);
    
    ///////////////////////////////////////////////////////////////////////
    /////////////////////////DEFINE GRID BOUNDS////////////////////////////
    ///////////////////////////////////////////////////////////////////////  
    
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
    
    //Range of physical parameters
    Double_t mingoodalpha = first_alpha;//1.0;
    Double_t maxgoodalpha = last_alpha;//7.0;
    Double_t mingoode0 = first_e0;//0.0;
    Double_t maxgoode0 = last_e0;//20.;
    
    std::cout << "The " << grid << " grid follows this definition:" << std::endl;
    std::cout << "    Alpha: [" << first_alpha << ", " << last_alpha << "] with " << step_alpha << " spacing" << std::endl;
    std::cout << "    E0: [" << first_e0 << ", " << last_e0 << "] with " << step_e0 << " spacing" << std::endl;
    std::cout << "    Luminosity: [" << first_lum << ", " << last_lum << "] with " << step_lum << " spacing" << std::endl;
    
    std::cout << "The parameter-fitting algorithm will consider the following physical range for alpha and E0:" << std::endl;
    std::cout << "    Alpha: [" << mingoodalpha << ", " << maxgoodalpha << "]" << std::endl;
    std::cout << "    E0: [" << mingoode0 << ", " << maxgoode0 << "]" << std::endl;
     
    ///////////////////////////////////////////////////////////////////////
    /////////////////////////DEFINE OUTPUT FILE////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //define variables to keep the event information
    Double_t tr_eventid;
    Double_t tr_alpha;
    Double_t tr_e0;
    Double_t tr_lum;
    Double_t tr_chi2;
    Double_t tr_dof;
    Double_t tr_prob;
    Double_t tr_int;
    Double_t tr_maxval;
    Double_t tr_peakenergy;
    
    TTree *tr = new TTree("data", "data");
    tr->Branch("eventid", &tr_eventid, "eventid/D");
    tr->Branch("alpha", &tr_alpha, "alpha/D");
    tr->Branch("e0", &tr_e0, "e0/D");
    tr->Branch("lum", &tr_lum, "lum/D");
    tr->Branch("chi2", &tr_chi2, "chi2/D");
    tr->Branch("dof", &tr_dof, "dof/D");
    tr->Branch("prob", &tr_prob, "prob/D");
    tr->Branch("int", &tr_int, "int/D");
    tr->Branch("maxval", &tr_maxval, "maxval/D");
    tr->Branch("peakenergy", &tr_peakenergy, "peakenergy/D");
    
    ///////////////////////////////////////////////////////////////////////
    ////////////////////////READ PINCHING PARAMETERS///////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    //open the pinched parameter file
    ifstream pin;
    pin.open(pinchedinfo);
    
    //define arrays, iterators to keep the parameters
    const Int_t maxflux = 350000;
    
    std::vector<Int_t> pnum; pnum.resize(maxflux);
    std::vector<Double_t> anue; anue.resize(maxflux);
    std::vector<Double_t> bnue; bnue.resize(maxflux);
    std::vector<Double_t> cnue; cnue.resize(maxflux);
    std::vector<Double_t> anuebar; anuebar.resize(maxflux);
    std::vector<Double_t> bnuebar; bnuebar.resize(maxflux);
    std::vector<Double_t> cnuebar; cnuebar.resize(maxflux);
    std::vector<Double_t> anux; anux.resize(maxflux);
    std::vector<Double_t> bnux; bnux.resize(maxflux);
    std::vector<Double_t> cnux; cnux.resize(maxflux);
    
    std::cout << "Reading flux parameters: " << std::endl;
    
    Int_t ip = 0;
    while(1){
        
        //Int_t p1;
        //Double_t a1, a2, a3, b1, b2, b3, c1, c2, c3;
        //pin >> p1 >> a1 >> a2 >> a3 >> b1 >> b2 >> b3 >> c1 >> c2 >> c3;
        
        pin >> pnum[ip] >> anue[ip] >> anuebar[ip] >> anux[ip] >> bnue[ip] >> bnuebar[ip] >> bnux[ip] >> cnue[ip] >> cnuebar[ip] >> cnux[ip];
        if(!pin.good()) break;
        
        //std::cout << pnum[ip]<<" "<<anue[ip]<<" "<<anuebar[ip]<<" "<<anux[ip]<<" "<<bnue[ip]<<" "<<bnuebar[ip]<<" "<<bnux[ip]<<" "<<cnue[ip]<<" "<<cnuebar[ip]<<" "<<cnux[ip]<<std::endl;
        
        
        ++ip;
        
    }
    
    pin.close();
    
    std::cout << "The " << grid << " grid contains " << ip << " fluxes." << std::endl;
    
    ///////////////////////////////////////////////////////////////////////
    ////////////////GENERATE FAKE SUPERNOVAE; PERFORM FIT//////////////////
    ///////////////////////////////////////////////////////////////////////
    //get test spectrum == specific observed spectrum
    //generated from real model or from pinching parameters
    TH1D *test_spectrum = fill_spect_hist(test_spect_filename, "test_spectrum", distanceSN, massfact);

    //std::cout << "Total events in test spectra: " << test_spectrum->Integral() << std::endl;
    
    //define some values up here
	const Int_t maxhist = maxflux;
	Int_t iFlux; //iterator over the grid
	Int_t numfluxes = ip; //total number of fluxes we care about
	Double_t chi2min;
	Int_t ifbest;
	Int_t dofbest;
        
	Double_t chi2;
	Int_t dof;
            
    for(int iEvent = 0; iEvent < numberOfSN; ++iEvent){
        tr_eventid = iEvent;
        
        TString event; event.Form("%d", iEvent);
        
        chi2min = 100000000000.0;
        ifbest = 0;
        
        chi2 = 0.;
        dof = 0;
        
        TString name = "test_spectrum_" + event;
        TH1D *test_spectrum2 = make_supernova_hist(test_spectrum, name);
        
        tr_int = test_spectrum2->Integral();
        tr_maxval = test_spectrum2->GetMaximum();
        tr_peakenergy = test_spectrum2->GetBinCenter(test_spectrum2->GetMaximumBin());
        
        for(iFlux = 0; iFlux < numfluxes; ++iFlux){
            
            //run a check to make sure E0 and alpha are reasonable
            if(anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0 ){
                //define strings for the histograms
                TString pargridfilename = smeareddir+"/pinched_"+TString::Format("%d",iFlux)+"_smeared_sum.dat";
                TString histname = "pinched_" + TString::Format("%d", iFlux) + "_" + event;
                
                //fill histogram
                TH1D* pargridhist = fill_spect_hist(pargridfilename, histname, distanceSN, massfact);
                
                //now determine the chi2 and dof for this histogram
                mychi2(test_spectrum2, pargridhist, &chi2, &dof);
                
                //Note: chi2 is not per dof
                //chi2 /= dof;
                
                if(chi2 < chi2min){
                    chi2min = chi2;
                    dofbest = dof;
                    ifbest = iFlux;
                }
                
                //std::cout << "Chi2 " << iFlux << " " << chi2 << " " << dof << std::endl;
                
                //pargridhist[iFlux]->Draw("SAME");
                delete pargridhist;
            }
            
        }//finish loop over fluxes
        
        std::cout << "Best-fit measurement: (" << anue[ifbest] << ", " << bnue[ifbest] << ", " << cnue[ifbest] << ") with reduced chi2 = " << chi2min/dofbest << std::endl;
        
        ///////////////////////////////////////////////////////////////////////
        //////////////////////////DETERMINE BEST FIT///////////////////////////
        ///////////////////////////////////////////////////////////////////////
        
        //Find required scale factor for the chi2 prob to be 10%
        //Double_t avg_chi2min_thresh = TMath::ChisquareQuantile(0.9, dof);
        //Double_t reqscale = (avg_chi2min_thresh - dof)/chi2min;
        
        tr_alpha = anue[ifbest];
        tr_e0 = bnue[ifbest];
        tr_lum = cnue[ifbest]/factor;
        tr_chi2 = chi2min;
        tr_dof = dof;
        tr_prob = TMath::Prob(chi2min,dof);
        tr->Fill();
        
        std::cout << "Event " << iEvent << " finished." << std::endl;
        
        
    }//finish loop over SN events
    
    TFile *fout = new TFile(outfile, "recreate");
    tr->Write();
    fout->Close();
    
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// DrawFakeSupernova: Draws supernova event rates with statistical fluctuations
//                    along with its flux parameter measurements. The supernova
//                    is distanceSN kpc from Earth and is defined by a test
//                    spectrum with SNOwGLoBES smearing ts, nue-Ar40 cross
//                    section ts_xscn, and post-smearing efficiency model
//                    ts_effic. The detector assumptions are defined by a grid
//                    described by its parameter bounds (grid), SNOwGLoBES
//                    smearing (grid_smear), nue-Ar40 cross section model
//                    (grid_xscn), and efficiency model (grid_xscn).
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DrawFakeSupernova(Double_t distanceSN, TString grid, TString grid_smear, TString grid_xscn, TString grid_effic, TString ts, TString ts_xscn, TString ts_effic){
	
	std::cout << "Drawing fake supernova spectrum + fit for a " << distanceSN << "kpc supernova and the following input parameters:" << std::endl;
    std::cout << "	Grid: " << grid << " with smearing " << grid_smear << ", xscn " << grid_xscn << ", and effic " << grid_effic << std::endl;  
    std::cout << "	Test spectra: " << ts << " with xscn " << ts_xscn << " and effic " << ts_effic << std::endl;
	
	///////////////////////////////////////////////////////////////////////
    ////////////////////////////DEFINE INPUTS//////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    TString test_spect_filename = "input/test_spectra/pinched_test_smeared_sum_" + ts + "_" + ts_xscn + "_" + ts_effic + ".dat";
    
    TString smeareddir = "input/" + grid + "/smear_" + grid_smear + "_" + grid_xscn + "_" + grid_effic;    
    TString pinchedinfo = "input/pinched_info/pinched_info_" + grid + ".dat";
    TString gridinfo = "input/grid_info/grid_info_" + grid + ".dat";
    
    TString super_title = "Example Best-Fit for Randomly Sampled Supernova Spectrum";
	
    Double_t massfact = 1.0;//2.4; //so that we see the full DUNE far detector response

    //styling options
    gRandom->SetSeed(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1,0);
    gStyle->SetTitleFontSize(0.05);
    
    ///////////////////////////////////////////////////////////////////////
    /////////////////////////DEFINE GRID BOUNDS////////////////////////////
    ///////////////////////////////////////////////////////////////////////  
    
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
    
    //Range of physical parameters
    Double_t mingoodalpha = first_alpha;//1.0;
    Double_t maxgoodalpha = last_alpha;//7.0;
    Double_t mingoode0 = first_e0;//0.0;
    Double_t maxgoode0 = last_e0;//20.;
    
    std::cout << "The " << grid << " grid follows this definition:" << std::endl;
    std::cout << "    Alpha: [" << first_alpha << ", " << last_alpha << "] with " << step_alpha << " spacing" << std::endl;
    std::cout << "    E0: [" << first_e0 << ", " << last_e0 << "] with " << step_e0 << " spacing" << std::endl;
    std::cout << "    Luminosity: [" << first_lum << ", " << last_lum << "] with " << step_lum << " spacing" << std::endl;
    
    std::cout << "The parameter-fitting algorithm will consider the following physical range for alpha and E0:" << std::endl;
    std::cout << "    Alpha: [" << mingoodalpha << ", " << maxgoodalpha << "]" << std::endl;
    std::cout << "    E0: [" << mingoode0 << ", " << maxgoode0 << "]" << std::endl;
    
    ///////////////////////////////////////////////////////////////////////
    ////////////////////////READ PINCHING PARAMETERS///////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    //open the pinched parameter file
    ifstream pin;
    pin.open(pinchedinfo);
    
    //define arrays, iterators to keep the parameters
    const Int_t maxflux = 350000;
    
    std::vector<Int_t> pnum; pnum.resize(maxflux);
    std::vector<Double_t> anue; anue.resize(maxflux);
    std::vector<Double_t> bnue; bnue.resize(maxflux);
    std::vector<Double_t> cnue; cnue.resize(maxflux);
    std::vector<Double_t> anuebar; anuebar.resize(maxflux);
    std::vector<Double_t> bnuebar; bnuebar.resize(maxflux);
    std::vector<Double_t> cnuebar; cnuebar.resize(maxflux);
    std::vector<Double_t> anux; anux.resize(maxflux);
    std::vector<Double_t> bnux; bnux.resize(maxflux);
    std::vector<Double_t> cnux; cnux.resize(maxflux);
    
    std::cout << "Reading flux parameters: " << std::endl;
    
    Int_t ip = 0;
    while(1){
        
        //Int_t p1;
        //Double_t a1, a2, a3, b1, b2, b3, c1, c2, c3;
        //pin >> p1 >> a1 >> a2 >> a3 >> b1 >> b2 >> b3 >> c1 >> c2 >> c3;
        
        pin >> pnum[ip] >> anue[ip] >> anuebar[ip] >> anux[ip] >> bnue[ip] >> bnuebar[ip] >> bnux[ip] >> cnue[ip] >> cnuebar[ip] >> cnux[ip];
        if(!pin.good()) break;
        
        //std::cout << pnum[ip]<<" "<<anue[ip]<<" "<<anuebar[ip]<<" "<<anux[ip]<<" "<<bnue[ip]<<" "<<bnuebar[ip]<<" "<<bnux[ip]<<" "<<cnue[ip]<<" "<<cnuebar[ip]<<" "<<cnux[ip]<<std::endl;
        
        
        ++ip;
        
    }
    
    pin.close();
    
    std::cout << "The " << grid << " grid contains " << ip << " fluxes." << std::endl;
    
    ///////////////////////////////////////////////////////////////////////
    ////////////////GENERATE FAKE SUPERNOVAE; PERFORM FIT//////////////////
    ///////////////////////////////////////////////////////////////////////
    //get test spectrum == specific observed spectrum
    //generated from real model or from pinching parameters
    TH1D *test_spectrum = fill_spect_hist(test_spect_filename, "test_spectrum", distanceSN, massfact);

    //std::cout << "Total events in test spectra: " << test_spectrum->Integral() << std::endl;
    
    //define some values up here
	const Int_t maxhist = maxflux;
	Int_t iFlux; //iterator over the grid
	Int_t numfluxes = ip; //total number of fluxes we care about
	Double_t chi2min;
	Int_t ifbest;
        
	Double_t chi2;
	Int_t dof;
	
	chi2min = 100000000000.0;
	ifbest = 0;
        
	chi2 = 0.;
	dof = 0;
        
	TString name = "test_spectrum_statistical";
	TH1D *test_spectrum2 = make_supernova_hist(test_spectrum, name);
	//test_spectrum2->Scale(test_spectrum->Integral()/test_spectrum2->Integral());
	//test_spectrum2->Smooth();
	
	std::vector<TH1D*> hists_grid; hists_grid.resize(maxhist);
        
	for(iFlux = 0; iFlux < numfluxes; ++iFlux){
            
		//run a check to make sure E0 and alpha are reasonable
		if(anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0 ){
			//define strings for the histograms
			TString pargridfilename = smeareddir+"/pinched_"+TString::Format("%d",iFlux)+"_smeared_sum.dat";
			TString histname = "pinched_" + TString::Format("%d", iFlux);
                
			//fill histogram
			//TH1D* pargridhist = fill_spect_hist(pargridfilename, histname, distanceSN, massfact);
			//now determine the chi2 and dof for this histogram
			//mychi2(test_spectrum2, pargridhist, &chi2, &dof);
			
			//fill histogram
			hists_grid[iFlux] = fill_spect_hist(pargridfilename, histname, distanceSN, massfact);
			//now determine the chi2 and dof for this histogram
			mychi2(test_spectrum2, hists_grid[iFlux], &chi2, &dof);
                                
			if(chi2 < chi2min){
				chi2min = chi2;
				ifbest = iFlux;
			}
                
			//std::cout << "Chi2 " << iFlux << " " << chi2 << " " << dof << std::endl;
                
			//pargridhist[iFlux]->Draw("SAME");
			//delete pargridhist;
		}
            
	}//finish loop over fluxes
        
    Double_t alphabest = anue[ifbest];
	Double_t e0best = bnue[ifbest];
	Double_t lumbest = cnue[ifbest]/factor;    
	
	test_spectrum2->SetTitle(" ");
	test_spectrum2->SetMarkerStyle(23); //8
	test_spectrum2->SetMarkerSize(1.5); 
	test_spectrum2->SetMarkerColor(2);
	test_spectrum2->SetXTitle("Observed energy (MeV)");
	test_spectrum2->SetYTitle("Number of events");
	test_spectrum2->GetXaxis()->SetTitleSize(0.04);
	test_spectrum2->GetXaxis()->SetTitleOffset(1.3);
	test_spectrum2->GetXaxis()->SetLabelSize(0.04);
	test_spectrum2->GetYaxis()->SetTitleSize(0.04);
	test_spectrum2->GetYaxis()->SetTitleOffset(1.3);
	test_spectrum2->GetYaxis()->SetLabelSize(0.04);
	test_spectrum2->GetXaxis()->SetRange(0,25);
	test_spectrum2->SetLineWidth(0.);//(2.);
	test_spectrum2->SetLineColor(2);
        
	TH1D* bestfit_spect = (TH1D*)hists_grid[ifbest]->Clone("bestfit_spect");
	bestfit_spect->SetMarkerStyle(7);
	bestfit_spect->SetMarkerColor(1);
	bestfit_spect->SetLineColor(1.);
	bestfit_spect->SetLineWidth(3);
	bestfit_spect->SetLineStyle(2);
	
	TLegend* leg = new TLegend(0.43,0.73,0.86,0.87);   
	TString legstring = "#alpha ="+TString::Format("%4.1f",alphabest)+", #LT E_{#nu} #GT="+TString::Format("%4.1f",e0best)+" MeV, #varepsilon = "+TString::Format("%4.1f",lumbest) + "e53 ergs";
        
	leg->AddEntry(bestfit_spect,legstring, "ple");
	leg->AddEntry(test_spectrum2,"Test spectrum", "pe");  
	leg->SetFillColor(10);
	leg->SetTextFont((Font_t) 62);
	leg->SetTextSize(0.03);
	leg->SetBorderSize(0);
	
	TCanvas *c = new TCanvas("c", "c", 900, 700);
	test_spectrum2->SetTitle(super_title);
    test_spectrum2->Draw("p");   
	//pargridhist[ifbest]->Draw("same");
	bestfit_spect->Draw("same");    
    leg->Draw("same");
	
	
	//deallocate vectors
	pnum = std::vector<Int_t>();
    anue = std::vector<Double_t>(); 
    bnue = std::vector<Double_t>(); 
    cnue = std::vector<Double_t>(); 
    anuebar = std::vector<Double_t>(); 
    bnuebar = std::vector<Double_t>();
    cnuebar = std::vector<Double_t>(); 
    anux = std::vector<Double_t>(); 
    bnux = std::vector<Double_t>();
    cnux = std::vector<Double_t>(); 
    hists_grid = std::vector<TH1D*>();
    
}

#endif /* Supernovae_h */
