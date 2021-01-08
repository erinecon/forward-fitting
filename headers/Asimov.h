//
//  Asimov.h
//  Erin Conley (erin.conley@duke.edu)
//

#ifndef Asimov_h
#define Asimov_h

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// InterpolateChi2Hist: Outputs an TH2D objects with number of row/column bins
//                      equal to numBins. Interpolates the bin content using
//                      histToInterpolate, and stores the bin content if it
//                      passes a cut defined by chi2cut
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TH2D* InterpolateChi2Hist(TH2D* histToInterpolate, Int_t numBins, Double_t chi2cut){
    
    TString name_interpolated = (TString)( histToInterpolate->GetName() ) + "_inter";
    TString title_interpolated = (TString)( histToInterpolate->GetTitle() ) + " (interpolated)";
    Double_t xmin = histToInterpolate->GetXaxis()->GetXmin();
    Double_t xmax = histToInterpolate->GetXaxis()->GetXmax();
    Double_t ymin = histToInterpolate->GetYaxis()->GetXmin();
    Double_t ymax = histToInterpolate->GetYaxis()->GetXmax();
    
    //hardcode the number of bins because I doubt they will ever need to change
    Int_t binCol = numBins;
    Int_t binRow = numBins;
    
    TH2D *h_interpolated = new TH2D(name_interpolated, title_interpolated, binCol, xmin, xmax, binRow, ymin, ymax);
    
    //loop over the bins
    for(int col = 1; col <= binCol; ++col){
        for(int row = 1; row <= binRow; ++row){
            
            //get the actual x, y for this col, row bin position
            double x = h_interpolated->GetXaxis()->GetBinCenter(col);
            double y = h_interpolated->GetYaxis()->GetBinCenter(row);
            
            //find the interpolated chi2 value
            double interbin = histToInterpolate->Interpolate(x, y);
            //if the interpolated chi2 value passes the cut, fill the hist
            if(interbin < chi2cut) h_interpolated->SetBinContent(col, row, interbin);
            
        }//end loop over rows
    }//end loop over columns
    
    return h_interpolated;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// findAreaOfContour: Calculates the area contained by the countor defined by
//                    the given TH2D object
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Double_t findAreaOfContour(TH2D *hist){
	Double_t area = -99.0;
	TCanvas *c = new TCanvas("c","c",1000,700);
	
	//specific contour set here for sensitivity study
	double contours[1];
	contours[0] = 0.001;
	
	hist->SetContour(1, contours);
	hist->Draw("CONT LIST");
	c->Update();
	TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	TList* contLevel = NULL;
		
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
            area = 0;
            TIter next(contLevel);
            TObject* obj =0;
            while(( obj = next() )){
                
                area += ((TGraph*)obj)->Integral();
            }
            
		}
	}
	delete c;
	
	return area;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// AsimovMethod: Produces chi2-minimized flux parameter measurements and plots
//               for a supernova located distanceSN from Earth. The supernova is
//               defined by a test spectrum with SNOwGLoBES smearing ts,
//               nue-Ar40 cross section ts_xscn, and post-smearing efficiency
//               model ts_effic. The detector assumptions are defined by a grid
//               described by its parameter bounds (grid), SNOwGLoBES smearing
//               (grid_smear), nue-Ar40 cross section model (grid_xscn), and
//               efficiency model (grid_xscn).
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void AsimovMethod(Double_t distanceSN, TString grid, TString grid_smear, TString grid_xscn, TString grid_effic, TString ts, TString ts_xscn, TString ts_effic){
    
    std::cout << "Producing sensitivity plots for a " << distanceSN << "kpc supernova and the following input parameters:" << std::endl;
    std::cout << "    Grid: " << grid << " with smearing " << grid_smear << ", xscn " << grid_xscn << ", and effic " << grid_effic << std::endl;
    std::cout << "    Test spectrum: " << ts << " with xscn " << ts_xscn << " and effic " << ts_effic << std::endl;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE INPUTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TString test_spect_filename = "input/test_spectra/pinched_test_smeared_sum_" + ts + "_" + ts_xscn + "_" + ts_effic + ".dat";
        
    TString smeareddir = "input/" + grid + "/smear_" + grid_smear + "_" + grid_xscn + "_" + grid_effic;
    TString pinchedinfo = "input/pinched_info/pinched_info_" + grid + ".dat";
    TString gridinfo = "input/grid_info/grid_info_" + grid + ".dat";
    //define output file
    TString distanceStr; distanceStr.Form("%.2lf", distanceSN);
    TString outfile = "out/chi2plots_" + grid + "_smear" + grid_smear + "_" + grid_xscn + "_" + grid_effic + "_spectra" + ts + "_" + ts_xscn + "_" + ts_effic + "_" + distanceStr + "kpc.root";
    
    Double_t alpha_true = 2.5;
    Double_t e0_true = 9.5; //mev
    Double_t lum_true = 5e52; //ergs
    
    Double_t massfact = 1.0;//2.4; //see the full DUNE far detector response
    Double_t chi2_90contour = 4.61;  //pdg 2018
    
    //use same style settings that Kate used
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1,0);
    //set random seed
    gRandom->SetSeed(0);
    
    //define tree for output file
    Double_t tr_alpha;
    Double_t tr_e0;
    Double_t tr_lum;
    Double_t tr_chi2;
    Double_t tr_dof;
    Double_t tr_numGridElementsAlphavsE0;
    Double_t tr_numGridElementsLumvsE0;
    Double_t tr_numGridElementsLumvsAlpha;
    Double_t tr_numGoodGridElementsAlphavsE0;
    Double_t tr_numGoodGridElementsLumvsE0;
    Double_t tr_numGoodGridElementsLumvsAlpha;
    TTree *tr = new TTree("data", "data");
    tr->Branch("alpha", &tr_alpha, "alpha/D");
    tr->Branch("e0", &tr_e0, "e0/D");
    tr->Branch("lum", &tr_lum, "lum/D");
    tr->Branch("chi2", &tr_chi2, "chi2/D");
    tr->Branch("dof", &tr_dof, "dof/D");
    tr->Branch("numgrid_alphavse0", &tr_numGridElementsAlphavsE0, "numgrid_alphavse0/D");
    tr->Branch("numgrid_lumvse0", &tr_numGridElementsLumvsE0, "numgrid_lumvse0/D");
    tr->Branch("numgrid_lumvsalpha", &tr_numGridElementsLumvsAlpha, "numgrid_lumvsalpha/D");
    tr->Branch("numgood_alphavse0", &tr_numGoodGridElementsAlphavsE0, "numgood_alphavse0/D");
    tr->Branch("numgood_lumvse0", &tr_numGoodGridElementsLumvsE0, "numgood_lumvse0/D");
    tr->Branch("numgood_lumvsalpha", &tr_numGoodGridElementsLumvsAlpha, "numgood_lumvsalpha/D");
    
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
    
    Double_t minalpha2 = first_alpha - step_alpha*2.0;
    Double_t maxalpha2 = last_alpha + step_alpha*2.0;
    Double_t mine02 = first_e0 - step_e0*2.0;
    Double_t maxe02 = last_e0 + step_e0*2.0;
    Double_t minlum2 = first_lum - step_lum*2.0;
    Double_t maxlum2 = last_lum + step_lum*2.0;
    
    Int_t numalphabins = int((last_alpha - first_alpha)/step_alpha)+1;
    Int_t nume0bins = int((last_e0 - first_e0)/step_e0)+1;
    Int_t numlumbins = (last_lum - first_lum)/(step_lum) + 1;
    
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

    //std::cout << numalphabins << " " << nume0bins << " " << numlumbins << std::endl;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // READ PINCHING PARAMETERS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
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
    
    std::cout << "Reading flux parameters: ";
    
    Int_t ip = 0;
    while(1){
        
        pin >> pnum[ip] >> anue[ip] >> anuebar[ip] >> anux[ip] >> bnue[ip] >> bnuebar[ip] >> bnux[ip] >> cnue[ip] >> cnuebar[ip] >> cnux[ip];
        if(!pin.good()) break;
        
        //std::cout << pnum[ip]<<" "<<anue[ip]<<" "<<anuebar[ip]<<" "<<anux[ip]<<" "<<bnue[ip]<<" "<<bnuebar[ip]<<" "<<bnux[ip]<<" "<<cnue[ip]<<" "<<cnuebar[ip]<<" "<<cnux[ip]<<std::endl;
        
        ++ip;
        
    }
    
    pin.close();
    
    std::cout << "The " << grid << " grid contains " << ip << " fluxes." << std::endl;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAKE CHI2 MAP
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    TH1D* test_spectrum = fill_spect_hist(test_spect_filename,"test_spectrum",distanceSN, massfact);
    
    //now we loop over templates in the grid
    Int_t iFlux; //iterator over the grid
    Int_t numfluxes = ip; //total number of fluxes we care about
    
    const Int_t maxhist = maxflux;
    
    Double_t chi2min = 100000000000.0;
    Int_t dofbest;
    Int_t ifbest;
    
    Double_t chi2;
    Int_t dof;
    
    //keep values of chi2 for each flux
    std::vector<Double_t> chi2values; chi2values.resize(maxhist);
    
    //test
    std::vector<Double_t> chi2values_tmp; std::vector<Int_t> iflux_tmp;
    
    for(iFlux = 0; iFlux < numfluxes; ++iFlux){
        
        // Do not go outside the limits
        if(anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0 ){
            
            //define filename
            TString pargridfilename = smeareddir+"/pinched_"+TString::Format("%d",iFlux)+"_smeared_sum.dat";
            TString histname = "pinched_"+grid + "_" + grid_smear + "_" + grid_xscn + "_" + grid_effic + "_" + ts + "_" + ts_xscn + "_" + ts_effic + "_" + distanceStr + "kpc_"+TString::Format("%d",iFlux);
            
            TH1D *pargridhist = fill_spect_hist(pargridfilename, histname, distanceSN, massfact);
            mychi2(test_spectrum,pargridhist,&chi2,&dof);
            
            // chi2 here is not per dof
            //      chi2 /= dof;
            
            //save value
            chi2values[iFlux] = chi2/dof;
            
            //cout << "Chi2 "<< iFlux<<" "<< chi2<<" "<<dof<<endl;
            
            if (chi2<chi2min) {
                chi2min = chi2;
                dofbest = dof;
                ifbest = iFlux;
            }
            
            chi2values_tmp.emplace_back(chi2);
            iflux_tmp.emplace_back(iFlux);
            
        }
    }
    
    std::cout << "Best-fit measurement: (" << anue[ifbest] << ", " << bnue[ifbest] << ", " << cnue[ifbest] << ") with reduced chi2 = " << chi2min/dofbest << std::endl;
    
    int index = std::distance(chi2values_tmp.begin(), std::min_element(chi2values_tmp.begin(), chi2values_tmp.end()) );
    
    //now I want to make plots of chi2 vs parameter for the other two parameters
    //at their fixed truth values!
    std::vector<double> alphachi2;
    std::vector<double> alphavals;
    TString titlealpha("Minimum reduced #chi^{2} vs. #alpha;#alpha;Reduced #chi^{2}");
    TString namealpha("Chi2VsAlpha");
    std::vector<double> e0chi2;
    std::vector<double> e0vals;
    TString titlee0("Minimum reduced #chi^{2} vs. #LT E_{#nu} #GT;#LT E_{#nu} #GT (MeV);Reduced #chi^{2}");
    TString namee0("Chi2VsE0");
    std::vector<double> lumchi2;
    std::vector<double> lumvals;
    TString titlelum("Minimum reduced #chi^{2} vs. #varepsilon;#varepsilon (10^{53} erg);Reduced #chi^{2}");
    TString namelum("Chi2VsLum");
    
    //we also want to keep the alpha, e0, lum values corresponding to "good" chi2 values
    TH2D *hAllowedRegion_AlphaVsE0 = new TH2D("alphavse0", ";#LT E_{#nu} #GT (MeV);#alpha", nume0bins, mine0, maxe0, numalphabins, minalpha, maxalpha );
    TH2D *hAllowedRegion_LumVsE0 = new TH2D("lumvse0", ";#LT E_{#nu} #GT (MeV);#varepsilon (10^{53} erg)", nume0bins, mine0, maxe0, numlumbins, minlum/factor, maxlum/factor );
    TH2D *hAllowedRegion_LumVsAlpha = new TH2D("lumvsalpha", ";#alpha;#varepsilon (10^{53} erg)", numalphabins, minalpha, maxalpha, numlumbins, minlum/factor, maxlum/factor );
    
    TH2D *hAllowedRegion_AlphaVsE0_90Cut = new TH2D("alphavse0_90cut", ";#LT E_{#nu} #GT (MeV);#alpha", nume0bins, mine0, maxe0, numalphabins, minalpha, maxalpha );
    TH2D *hAllowedRegion_LumVsE0_90Cut = new TH2D("lumvse0_90cut", ";#LT E_{#nu} #GT (MeV);#varepsilon (10^{53} erg)", nume0bins, mine0, maxe0, numlumbins, minlum/factor, maxlum/factor );
    TH2D *hAllowedRegion_LumVsAlpha_90Cut = new TH2D("lumvsalpha_90cut", ";#alpha;#varepsilon (10^{53} erg)", numalphabins, minalpha, maxalpha, numlumbins, minlum/factor, maxlum/factor );
    
    //set to zero
    tr_numGoodGridElementsAlphavsE0 = 0;
    tr_numGoodGridElementsLumvsE0 = 0;
    tr_numGoodGridElementsLumvsAlpha = 0;
    
    for(Double_t iAlpha = first_alpha; iAlpha <= last_alpha+step_alpha; iAlpha += step_alpha){
        
        //iAlpha tells us the alpha value we care about
        std::vector<double> chi2vals;
        
        for(iFlux = 0; iFlux < numfluxes; ++iFlux){
            //check to make sure we're at the right alpha value, and that alpha/E0 are physical
            
            Int_t diff = ( anue[iFlux] - iAlpha )*100.0;
            
            //if(diff == 0.0){
            if(diff == 0.0 && anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0){
                
                chi2vals.emplace_back(chi2values[iFlux]);
            }
        }
        
        //std::cout << iAlpha << " " << chi2vals.size() << std::endl;
        
        //now we find the minimum chi2 from these values
        if(chi2vals.size() > 0){
            size_t iSmallest = std::distance(chi2vals.begin(), std::min_element( chi2vals.begin(), chi2vals.end() ));
            
            alphavals.emplace_back( iAlpha );
            
            alphachi2.emplace_back( chi2vals[iSmallest] );
            
            //std::cout << iAlpha << " " << chi2vals[iSmallest] << std::endl;
            
        }
        
        //2D plot
        for(Double_t iE0 = first_e0; iE0 <= last_e0; iE0 += step_e0){
            
            std::vector<double> chi2vals2d;
            
            for(iFlux = 0; iFlux < numfluxes; ++iFlux){
                //check to make sure we're at the right alpha value, and that alpha/E0 are physical
                
                Int_t diffalpha = ( anue[iFlux] - iAlpha )*100.0;
                Int_t diffe0 = ( bnue[iFlux] - iE0 )*100.0;
                
                //if(diffalpha == 0.0 && diffe0 == 0.0){
                if(diffalpha == 0.0 && diffe0 == 0.0 && anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0){
                    chi2vals2d.emplace_back(chi2values[iFlux]);
                }
                
            }
            
            //std::cout << iAlpha << " " << iE0 << " " << chi2vals2d.size() << std::endl;
            
            if(chi2vals2d.size() > 0){
                size_t iSmallest = std::distance(chi2vals2d.begin(), std::min_element( chi2vals2d.begin(), chi2vals2d.end() ));
                
                //std::cout << iAlpha << " " << iE0 << " " << chi2vals2d[iSmallest] << std::endl;
                Double_t toFill = chi2vals2d[iSmallest];
                if(toFill == 0.0) toFill = 1e-20;
                
                if(toFill <= chi2_90contour){
                    ++tr_numGoodGridElementsAlphavsE0;
                    hAllowedRegion_AlphaVsE0_90Cut->Fill(iE0, iAlpha, toFill);
                }
                
                hAllowedRegion_AlphaVsE0->Fill(iE0, iAlpha, toFill);
                
            }
            
            chi2vals2d = std::vector<double>();
            
        }
        
        //2D plot:
        for(Double_t iLum = first_lum; iLum <= last_lum+step_lum; iLum += step_lum){
            std::vector<double> chi2vals2d;
            
            for(iFlux = 0; iFlux < numfluxes; ++iFlux){
                //check to make sure we're at the right alpha value, and that alpha/E0 are physical
                Int_t diffalpha = ( anue[iFlux] - iAlpha )*100.0;
                Int_t difflum = ( cnue[iFlux]/factor - iLum/factor )*100.0;
                
                //if(diffalpha == 0.0 && difflum == 0.0){
                if(diffalpha == 0.0 && difflum == 0.0 && anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0){
                    chi2vals2d.emplace_back(chi2values[iFlux]);
                }
                
            }
            
            //std::cout << iAlpha << " " << iLum << " " << chi2vals2d.size() << std::endl;
            
            if(chi2vals2d.size() > 0){
                size_t iSmallest = std::distance(chi2vals2d.begin(), std::min_element( chi2vals2d.begin(), chi2vals2d.end() ));
                
                //std::cout << iAlpha << " " << iLum << " " << chi2vals2d[iSmallest] << std::endl;
                
                Double_t toFill = chi2vals2d[iSmallest];
                if(toFill == 0.0) toFill = 1e-20;
                
                if(toFill <= chi2_90contour){
                    ++tr_numGoodGridElementsLumvsAlpha;
                    hAllowedRegion_LumVsAlpha_90Cut->Fill(iAlpha, iLum/factor, toFill);
                }
                
                hAllowedRegion_LumVsAlpha->Fill(iAlpha, iLum/factor, toFill);
                
            }
            
            chi2vals2d = std::vector<double>();
        }
        
        chi2vals = std::vector<double>();
        
    }
    
    for(Double_t iE0 = first_e0; iE0 <= last_e0; iE0 += step_e0){
        
        //iE0 tells us the E0 value we care about
        std::vector<double> chi2vals;
        
        for(iFlux = 0; iFlux < numfluxes; ++iFlux){
            //check to make sure we're at the right alpha value, and that alpha/E0 are physical
            Int_t diffe0 = ( bnue[iFlux] - iE0 )*100.0;
            
            //if(diffe0 == 0.0){
            if(diffe0 == 0.0 && anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0){
                chi2vals.emplace_back(chi2values[iFlux]);
            }
        }
        
        //now we find the minimum chi2 from these values
        if(chi2vals.size() > 0){
            size_t iSmallest = std::distance(chi2vals.begin(), std::min_element( chi2vals.begin(), chi2vals.end() ));
            
            e0vals.emplace_back( iE0 );
            
            e0chi2.emplace_back( chi2vals[iSmallest] );
            
            //std::cout << iE0 << " " << chi2vals[iSmallest] << std::endl;
            
        }
        
        //2D plot
        for(Double_t iLum = first_lum; iLum <= last_lum+step_lum; iLum += step_lum){
            std::vector<double> chi2vals2d;
            
            for(iFlux = 0; iFlux < numfluxes; ++iFlux){
                //check to make sure we're at the right e0 value, and that alpha/E0 are physical
                Int_t diffe0 = ( bnue[iFlux] - iE0 )*100.0;
                Int_t difflum = ( cnue[iFlux]/factor - iLum/factor )*100.0;
                
                //if(diffe0 == 0.0 && difflum == 0.0){
                if(diffe0 == 0.0 && difflum == 0.0 && anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0){
                    chi2vals2d.emplace_back(chi2values[iFlux]);
                }
            }
            
            if(chi2vals2d.size() > 0){
                size_t iSmallest = std::distance(chi2vals2d.begin(), std::min_element( chi2vals2d.begin(), chi2vals2d.end() ));
                
                //std::cout << iE0 << " " << iLum << " " << chi2vals2d[iSmallest] << std::endl;
                
                Double_t toFill = chi2vals2d[iSmallest];
                if(toFill == 0.0) toFill = 1e-20;
                
                if(toFill <= chi2_90contour){
                    ++tr_numGoodGridElementsLumvsE0;
                    hAllowedRegion_LumVsE0_90Cut->Fill(iE0, iLum/factor, toFill);
                }
                
                hAllowedRegion_LumVsE0->Fill(iE0, iLum/factor, toFill);
                
            }
            
            chi2vals2d = std::vector<double>();
        }
        
        chi2vals = std::vector<double>();
        
    }
    
    for(Double_t iLum = first_lum; iLum <= last_lum+step_lum; iLum += step_lum){
        
        //iLum tells us the luminosity value we care about
        std::vector<double> chi2vals;
        
        for(iFlux = 0; iFlux < numfluxes; ++iFlux){
            
            //check to make sure we're at the right lum value, and that alpha/E0 are physical
            Int_t difflum = ( cnue[iFlux]/factor - iLum/factor )*100.0;
            //if(cnue[iFlux] == iLum){
            if(difflum == 0.0){
                
                if(anue[iFlux] < maxgoodalpha && bnue[iFlux] < maxgoode0 && anue[iFlux] > mingoodalpha && bnue[iFlux] > mingoode0) chi2vals.emplace_back(chi2values[iFlux]);
            }
        }
        
        //now we find the minimum chi2 from these values
        if(chi2vals.size() > 0){
            size_t iSmallest = std::distance(chi2vals.begin(), std::min_element( chi2vals.begin(), chi2vals.end() ));
            
            lumvals.emplace_back( iLum/factor );
            
            lumchi2.emplace_back( chi2vals[iSmallest] );
            
            //std::cout << iLum << " " << chi2vals[iSmallest] << std::endl;
            
        }
        
        chi2vals = std::vector<double>();
    }
    
    TGraph *g_Chi2VsAlpha = makeTGraphFromVectors(alphavals, alphachi2, titlealpha);
    TGraph *g_Chi2VsE0 = makeTGraphFromVectors(e0vals, e0chi2, titlee0);
    TGraph *g_Chi2VsLum = makeTGraphFromVectors(lumvals, lumchi2, titlelum);
    
    hAllowedRegion_AlphaVsE0->GetXaxis()->CenterTitle();
    hAllowedRegion_AlphaVsE0->GetXaxis()->SetTitleSize(0.05);
    hAllowedRegion_AlphaVsE0->GetYaxis()->CenterTitle();
    hAllowedRegion_AlphaVsE0->GetYaxis()->SetTitleSize(0.05);
    
    hAllowedRegion_LumVsE0->GetXaxis()->CenterTitle();
    hAllowedRegion_LumVsE0->GetXaxis()->SetTitleSize(0.04);
    hAllowedRegion_LumVsE0->GetYaxis()->CenterTitle();
    hAllowedRegion_LumVsE0->GetYaxis()->SetTitleSize(0.04);
    
    hAllowedRegion_LumVsAlpha->GetXaxis()->CenterTitle();
    hAllowedRegion_LumVsAlpha->GetXaxis()->SetTitleSize(0.05);
    hAllowedRegion_LumVsAlpha->GetYaxis()->CenterTitle();
    hAllowedRegion_LumVsAlpha->GetYaxis()->SetTitleSize(0.05);
    
    //change boundaries
    
    hAllowedRegion_AlphaVsE0->GetXaxis()->SetRangeUser(mine02, maxe02);
    hAllowedRegion_AlphaVsE0->GetYaxis()->SetRangeUser(minalpha2, maxalpha2);
    
    hAllowedRegion_LumVsE0->GetXaxis()->SetRangeUser(mine02, maxe02);
    hAllowedRegion_LumVsE0->GetYaxis()->SetRangeUser(minlum2/factor, maxlum2/factor);
    
    hAllowedRegion_LumVsAlpha->GetXaxis()->SetRangeUser(minalpha2, maxalpha2);
    hAllowedRegion_LumVsAlpha->GetYaxis()->SetRangeUser(minlum2/factor, maxlum2/factor);
    
    hAllowedRegion_AlphaVsE0_90Cut->GetXaxis()->SetRangeUser(mine02, maxe02);
    hAllowedRegion_AlphaVsE0_90Cut->GetYaxis()->SetRangeUser(minalpha2, maxalpha2);
    
    hAllowedRegion_LumVsE0_90Cut->GetXaxis()->SetRangeUser(mine02, maxe02);
    hAllowedRegion_LumVsE0_90Cut->GetYaxis()->SetRangeUser(minlum2/factor, maxlum2/factor);
    
    hAllowedRegion_LumVsAlpha_90Cut->GetXaxis()->SetRangeUser(minalpha2, maxalpha2);
    hAllowedRegion_LumVsAlpha_90Cut->GetYaxis()->SetRangeUser(minlum2/factor, maxlum2/factor);
    
    //limit 1D plots to ~5 sigma
    //if(g_Chi2VsAlpha->GetYaxis()->GetXmax() > 25.0) g_Chi2VsAlpha->SetMaximum(25.0);
    //if(g_Chi2VsE0->GetYaxis()->GetXmax() > 25.0) g_Chi2VsE0->SetMaximum(25.0);
    //if(g_Chi2VsLum->GetYaxis()->GetXmax() > 25.0) g_Chi2VsLum->SetMaximum(25.0);
    
    //set names
    g_Chi2VsAlpha->SetName(namealpha);
    g_Chi2VsE0->SetName(namee0);
    g_Chi2VsLum->SetName(namelum);
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // SAVE PLOTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //fill tree here
    tr_alpha = anue[ifbest];
    tr_e0 = bnue[ifbest];
    tr_lum = cnue[ifbest]/factor;
    tr_chi2 = chi2min;
    tr_dof = dofbest;
    tr_numGridElementsAlphavsE0 = numalphabins*nume0bins;
    tr_numGridElementsLumvsE0 = numlumbins*nume0bins;
    tr_numGridElementsLumvsAlpha = numalphabins*numlumbins;
    tr->Fill();
    
    TFile *fout = new TFile(outfile, "recreate");
    tr->Write();
    g_Chi2VsAlpha->Write();
    g_Chi2VsE0->Write();
    g_Chi2VsLum->Write();
    hAllowedRegion_AlphaVsE0->Write();
    hAllowedRegion_AlphaVsE0_90Cut->Write();
    hAllowedRegion_LumVsE0->Write();
    hAllowedRegion_LumVsE0_90Cut->Write();
    hAllowedRegion_LumVsAlpha->Write();
    hAllowedRegion_LumVsAlpha_90Cut->Write();
    fout->Close();
    
    //delete
    delete test_spectrum;
    delete hAllowedRegion_AlphaVsE0;
    delete hAllowedRegion_AlphaVsE0_90Cut;
    delete hAllowedRegion_LumVsE0;
    delete hAllowedRegion_LumVsE0_90Cut;
    delete hAllowedRegion_LumVsAlpha;
    delete hAllowedRegion_LumVsAlpha_90Cut;
    delete g_Chi2VsAlpha;
    delete g_Chi2VsE0;
    delete g_Chi2VsLum;
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
    chi2values = std::vector<Double_t>();
    alphavals = std::vector<double>();
    alphachi2 = std::vector<double>();
    e0vals = std::vector<double>();
    e0chi2 = std::vector<double>();
    lumvals = std::vector<double>();
    lumchi2 = std::vector<double>();
    
    //add line to separate the terminal output(s)
    std::cout << "--------------------------------" << std::endl;
    
}//end function

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// AsimovMethodStepEfficStudy: Produces chi2-minimized flux parameter
//                             measurements and plots for a supernova located
//                             distanceSN from Earth. The grid and test spectrum
//                             are both described using the same smearing
//                             (smear), nue-Ar40 cross section model (xscn), and
//                             post-smearing efficiency model (effic) which must
//                             be a "step efficiency" model, i.e., 100%
//                             detection efficiency above some threshold. The
//                             grid detection threshold (grid_shift) and test
//                             spectrum threshold (ts_shift) can be different
//                             values.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void AsimovMethodStepEfficStudy(Double_t distanceSN, TString grid, TString smear, TString xscn, TString effic, Double_t grid_shift, Double_t ts_shift){
	
    std::cout << "Step efficiency study for a " << distanceSN << "kpc supernova, grid " << grid << ", and the following:" << std::endl;
    std::cout << "	Smearing " << smear << ", xscn " << xscn << ", and effic " << effic << std::endl;  
    std::cout << "	Grid efficiency threshold " << grid_shift << " MeV; test spectrum efficiency threshold " << ts_shift << " MeV." << std::endl;
	
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE INPUTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TString test_spect_filename = "input/test_spectra/pinched_test_smeared_sum_" + smear + "_" + xscn + "_" + effic + ".dat";
        
    TString smeareddir = "input/" + grid + "/smear_" + smear + "_" + xscn + "_" + effic;    
    TString pinchedinfo = "input/pinched_info/pinched_info_" + grid + ".dat";
    TString gridinfo = "input/grid_info/grid_info_" + grid + ".dat";
    //define output file
    TString distanceStr; distanceStr.Form("%.2lf", distanceSN);
    TString gridShiftStr; gridShiftStr.Form("%.1lf", grid_shift);
    TString tsShiftStr; tsShiftStr.Form("%.1lf", ts_shift);
    TString outfile = "out/chi2plots_" + grid + "_smear" + smear + "_" + xscn + "_" + effic + gridShiftStr + "MeV_spectra" + smear + "_" + xscn + "_" + effic + tsShiftStr + "MeV_" + distanceStr + "kpc.root";    
    
    Double_t alpha_true = 2.5;
    Double_t e0_true = 9.5; //mev
    Double_t lum_true = 5e52; //ergs
    
    Double_t massfact = 1.0;//2.4; //see the full DUNE far detector response
    Double_t chi2_90contour = 4.61;  //pdg 2018
    
    //use same style settings that Kate used
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1,0);
    //set random seed
    gRandom->SetSeed(0);
    
    //define tree for output file
    Double_t tr_alpha;
    Double_t tr_e0;
    Double_t tr_lum;
    Double_t tr_chi2;
    Double_t tr_dof;
    TTree *tr = new TTree("data", "data");
    tr->Branch("alpha", &tr_alpha, "alpha/D");
    tr->Branch("e0", &tr_e0, "e0/D");
    tr->Branch("lum", &tr_lum, "lum/D");
    tr->Branch("chi2", &tr_chi2, "chi2/D");
    tr->Branch("dof", &tr_dof, "dof/D");
    
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
    
    Double_t minalpha2 = first_alpha - step_alpha*2.0;
    Double_t maxalpha2 = last_alpha + step_alpha*2.0;
    Double_t mine02 = first_e0 - step_e0*2.0;
    Double_t maxe02 = last_e0 + step_e0*2.0;
    Double_t minlum2 = first_lum - step_lum*2.0;
    Double_t maxlum2 = last_lum + step_lum*2.0;
    
    Int_t numalphabins = int(last_alpha - first_alpha)/step_alpha+1;
    Int_t nume0bins = int(last_e0 - first_e0)/step_e0+1;
    Int_t numlumbins = (last_lum - first_lum)/(step_lum) + 1;
    
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
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // READ PINCHING PARAMETERS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
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
    
    std::cout << "Reading flux parameters: ";
    
    Int_t ip = 0;
    while(1){
        
        pin >> pnum[ip] >> anue[ip] >> anuebar[ip] >> anux[ip] >> bnue[ip] >> bnuebar[ip] >> bnux[ip] >> cnue[ip] >> cnuebar[ip] >> cnux[ip];
        if(!pin.good()) break;
        
        //std::cout << pnum[ip]<<" "<<anue[ip]<<" "<<anuebar[ip]<<" "<<anux[ip]<<" "<<bnue[ip]<<" "<<bnuebar[ip]<<" "<<bnux[ip]<<" "<<cnue[ip]<<" "<<cnuebar[ip]<<" "<<cnux[ip]<<std::endl;
        
        ++ip;
        
    }
    
    pin.close();
    
    std::cout << "The " << grid << " grid contains " << ip << " fluxes." << std::endl;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAKE CHI2 MAP
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    TH1D* test_spectrum = shift_spect_hist(test_spect_filename,"test_spectrum",distanceSN, massfact, ts_shift); 
    
    //now we loop over templates in the grid
    Int_t iFlux; //iterator over the grid
    Int_t numfluxes = ip; //total number of fluxes we care about
    
    const Int_t maxhist = maxflux;
    
    Double_t chi2min = 100000000000.0;
    Int_t dofbest;
    Int_t ifbest;
    
    Double_t chi2;
    Int_t dof;
    
    //keep values of chi2 for each flux
    std::vector<Double_t> chi2values; chi2values.resize(maxhist);
    
    for(iFlux = 0; iFlux < numfluxes; ++iFlux){
        
        // Do not go outside the limits
        if(anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0 ){
            
            //define filename
            TString pargridfilename = smeareddir+"/pinched_"+TString::Format("%d",iFlux)+"_smeared_sum.dat";
            TString histname = "pinched_"+grid + "_" + smear + "_" + xscn + "_" + effic + gridShiftStr + "_" + smear + "_" + xscn + "_" + effic + tsShiftStr + "_" + distanceStr + "kpc_"+TString::Format("%d",iFlux);
            
            //TH1D *pargridhist = fill_spect_hist(pargridfilename, histname, distanceSN, massfact);
            TH1D *pargridhist = shift_spect_hist(pargridfilename, histname, distanceSN, massfact, grid_shift);
            
            mychi2(test_spectrum,pargridhist,&chi2,&dof);
            
            // chi2 here is not per dof
            //      chi2 /= dof;
            
            //save value
            chi2values[iFlux] = chi2/dof;
            
            if(chi2 == 0.0){ 
				//Double_t chi2test; Int_t doftest;				
				//switch the order of pargridhist, test_spectrum
				//if they are truly identical, then chi2 == 0
				mychi2(pargridhist, test_spectrum, &chi2, &dof);
				
			}
            
            //std::cout << anue[iFlux] << " " << bnue[iFlux] << " " << cnue[iFlux] << " " << chi2/dof << std::endl;
            
            //cout << "Chi2 "<< iFlux<<" "<< chi2<<" "<<dof<<endl;
            
            if (chi2<chi2min) {
                chi2min = chi2;
                dofbest = dof;
                ifbest = iFlux;
            }
            
        }
    }
    
    //fill tree here
    tr_alpha = anue[ifbest];
    tr_e0 = bnue[ifbest];
    tr_lum = cnue[ifbest]/factor;
    tr_chi2 = chi2min;
    tr_dof = dofbest;
    tr->Fill();
    
	std::cout << "Best-fit measurement: (" << anue[ifbest] << ", " << bnue[ifbest] << ", " << cnue[ifbest] << ") with reduced chi2 = " << chi2min/dofbest << std::endl;
    
    //now I want to make plots of chi2 vs parameter for the other two parameters
    //at their fixed truth values!
    std::vector<double> alphachi2;
    std::vector<double> alphavals;
    TString titlealpha("Minimum reduced #chi^{2} vs. #alpha;#alpha;Reduced #chi^{2}");
    TString namealpha("Chi2VsAlpha");
    std::vector<double> e0chi2;
    std::vector<double> e0vals;
    TString titlee0("Minimum reduced #chi^{2} vs. #LT E_{#nu} #GT;#LT E_{#nu} #GT (MeV);Reduced #chi^{2}");
    TString namee0("Chi2VsE0");
    std::vector<double> lumchi2;
    std::vector<double> lumvals;
    TString titlelum("Minimum reduced #chi^{2} vs. #epsilon;#epsilon (10^{53} erg);Reduced #chi^{2}");
    TString namelum("Chi2VsLum");
    
    //we also want to keep the alpha, e0, lum values corresponding to "good" chi2 values
    TH2D *hAllowedRegion_AlphaVsE0 = new TH2D("alphavse0", ";#LT E_{#nu} #GT (MeV);#alpha", nume0bins, mine0, maxe0, numalphabins, minalpha, maxalpha );
    TH2D *hAllowedRegion_LumVsE0 = new TH2D("lumvse0", ";#LT E_{#nu} #GT (MeV);#epsilon (10^{53} erg)", nume0bins, mine0, maxe0, numlumbins, minlum/factor, maxlum/factor );
    TH2D *hAllowedRegion_LumVsAlpha = new TH2D("lumvsalpha", ";#alpha;#epsilon (10^{53} erg)", numalphabins, minalpha, maxalpha, numlumbins, minlum/factor, maxlum/factor );
    
    TH2D *hAllowedRegion_AlphaVsE0_90Cut = new TH2D("alphavse0_90cut", ";#LT E_{#nu} #GT (MeV);#alpha", nume0bins, mine0, maxe0, numalphabins, minalpha, maxalpha );
    TH2D *hAllowedRegion_LumVsE0_90Cut = new TH2D("lumvse0_90cut", ";#LT E_{#nu} #GT (MeV);#epsilon (10^{53} erg)", nume0bins, mine0, maxe0, numlumbins, minlum/factor, maxlum/factor );
    TH2D *hAllowedRegion_LumVsAlpha_90Cut = new TH2D("lumvsalpha_90cut", ";#alpha;#epsilon (10^{53} erg)", numalphabins, minalpha, maxalpha, numlumbins, minlum/factor, maxlum/factor );
    
    for(Double_t iAlpha = first_alpha; iAlpha <= last_alpha+step_alpha; iAlpha += step_alpha){
        
        //iAlpha tells us the alpha value we care about
        std::vector<double> chi2vals;
        
        for(iFlux = 0; iFlux < numfluxes; ++iFlux){
            //check to make sure we're at the right alpha value, and that alpha/E0 are physical
            
            Int_t diff = ( anue[iFlux] - iAlpha )*100.0;
            
            if(diff == 0.0 && anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0){
                
                chi2vals.emplace_back(chi2values[iFlux]);
            }
        }
        
        //std::cout << iAlpha << " " << chi2vals.size() << std::endl;
        
        //now we find the minimum chi2 from these values
        if(chi2vals.size() > 0){
            size_t iSmallest = std::distance(chi2vals.begin(), std::min_element( chi2vals.begin(), chi2vals.end() ));
            
            alphavals.emplace_back( iAlpha );
            
            alphachi2.emplace_back( chi2vals[iSmallest] );
            
            //std::cout << iAlpha << " " << chi2vals[iSmallest] << std::endl;
            
        }
        
        //2D plot
        for(Double_t iE0 = first_e0; iE0 <= last_e0; iE0 += step_e0){
            
            std::vector<double> chi2vals2d;
            
            for(iFlux = 0; iFlux < numfluxes; ++iFlux){
                //check to make sure we're at the right alpha value, and that alpha/E0 are physical
                
                Int_t diffalpha = ( anue[iFlux] - iAlpha )*100.0;
                Int_t diffe0 = ( bnue[iFlux] - iE0 )*100.0;
                
                if(diffalpha == 0.0 && diffe0 == 0.0 && anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0){
                    chi2vals2d.emplace_back(chi2values[iFlux]);
                }
                
            }
            
            //std::cout << iAlpha << " " << iE0 << " " << chi2vals2d.size() << std::endl;
            
            if(chi2vals2d.size() > 0){
                size_t iSmallest = std::distance(chi2vals2d.begin(), std::min_element( chi2vals2d.begin(), chi2vals2d.end() ));
                
                //std::cout << iAlpha << " " << iE0 << " " << chi2vals2d[iSmallest] << std::endl;
                Double_t toFill = chi2vals2d[iSmallest];
                if(toFill == 0.0) toFill = 1e-20;
                
                if(toFill <= chi2_90contour) hAllowedRegion_AlphaVsE0_90Cut->Fill(iE0, iAlpha, toFill);
                
                hAllowedRegion_AlphaVsE0->Fill(iE0, iAlpha, toFill);
                
            }
            
            chi2vals2d = std::vector<double>();
            
        }
        
        //2D plot:
        for(Double_t iLum = first_lum; iLum <= last_lum+step_lum; iLum += step_lum){
            std::vector<double> chi2vals2d;
            
            for(iFlux = 0; iFlux < numfluxes; ++iFlux){
                //check to make sure we're at the right alpha value, and that alpha/E0 are physical
                Int_t diffalpha = ( anue[iFlux] - iAlpha )*100.0;
                Int_t difflum = ( cnue[iFlux]/factor - iLum/factor )*100.0;
                
                if(diffalpha == 0.0 && difflum == 0.0 && anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0){
                    chi2vals2d.emplace_back(chi2values[iFlux]);
                }
                
            }
            
            //std::cout << iAlpha << " " << iLum << " " << chi2vals2d.size() << std::endl;
            
            if(chi2vals2d.size() > 0){
                size_t iSmallest = std::distance(chi2vals2d.begin(), std::min_element( chi2vals2d.begin(), chi2vals2d.end() ));
                
                //std::cout << iAlpha << " " << iLum << " " << chi2vals2d[iSmallest] << std::endl;
                
                Double_t toFill = chi2vals2d[iSmallest];
                if(toFill == 0.0) toFill = 1e-20;
                
                if(toFill <= chi2_90contour) hAllowedRegion_LumVsAlpha_90Cut->Fill(iAlpha, iLum/factor, toFill);
                
                hAllowedRegion_LumVsAlpha->Fill(iAlpha, iLum/factor, toFill);
                
            }
            
            chi2vals2d = std::vector<double>();
        }
        
        chi2vals = std::vector<double>();
        
    }
    
    for(Double_t iE0 = first_e0; iE0 <= last_e0; iE0 += step_e0){
        
        //iE0 tells us the E0 value we care about
        std::vector<double> chi2vals;
        
        for(iFlux = 0; iFlux < numfluxes; ++iFlux){
            //check to make sure we're at the right alpha value, and that alpha/E0 are physical
            Int_t diffe0 = ( bnue[iFlux] - iE0 )*100.0;
            
            if(diffe0 == 0.0 && anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0){
                //if(bnue[iFlux] == iE0 && anue[iFlux] < maxgoodalpha && bnue[iFlux] < maxgoode0 && anue[iFlux] > mingoodalpha && bnue[iFlux] > mingoode0){
                chi2vals.emplace_back(chi2values[iFlux]);
            }
        }
        
        //now we find the minimum chi2 from these values
        if(chi2vals.size() > 0){
            size_t iSmallest = std::distance(chi2vals.begin(), std::min_element( chi2vals.begin(), chi2vals.end() ));
            
            e0vals.emplace_back( iE0 );
            
            e0chi2.emplace_back( chi2vals[iSmallest] );
            
            //std::cout << iE0 << " " << chi2vals[iSmallest] << std::endl;
            
        }
        
        //2D plot
        for(Double_t iLum = first_lum; iLum <= last_lum+step_lum; iLum += step_lum){
            std::vector<double> chi2vals2d;
            
            for(iFlux = 0; iFlux < numfluxes; ++iFlux){
                //check to make sure we're at the right e0 value, and that alpha/E0 are physical
                Int_t diffe0 = ( bnue[iFlux] - iE0 )*100.0;
                Int_t difflum = ( cnue[iFlux]/factor - iLum/factor )*100.0;
                if(diffe0 == 0.0 && difflum == 0.0 && anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0){
                    chi2vals2d.emplace_back(chi2values[iFlux]);
                }
            }
            
            if(chi2vals2d.size() > 0){
                size_t iSmallest = std::distance(chi2vals2d.begin(), std::min_element( chi2vals2d.begin(), chi2vals2d.end() ));
                
                //std::cout << iE0 << " " << iLum << " " << chi2vals2d[iSmallest] << std::endl;
                
                Double_t toFill = chi2vals2d[iSmallest];
                if(toFill == 0.0) toFill = 1e-20;
                
                if(toFill <= chi2_90contour) hAllowedRegion_LumVsE0_90Cut->Fill(iE0, iLum/factor, toFill);
                
                hAllowedRegion_LumVsE0->Fill(iE0, iLum/factor, toFill);
                
            }
            
            chi2vals2d = std::vector<double>();
        }
        
        chi2vals = std::vector<double>();
        
    }
    
    for(Double_t iLum = first_lum; iLum <= last_lum+step_lum; iLum += step_lum){
        
        //iLum tells us the luminosity value we care about
        std::vector<double> chi2vals;
        
        for(iFlux = 0; iFlux < numfluxes; ++iFlux){
            
            //check to make sure we're at the right lum value, and that alpha/E0 are physical
            Int_t difflum = ( cnue[iFlux]/factor - iLum/factor )*100.0;
            //if(cnue[iFlux] == iLum){
            if(difflum == 0.0){
                
                if(anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0) chi2vals.emplace_back(chi2values[iFlux]);
            }
        }
        
        //now we find the minimum chi2 from these values
        if(chi2vals.size() > 0){
            size_t iSmallest = std::distance(chi2vals.begin(), std::min_element( chi2vals.begin(), chi2vals.end() ));
            
            lumvals.emplace_back( iLum/factor );
            
            lumchi2.emplace_back( chi2vals[iSmallest] );
            
            //std::cout << iLum << " " << chi2vals[iSmallest] << std::endl;
            
        }
        
        chi2vals = std::vector<double>();
    }
    
    TGraph *g_Chi2VsAlpha = makeTGraphFromVectors(alphavals, alphachi2, titlealpha);
    TGraph *g_Chi2VsE0 = makeTGraphFromVectors(e0vals, e0chi2, titlee0);
    TGraph *g_Chi2VsLum = makeTGraphFromVectors(lumvals, lumchi2, titlelum);
    
    hAllowedRegion_AlphaVsE0->GetXaxis()->CenterTitle();
    hAllowedRegion_AlphaVsE0->GetXaxis()->SetTitleSize(0.05);
    hAllowedRegion_AlphaVsE0->GetYaxis()->CenterTitle();
    hAllowedRegion_AlphaVsE0->GetYaxis()->SetTitleSize(0.05);
    
    hAllowedRegion_LumVsE0->GetXaxis()->CenterTitle();
    hAllowedRegion_LumVsE0->GetXaxis()->SetTitleSize(0.04);
    hAllowedRegion_LumVsE0->GetYaxis()->CenterTitle();
    hAllowedRegion_LumVsE0->GetYaxis()->SetTitleSize(0.04);
    
    hAllowedRegion_LumVsAlpha->GetXaxis()->CenterTitle();
    hAllowedRegion_LumVsAlpha->GetXaxis()->SetTitleSize(0.05);
    hAllowedRegion_LumVsAlpha->GetYaxis()->CenterTitle();
    hAllowedRegion_LumVsAlpha->GetYaxis()->SetTitleSize(0.05);
    
    //change boundaries
    
    hAllowedRegion_AlphaVsE0->GetXaxis()->SetRangeUser(mine02, maxe02);
    hAllowedRegion_AlphaVsE0->GetYaxis()->SetRangeUser(minalpha2, maxalpha2);
    
    hAllowedRegion_LumVsE0->GetXaxis()->SetRangeUser(mine02, maxe02);
    hAllowedRegion_LumVsE0->GetYaxis()->SetRangeUser(minlum2/factor, maxlum2/factor);
    
    hAllowedRegion_LumVsAlpha->GetXaxis()->SetRangeUser(minalpha2, maxalpha2);
    hAllowedRegion_LumVsAlpha->GetYaxis()->SetRangeUser(minlum2/factor, maxlum2/factor);
    
    hAllowedRegion_AlphaVsE0_90Cut->GetXaxis()->SetRangeUser(mine02, maxe02);
    hAllowedRegion_AlphaVsE0_90Cut->GetYaxis()->SetRangeUser(minalpha2, maxalpha2);
    
    hAllowedRegion_LumVsE0_90Cut->GetXaxis()->SetRangeUser(mine02, maxe02);
    hAllowedRegion_LumVsE0_90Cut->GetYaxis()->SetRangeUser(minlum2/factor, maxlum2/factor);
    
    hAllowedRegion_LumVsAlpha_90Cut->GetXaxis()->SetRangeUser(minalpha2, maxalpha2);
    hAllowedRegion_LumVsAlpha_90Cut->GetYaxis()->SetRangeUser(minlum2/factor, maxlum2/factor);
    
    //limit 1D plots to ~5 sigma
    //if(g_Chi2VsAlpha->GetYaxis()->GetXmax() > 25.0) g_Chi2VsAlpha->SetMaximum(25.0);
    //if(g_Chi2VsE0->GetYaxis()->GetXmax() > 25.0) g_Chi2VsE0->SetMaximum(25.0);
    //if(g_Chi2VsLum->GetYaxis()->GetXmax() > 25.0) g_Chi2VsLum->SetMaximum(25.0);
    
    //set names
    g_Chi2VsAlpha->SetName(namealpha);
    g_Chi2VsE0->SetName(namee0);
    g_Chi2VsLum->SetName(namelum);
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // SAVE PLOTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    TFile *fout = new TFile(outfile, "recreate");
    tr->Write();
    g_Chi2VsAlpha->Write();
    g_Chi2VsE0->Write();
    g_Chi2VsLum->Write();
    hAllowedRegion_AlphaVsE0->Write();
    hAllowedRegion_AlphaVsE0_90Cut->Write();
    hAllowedRegion_LumVsE0->Write();
    hAllowedRegion_LumVsE0_90Cut->Write();
    hAllowedRegion_LumVsAlpha->Write();
    hAllowedRegion_LumVsAlpha_90Cut->Write();
    fout->Close();
    
    //delete
    delete test_spectrum;
    delete hAllowedRegion_AlphaVsE0;
    delete hAllowedRegion_AlphaVsE0_90Cut;
    delete hAllowedRegion_LumVsE0;
    delete hAllowedRegion_LumVsE0_90Cut;
    delete hAllowedRegion_LumVsAlpha;
    delete hAllowedRegion_LumVsAlpha_90Cut;
	delete g_Chi2VsAlpha;
    delete g_Chi2VsE0;
    delete g_Chi2VsLum;
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
    chi2values = std::vector<Double_t>();
    alphavals = std::vector<double>();
    alphachi2 = std::vector<double>();
	e0vals = std::vector<double>();
	e0chi2 = std::vector<double>(); 
	lumvals = std::vector<double>(); 
	lumchi2 = std::vector<double>();
	
	std::cout << "------------------------------------------------------" << std::endl;
	
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// AsimovMethodScaled: Produces chi2-minimized flux parameter measurements and
//                     plots for a supernova located distanceSN from Earth. The
//                     supernova is defined by a test spectrum with SNOwGLoBES
//                     smearing ts, nue-Ar40 cross section ts_xscn, and post-
//                     smearing efficiency model ts_effic. The energy scale of
//                     the test spectrum can be shifted by setting the ts_scale
//                     parameter in the range (-1, 1). The detector assumptions
//                     are defined by a grid described by its parameter bounds
//                     (grid), SNOwGLoBES smearing (grid_smear), nue-Ar40 cross
//                     section model (grid_xscn), and efficiency model
//                     (grid_xscn). The energy scale for each of the grid
//                     elements can be shifted by setting the grid_scale
//                     parameter in the range (-1, 1).
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void AsimovMethodScaled(Double_t distanceSN, TString grid, TString grid_smear, Double_t grid_scale, TString grid_xscn, TString grid_effic, TString ts, Double_t ts_scale, TString ts_xscn, TString ts_effic){
	
    std::cout << "Producing sensitivity plots for a " << distanceSN << "kpc supernova and the following input parameters:" << std::endl;
    std::cout << "	Grid: " << grid << " with smearing " << grid_smear << ", scaling " << grid_scale << "%, xscn " << grid_xscn << ", and effic " << grid_effic << std::endl;  
    std::cout << "	Test spectrum: " << ts << " with scaling " << ts_scale << "%, xscn " << ts_xscn << " and effic " << ts_effic << std::endl;
	
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE INPUTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TString test_spect_filename = "input/test_spectra/pinched_test_smeared_sum_" + ts + "_" + ts_xscn + "_" + ts_effic + ".dat";
    TString smeareddir = "input/" + grid + "/smear_" + grid_smear + "_" + grid_xscn + "_" + grid_effic;    
    TString pinchedinfo = "input/pinched_info/pinched_info_" + grid + ".dat";
    TString gridinfo = "input/grid_info/grid_info_" + grid + ".dat";
    //define output file
    TString distanceStr; distanceStr.Form("%.2lf", distanceSN);
    //need to define the percent;
    TString scaleGridStr; scaleGridStr.Form("%.0lf",grid_scale);
    TString scaleGridPercent;
    if(grid_scale > 0) scaleGridPercent = "+" + scaleGridStr + "Percent";
    else if(grid_scale < 0) scaleGridPercent = scaleGridStr + "Percent";
    TString scaleTSStr; scaleTSStr.Form("%.0lf",ts_scale);
    TString scaleTSPercent;
    if(ts_scale > 0) scaleTSPercent = "+" + scaleTSStr + "Percent";
    else if(ts_scale < 0) scaleTSPercent = scaleTSStr + "Percent";
    //then scalePercent should be nothing if scale_dc == 0.0
    TString outfile = "out/chi2plots_" + grid + "_smear" + grid_smear + scaleGridPercent + "_" + grid_xscn + "_" + grid_effic + "_spectra" + ts + scaleTSPercent + "_" + ts_xscn + "_" + ts_effic + "_" + distanceStr + "kpc.root";    
    
    std::cout << "Output file: " << outfile << std::endl;
    
    //define true parameters
	Double_t alpha_true = 2.5;
	Double_t e0_true = 9.5; //mev
	Double_t lum_true = 5e52; //ergs
	    
	Double_t massfact = 1.0;//2.4; //see the full DUNE far detector response
	Double_t chi2_90contour = 4.61;  //pdg 2018
    
    //use same style settings that Kate used
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1,0);
    //set random seed
    gRandom->SetSeed(0);
    
    //define tree for output file
    Double_t tr_alpha;
    Double_t tr_e0;
    Double_t tr_lum;
    Double_t tr_chi2;
    Double_t tr_dof;
    TTree *tr = new TTree("data", "data");
    tr->Branch("alpha", &tr_alpha, "alpha/D");
    tr->Branch("e0", &tr_e0, "e0/D");
    tr->Branch("lum", &tr_lum, "lum/D");
    tr->Branch("chi2", &tr_chi2, "chi2/D");
    tr->Branch("dof", &tr_dof, "dof/D");
       
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
    
    Double_t minalpha2 = first_alpha - step_alpha*2.0;
    Double_t maxalpha2 = last_alpha + step_alpha*2.0;
    Double_t mine02 = first_e0 - step_e0*2.0;
    Double_t maxe02 = last_e0 + step_e0*2.0;
    Double_t minlum2 = first_lum - step_lum*2.0;
    Double_t maxlum2 = last_lum + step_lum*2.0;
    
    Int_t numalphabins2 = int(maxalpha2 - minalpha2)/step_alpha+1;
    Int_t nume0bins2 = int(maxe02 - mine02)/step_e0+1;
    Int_t numlumbins2 = (maxlum2 - minlum2)/(step_lum) + 1;

	Int_t numalphabins = int(last_alpha - first_alpha)/step_alpha+1;
	Int_t nume0bins = int(last_e0 - first_e0)/step_e0+1;
	Int_t numlumbins = (last_lum - first_lum)/(step_lum) + 1;
    
    //Range of physical parameters
    Double_t mingoodalpha = first_alpha;//1.0;
    Double_t maxgoodalpha = last_alpha;//7.0;
    Double_t mingoode0 = first_e0;//0.0;
    Double_t maxgoode0 = last_e0;//20.;
    
    std::cout << "The " << grid << " grid follows this definition:" << std::endl;
    std::cout << "	Alpha: [" << first_alpha << ", " << last_alpha << "] with " << step_alpha << " spacing" << std::endl;
    std::cout << "	E0: [" << first_e0 << ", " << last_e0 << "] with " << step_e0 << " spacing" << std::endl;
    std::cout << "	Luminosity: [" << first_lum << ", " << last_lum << "] with " << step_lum << " spacing" << std::endl;
    
    std::cout << "The parameter-fitting algorithm will consider the following physical range for alpha and E0:" << std::endl;
    std::cout << "    Alpha: [" << mingoodalpha << ", " << maxgoodalpha << "]" << std::endl;
    std::cout << "    E0: [" << mingoode0 << ", " << maxgoode0 << "]" << std::endl;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // READ PINCHING PARAMETERS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
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
    
    std::cout << "Reading flux parameters: ";
    
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
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FIND BEST-FIT PARAMETERS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	//test spectrum already made with the appropriate scaling, so don't need to change this
	//TH1D* test_spectrum = fill_spect_hist(test_spect_filename,"test_spectrum",distanceSN, massfact); 
	TH1D* test_spectrum = scale_spect_hist(test_spect_filename, ts_scale, "test_spectrum", distanceSN, massfact);
	
	//now we loop over templates in the grid
	Int_t iFlux; //iterator over the grid
	Int_t numfluxes = ip; //total number of fluxes we care about
        
	const Int_t maxhist = maxflux;
        
	//define array to hold the hists
	//TH1D** pargridhist;
	//pargridhist = new TH1D*[maxhist];
        
	Double_t chi2min = 100000000000.0;
	Int_t dofbest;
	Int_t ifbest;
        
	Double_t chi2;
	Int_t dof;
	
	//keep values of chi2 for each flux
	//Double_t* chi2values = new Double_t[maxhist];
	//Double_t chi2values[maxhist];
	std::vector<Double_t> chi2values; chi2values.resize(maxhist);
	
	for(iFlux = 0; iFlux < numfluxes; ++iFlux){
		
		
		// Do not go outside the limits 
		if(anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0 ){
			
			//define filename 
			TString pargridfilename = smeareddir+"/pinched_"+TString::Format("%d",iFlux)+"_smeared_sum.dat";
			
			//TString histname = "pinched_"+TString::Format("%d",iFlux);
			TString histname = "pinched_"+grid + "_" + grid_smear + TString::Format("%.2lf",grid_scale)+ "_" + ts + "_"+TString::Format("%d",iFlux);
			
			
			TH1D* pargridhist = scale_spect_hist(pargridfilename, grid_scale, histname, distanceSN, massfact);
			
			//mychi2(test_spectrum,pargridhist[iFlux],&chi2,&dof);
			mychi2(test_spectrum,pargridhist,&chi2,&dof);
			// chi2 here is not per dof
			//      chi2 /= dof;
				
			//save value
			chi2values[iFlux] = chi2/dof;
				
			//std::cout << anue[iFlux] << " " << bnue[iFlux] << " " << cnue[iFlux] << " " << chi2/dof << std::endl;
			
			//cout << "Chi2 "<< iFlux<<" "<< chi2<<" "<<dof<<endl;
				
			if (chi2<chi2min) {
				chi2min = chi2;
				dofbest = dof;
				ifbest = iFlux;
			}
			
			//pargridhist[iFlux]->Draw("SAME");
			delete pargridhist;				
			
		}
	}
	
	//fill tree here
    tr_alpha = anue[ifbest];
    tr_e0 = bnue[ifbest];
    tr_lum = cnue[ifbest]/factor;
    tr_chi2 = chi2min;
    tr_dof = dofbest;
    tr->Fill();
	
	std::cout << "Best-fit measurement: (" << anue[ifbest] << ", " << bnue[ifbest] << ", " << cnue[ifbest] << ") with reduced chi2 = " << chi2min/dofbest << std::endl;
	
	//now I want to make plots of chi2 vs parameter for the other two parameters 
	//at their fixed truth values!
	std::vector<double> alphachi2;
	std::vector<double> alphavals;
	TString titlealpha("Minimum reduced #chi^{2} vs. #alpha;#alpha;Reduced #chi^{2}");
	TString namealpha("Chi2VsAlpha");
	std::vector<double> e0chi2;
	std::vector<double> e0vals;
	TString titlee0("Minimum reduced #chi^{2} vs. #LT E_{#nu} #GT;#LT E_{#nu} #GT (MeV);Reduced #chi^{2}");
	TString namee0("Chi2VsE0");
	std::vector<double> lumchi2;
	std::vector<double> lumvals;
	TString titlelum("Minimum reduced #chi^{2} vs. #epsilon;#epsilon (10^{53} erg);Reduced #chi^{2}");
	TString namelum("Chi2VsLum");
	
	//we also want to keep the alpha, e0, lum values corresponding to "good" chi2 values
	TH2D *hAllowedRegion_AlphaVsE0 = new TH2D("alphavse0", ";#LT E_{#nu} #GT (MeV);#alpha", nume0bins, mine0, maxe0, numalphabins, minalpha, maxalpha );
	TH2D *hAllowedRegion_LumVsE0 = new TH2D("lumvse0", ";#LT E_{#nu} #GT (MeV);#epsilon (10^{53} erg)", nume0bins, mine0, maxe0, numlumbins, minlum/factor, maxlum/factor );
	TH2D *hAllowedRegion_LumVsAlpha = new TH2D("lumvsalpha", ";#alpha;#epsilon (10^{53} erg)", numalphabins, minalpha, maxalpha, numlumbins, minlum/factor, maxlum/factor );
    
    TH2D *hAllowedRegion_AlphaVsE0_90Cut = new TH2D("alphavse0_90cut", ";#LT E_{#nu} #GT (MeV);#alpha", nume0bins, mine0, maxe0, numalphabins, minalpha, maxalpha );
    TH2D *hAllowedRegion_LumVsE0_90Cut = new TH2D("lumvse0_90cut", ";#LT E_{#nu} #GT (MeV);#epsilon (10^{53} erg)", nume0bins, mine0, maxe0, numlumbins, minlum/factor, maxlum/factor );
	TH2D *hAllowedRegion_LumVsAlpha_90Cut = new TH2D("lumvsalpha_90cut", ";#alpha;#epsilon (10^{53} erg)", numalphabins, minalpha, maxalpha, numlumbins, minlum/factor, maxlum/factor );
	
	for(Double_t iAlpha = first_alpha; iAlpha <= last_alpha+step_alpha; iAlpha += step_alpha){
		
		//iAlpha tells us the alpha value we care about
		std::vector<double> chi2vals;
		
		for(iFlux = 0; iFlux < numfluxes; ++iFlux){
			//check to make sure we're at the right alpha value, and that alpha/E0 are physical
			
			Int_t diff = ( anue[iFlux] - iAlpha )*100.0;
			
			if(diff == 0.0 && anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0){
				
				chi2vals.emplace_back(chi2values[iFlux]);
			}
		}
		
		//std::cout << iAlpha << " " << chi2vals.size() << std::endl;
		
		//now we find the minimum chi2 from these values 
		if(chi2vals.size() > 0){
			size_t iSmallest = std::distance(chi2vals.begin(), std::min_element( chi2vals.begin(), chi2vals.end() ));
			
			alphavals.emplace_back( iAlpha );
			
			alphachi2.emplace_back( chi2vals[iSmallest] );
			
			//std::cout << iAlpha << " " << chi2vals[iSmallest] << std::endl;
			
		}
		
		//2D plot
		for(Double_t iE0 = first_e0; iE0 <= last_e0; iE0 += step_e0){
			
			std::vector<double> chi2vals2d;
			
			for(iFlux = 0; iFlux < numfluxes; ++iFlux){
				//check to make sure we're at the right alpha value, and that alpha/E0 are physical
				
				Int_t diffalpha = ( anue[iFlux] - iAlpha )*100.0;
				Int_t diffe0 = ( bnue[iFlux] - iE0 )*100.0;
				
				if(diffalpha == 0.0 && diffe0 == 0.0 && anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0){
					chi2vals2d.emplace_back(chi2values[iFlux]);
				}
				
				
			}
			
			//std::cout << iAlpha << " " << iE0 << " " << chi2vals2d.size() << std::endl;
			
			if(chi2vals2d.size() > 0){
				size_t iSmallest = std::distance(chi2vals2d.begin(), std::min_element( chi2vals2d.begin(), chi2vals2d.end() ));
				
				//std::cout << iAlpha << " " << iE0 << " " << chi2vals2d[iSmallest] << std::endl;
				
				if(chi2vals2d[iSmallest] <= chi2_90contour) hAllowedRegion_AlphaVsE0_90Cut->Fill(iE0, iAlpha, chi2vals2d[iSmallest]);
				
				hAllowedRegion_AlphaVsE0->Fill(iE0, iAlpha, chi2vals2d[iSmallest]);
				
			}
			
		}
		
		//2D plot: 
		for(Double_t iLum = first_lum; iLum <= last_lum+step_lum; iLum += step_lum){
			std::vector<double> chi2vals2d;
			
			for(iFlux = 0; iFlux < numfluxes; ++iFlux){
				//check to make sure we're at the right alpha value, and that alpha/E0 are physical
				Int_t diffalpha = ( anue[iFlux] - iAlpha )*100.0;
				Int_t difflum = ( cnue[iFlux]/factor - iLum/factor )*100.0;
				
				if(diffalpha == 0.0 && difflum == 0.0 && anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0){
					chi2vals2d.emplace_back(chi2values[iFlux]);
				}
				
			}
			
			//std::cout << iAlpha << " " << iLum << " " << chi2vals2d.size() << std::endl;
			
			if(chi2vals2d.size() > 0){
				size_t iSmallest = std::distance(chi2vals2d.begin(), std::min_element( chi2vals2d.begin(), chi2vals2d.end() ));
				
				//std::cout << iAlpha << " " << iLum << " " << chi2vals2d[iSmallest] << std::endl;
								
				if(chi2vals2d[iSmallest] <= chi2_90contour) hAllowedRegion_LumVsAlpha_90Cut->Fill(iAlpha, iLum/factor, chi2vals2d[iSmallest]);
				
				hAllowedRegion_LumVsAlpha->Fill(iAlpha, iLum/factor, chi2vals2d[iSmallest]);
				
			}
		}
		
		//de-allocate vectors here
		chi2vals = std::vector<double>();
	}
	
	for(Double_t iE0 = first_e0; iE0 <= last_e0; iE0 += step_e0){
		
		//iE0 tells us the E0 value we care about
		std::vector<double> chi2vals;
		
		for(iFlux = 0; iFlux < numfluxes; ++iFlux){
			//check to make sure we're at the right alpha value, and that alpha/E0 are physical
			Int_t diffe0 = ( bnue[iFlux] - iE0 )*100.0;
			
			if(diffe0 == 0.0 && anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0){
				chi2vals.emplace_back(chi2values[iFlux]);
			}
		}
		
		//now we find the minimum chi2 from these values 
		if(chi2vals.size() > 0){
			size_t iSmallest = std::distance(chi2vals.begin(), std::min_element( chi2vals.begin(), chi2vals.end() ));
			
			e0vals.emplace_back( iE0 );
			
			e0chi2.emplace_back( chi2vals[iSmallest] );
			
			
			//std::cout << iE0 << " " << chi2vals[iSmallest] << std::endl;
			
		}
		
		//2D plot
		for(Double_t iLum = first_lum; iLum <= last_lum+step_lum; iLum += step_lum){
			std::vector<double> chi2vals2d;
			
			for(iFlux = 0; iFlux < numfluxes; ++iFlux){
				//check to make sure we're at the right e0 value, and that alpha/E0 are physical
				Int_t diffe0 = ( bnue[iFlux] - iE0 )*100.0;
				Int_t difflum = ( cnue[iFlux]/factor - iLum/factor )*100.0;
				if(diffe0 == 0.0 && difflum == 0.0 && anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0){
					chi2vals2d.emplace_back(chi2values[iFlux]);
				}
			}
			
			if(chi2vals2d.size() > 0){
				size_t iSmallest = std::distance(chi2vals2d.begin(), std::min_element( chi2vals2d.begin(), chi2vals2d.end() ));
				
				//std::cout << iE0 << " " << iLum << " " << chi2vals2d[iSmallest] << std::endl;
								
				if(chi2vals2d[iSmallest] <= chi2_90contour) hAllowedRegion_LumVsE0_90Cut->Fill(iE0, iLum/factor, chi2vals2d[iSmallest]);
				
				hAllowedRegion_LumVsE0->Fill(iE0, iLum/factor, chi2vals2d[iSmallest]);
				
			}
		}
		
		//de-allocate vectors here
		chi2vals = std::vector<double>();
		
	}
	
	for(Double_t iLum = first_lum; iLum <= last_lum+step_lum; iLum += step_lum){
		
		//iLum tells us the luminosity value we care about
		std::vector<double> chi2vals;
		
		for(iFlux = 0; iFlux < numfluxes; ++iFlux){
						
			//check to make sure we're at the right lum value, and that alpha/E0 are physical
			Int_t difflum = ( cnue[iFlux]/factor - iLum/factor )*100.0;
			//if(cnue[iFlux] == iLum){
			if(difflum == 0.0){
				
				if(anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0) chi2vals.emplace_back(chi2values[iFlux]);
			}
		}
		
		
		
		//now we find the minimum chi2 from these values 
		if(chi2vals.size() > 0){
			size_t iSmallest = std::distance(chi2vals.begin(), std::min_element( chi2vals.begin(), chi2vals.end() ));
			
			lumvals.emplace_back( iLum/factor );
			
			lumchi2.emplace_back( chi2vals[iSmallest] );
			
			//std::cout << iLum << " " << chi2vals[iSmallest] << std::endl;
			
		}
		
		//de-allocate vectors here
		chi2vals = std::vector<double>();
	}
		
	TGraph *g_Chi2VsAlpha = makeTGraphFromVectors(alphavals, alphachi2, titlealpha);
	TGraph *g_Chi2VsE0 = makeTGraphFromVectors(e0vals, e0chi2, titlee0);
	TGraph *g_Chi2VsLum = makeTGraphFromVectors(lumvals, lumchi2, titlelum);
	
	hAllowedRegion_AlphaVsE0->GetXaxis()->CenterTitle();
    hAllowedRegion_AlphaVsE0->GetXaxis()->SetTitleSize(0.05);
	hAllowedRegion_AlphaVsE0->GetYaxis()->CenterTitle();
    hAllowedRegion_AlphaVsE0->GetYaxis()->SetTitleSize(0.05);
    
    hAllowedRegion_LumVsE0->GetXaxis()->CenterTitle();
    hAllowedRegion_LumVsE0->GetXaxis()->SetTitleSize(0.04);
    hAllowedRegion_LumVsE0->GetYaxis()->CenterTitle();
    hAllowedRegion_LumVsE0->GetYaxis()->SetTitleSize(0.04);
    
    hAllowedRegion_LumVsAlpha->GetXaxis()->CenterTitle();
    hAllowedRegion_LumVsAlpha->GetXaxis()->SetTitleSize(0.05);
    hAllowedRegion_LumVsAlpha->GetYaxis()->CenterTitle();
    hAllowedRegion_LumVsAlpha->GetYaxis()->SetTitleSize(0.05);
    
    //change boundaries
    
    hAllowedRegion_AlphaVsE0->GetXaxis()->SetRangeUser(mine02, maxe02);
    hAllowedRegion_AlphaVsE0->GetYaxis()->SetRangeUser(minalpha2, maxalpha2);
    
    hAllowedRegion_LumVsE0->GetXaxis()->SetRangeUser(mine02, maxe02);
    hAllowedRegion_LumVsE0->GetYaxis()->SetRangeUser(minlum2/factor, maxlum2/factor);
    
    hAllowedRegion_LumVsAlpha->GetXaxis()->SetRangeUser(minalpha2, maxalpha2);
    hAllowedRegion_LumVsAlpha->GetYaxis()->SetRangeUser(minlum2/factor, maxlum2/factor);
    
    hAllowedRegion_AlphaVsE0_90Cut->GetXaxis()->SetRangeUser(mine02, maxe02);
    hAllowedRegion_AlphaVsE0_90Cut->GetYaxis()->SetRangeUser(minalpha2, maxalpha2);
    
    hAllowedRegion_LumVsE0_90Cut->GetXaxis()->SetRangeUser(mine02, maxe02);
    hAllowedRegion_LumVsE0_90Cut->GetYaxis()->SetRangeUser(minlum2/factor, maxlum2/factor);
    
    hAllowedRegion_LumVsAlpha_90Cut->GetXaxis()->SetRangeUser(minalpha2, maxalpha2);
    hAllowedRegion_LumVsAlpha_90Cut->GetYaxis()->SetRangeUser(minlum2/factor, maxlum2/factor);
	
	//set names
	g_Chi2VsAlpha->SetName(namealpha);
	g_Chi2VsE0->SetName(namee0);
	g_Chi2VsLum->SetName(namelum);
	
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // SAVE PLOTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	TFile *fout = new TFile(outfile, "recreate");
	tr->Write();
    g_Chi2VsAlpha->Write();
    g_Chi2VsE0->Write();
    g_Chi2VsLum->Write();
    hAllowedRegion_AlphaVsE0->Write();
    hAllowedRegion_AlphaVsE0_90Cut->Write();
    hAllowedRegion_LumVsE0->Write();
    hAllowedRegion_LumVsE0_90Cut->Write();
    hAllowedRegion_LumVsAlpha->Write();
    hAllowedRegion_LumVsAlpha_90Cut->Write();
    fout->Close();
    
    //delete
    delete test_spectrum;
    delete hAllowedRegion_AlphaVsE0;
    delete hAllowedRegion_AlphaVsE0_90Cut;
    delete hAllowedRegion_LumVsE0;
    delete hAllowedRegion_LumVsE0_90Cut;
    delete hAllowedRegion_LumVsAlpha;
    delete hAllowedRegion_LumVsAlpha_90Cut;
    delete g_Chi2VsAlpha;
    delete g_Chi2VsE0;
    delete g_Chi2VsLum;
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
    chi2values = std::vector<Double_t>();
    alphavals = std::vector<double>();
    alphachi2 = std::vector<double>();
	e0vals = std::vector<double>();
	e0chi2 = std::vector<double>(); 
	lumvals = std::vector<double>(); 
	lumchi2 = std::vector<double>();
    
    //add line to separate the terminal output(s)
    std::cout << "--------------------------------" << std::endl;
       
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// AsimovMethodSNDistance: Produces chi2-minimized flux parameter measurements
//                         and plots for a supernova located distanceSN +/-
//                         distanceUncertainty from Earth. The supernova is
//                         defined by a test spectrum with SNOwGLoBES smearing
//                         ts, nue-Ar40 cross section ts_xscn, and post-smearing
//                         efficiency model ts_effic. The detector assumptions
//                         are defined by a grid described by its parameter
//                         bounds (grid), SNOwGLoBES smearing (grid_smear),
//                         nue-Ar40 cross section model (grid_xscn), and
//                         efficiency model (grid_effic).
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void AsimovMethodSNDistance(Double_t distanceSN, Double_t distanceUncertainty, TString grid, TString grid_smear, TString grid_xscn, TString grid_effic, TString ts, TString ts_xscn, TString ts_effic){
    
    std::cout << "Producing sensitivity plots for a " << distanceSN << " +/- " << distanceSN*distanceUncertainty << " kpc supernova and the following input parameters:" << std::endl;
    std::cout << "	Grid: " << grid << " with smearing " << grid_smear << ", xscn " << grid_xscn << ", and effic " << grid_effic << std::endl;  
    std::cout << "	Test spectrum: " << ts << " with xscn " << ts_xscn << " and effic " << ts_effic << std::endl;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE INPUTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TString test_spect_filename = "input/test_spectra/pinched_test_smeared_sum_" + ts + "_" + ts_xscn + "_" + ts_effic + ".dat";
        
    TString smeareddir = "input/" + grid + "/smear_" + grid_smear + "_" + grid_xscn + "_" + grid_effic;    
    TString pinchedinfo = "input/pinched_info/pinched_info_" + grid + ".dat";
    TString gridinfo = "input/grid_info/grid_info_" + grid + ".dat";
    //define output file
    TString distanceStr; distanceStr.Form("%.2lf", distanceSN);
    TString distanceUncertaintyStr; distanceUncertaintyStr.Form("%.2lf", distanceUncertainty);
    TString outfile = "out/chi2plots_" + grid + "_smear" + grid_smear + "_" + grid_xscn + "_" + grid_effic + "_spectra" + ts + "_" + ts_xscn + "_" + ts_effic + "_" + distanceStr + "kpc_" + distanceUncertaintyStr + "uncertain.root";    
    
    Double_t alpha_true = 2.5;
    Double_t e0_true = 9.5; //mev
    Double_t lum_true = 5e52; //ergs
    
    Double_t massfact = 1.0;//2.4; //see the full DUNE far detector response
    Double_t chi2_90contour = 4.61;  //pdg 2018
    
    //use same style settings that Kate used
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1,0);
    //set random seed
    gRandom->SetSeed(0);
    
    //define tree for output file
    Double_t tr_alpha;
    Double_t tr_e0;
    Double_t tr_lum;
    Double_t tr_chi2;
    Double_t tr_dof;
    Double_t tr_numGridElementsAlphavsE0;
    Double_t tr_numGridElementsLumvsE0;
    Double_t tr_numGridElementsLumvsAlpha;
    Double_t tr_numGoodGridElementsAlphavsE0;
    Double_t tr_numGoodGridElementsLumvsE0;
    Double_t tr_numGoodGridElementsLumvsAlpha;
    TTree *tr = new TTree("data", "data");
    tr->Branch("alpha", &tr_alpha, "alpha/D");
    tr->Branch("e0", &tr_e0, "e0/D");
    tr->Branch("lum", &tr_lum, "lum/D");
    tr->Branch("chi2", &tr_chi2, "chi2/D");
    tr->Branch("dof", &tr_dof, "dof/D");
    tr->Branch("numgrid_alphavse0", &tr_numGridElementsAlphavsE0, "numgrid_alphavse0/D");
    tr->Branch("numgrid_lumvse0", &tr_numGridElementsLumvsE0, "numgrid_lumvse0/D");
    tr->Branch("numgrid_lumvsalpha", &tr_numGridElementsLumvsAlpha, "numgrid_lumvsalpha/D");
    tr->Branch("numgood_alphavse0", &tr_numGoodGridElementsAlphavsE0, "numgood_alphavse0/D");
    tr->Branch("numgood_lumvse0", &tr_numGoodGridElementsLumvsE0, "numgood_lumvse0/D");
    tr->Branch("numgood_lumvsalpha", &tr_numGoodGridElementsLumvsAlpha, "numgood_lumvsalpha/D");
    
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
    
    Double_t minalpha2 = first_alpha - step_alpha*2.0;
    Double_t maxalpha2 = last_alpha + step_alpha*2.0;
    Double_t mine02 = first_e0 - step_e0*2.0;
    Double_t maxe02 = last_e0 + step_e0*2.0;
    Double_t minlum2 = first_lum - step_lum*2.0;
    Double_t maxlum2 = last_lum + step_lum*2.0;
    
    Int_t numalphabins = int(last_alpha - first_alpha)/step_alpha+1;
    Int_t nume0bins = int(last_e0 - first_e0)/step_e0+1;
    Int_t numlumbins = (last_lum - first_lum)/(step_lum) + 1;
    
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

    //std::cout << numalphabins << " " << nume0bins << " " << numlumbins << std::endl;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // READ PINCHING PARAMETERS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
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
    
    std::cout << "Reading flux parameters: ";
    
    Int_t ip = 0;
    while(1){
        
        pin >> pnum[ip] >> anue[ip] >> anuebar[ip] >> anux[ip] >> bnue[ip] >> bnuebar[ip] >> bnux[ip] >> cnue[ip] >> cnuebar[ip] >> cnux[ip];
        if(!pin.good()) break;
        
        ++ip;
        
    }
    
    pin.close();
    
    std::cout << "The " << grid << " grid contains " << ip << " fluxes." << std::endl;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAKE CHI2 MAP
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    TH1D* test_spectrum = fill_spect_hist(test_spect_filename,"test_spectrum",distanceSN, massfact);
    
    //now we loop over templates in the grid
    Int_t iFlux; //iterator over the grid
    Int_t numfluxes = ip; //total number of fluxes we care about
    
    const Int_t maxhist = maxflux;
    
    Double_t chi2min = 100000000000.0;
    Int_t dofbest;
    Int_t ifbest;
    
    Double_t chi2;
    Int_t dof;
    
    //keep values of chi2 for each flux
    std::vector<Double_t> chi2values; chi2values.resize(maxhist);
    
    //test
    std::vector<Double_t> chi2values_tmp; std::vector<Int_t> iflux_tmp;
    
    for(iFlux = 0; iFlux < numfluxes; ++iFlux){
        
        // Do not go outside the limits
        if(anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0 ){
            
            //define filename
            TString pargridfilename = smeareddir+"/pinched_"+TString::Format("%d",iFlux)+"_smeared_sum.dat";
            TString histname = "pinched_"+grid + "_" + grid_smear + "_" + grid_xscn + "_" + grid_effic + "_" + ts + "_" + ts_xscn + "_" + ts_effic + "_" + distanceStr + "kpc_" + distanceUncertaintyStr + "uncertain_"+TString::Format("%d",iFlux);
            
            TH1D *pargridhist = fill_spect_hist(pargridfilename, histname, distanceSN, massfact);
            mychi2_distuncertainty(test_spectrum,pargridhist,distanceSN*distanceUncertainty,&chi2,&dof);
            
            // chi2 here is not per dof
            //      chi2 /= dof;
            
            //save value
            chi2values[iFlux] = chi2/dof;
			
            
            //cout << "Chi2 "<< iFlux<<" "<< chi2<<" "<<dof<<endl;
            
            if (chi2<chi2min) {
                chi2min = chi2;
                dofbest = dof;
                ifbest = iFlux;
            }
            
            chi2values_tmp.emplace_back(chi2);
            iflux_tmp.emplace_back(iFlux);
            
	    }
    }
    
    std::cout << "Best-fit measurement: (" << anue[ifbest] << ", " << bnue[ifbest] << ", " << cnue[ifbest] << ") with reduced chi2 = " << chi2min/dofbest << std::endl;
    
    int index = std::distance(chi2values_tmp.begin(), std::min_element(chi2values_tmp.begin(), chi2values_tmp.end()) );
    //std::cout << "Lowest chi2 = " << chi2values_tmp[index] << " " << iflux_tmp[index] << ": (" << anue[iflux_tmp[index]] << ", " << bnue[iflux_tmp[index]] << ", " << cnue[iflux_tmp[index]] << ")" << std::endl;
    
    //now I want to make plots of chi2 vs parameter for the other two parameters
    //at their fixed truth values!
    std::vector<double> alphachi2;
    std::vector<double> alphavals;
    TString titlealpha("Minimum reduced #chi^{2} vs. #alpha;#alpha;Reduced #chi^{2}");
    TString namealpha("Chi2VsAlpha");
    std::vector<double> e0chi2;
    std::vector<double> e0vals;
    TString titlee0("Minimum reduced #chi^{2} vs. #LT E_{#nu} #GT;#LT E_{#nu} #GT (MeV);Reduced #chi^{2}");
    TString namee0("Chi2VsE0");
    std::vector<double> lumchi2;
    std::vector<double> lumvals;
    TString titlelum("Minimum reduced #chi^{2} vs. #varepsilon;#varepsilon (10^{53} erg);Reduced #chi^{2}");
    TString namelum("Chi2VsLum");
    
    //we also want to keep the alpha, e0, lum values corresponding to "good" chi2 values
    TH2D *hAllowedRegion_AlphaVsE0 = new TH2D("alphavse0", ";#LT E_{#nu} #GT (MeV);#alpha", nume0bins, mine0, maxe0, numalphabins, minalpha, maxalpha );
    TH2D *hAllowedRegion_LumVsE0 = new TH2D("lumvse0", ";#LT E_{#nu} #GT (MeV);#varepsilon (10^{53} erg)", nume0bins, mine0, maxe0, numlumbins, minlum/factor, maxlum/factor );
    TH2D *hAllowedRegion_LumVsAlpha = new TH2D("lumvsalpha", ";#alpha;#varepsilon (10^{53} erg)", numalphabins, minalpha, maxalpha, numlumbins, minlum/factor, maxlum/factor );
    
    TH2D *hAllowedRegion_AlphaVsE0_90Cut = new TH2D("alphavse0_90cut", ";#LT E_{#nu} #GT (MeV);#alpha", nume0bins, mine0, maxe0, numalphabins, minalpha, maxalpha );
    TH2D *hAllowedRegion_LumVsE0_90Cut = new TH2D("lumvse0_90cut", ";#LT E_{#nu} #GT (MeV);#varepsilon (10^{53} erg)", nume0bins, mine0, maxe0, numlumbins, minlum/factor, maxlum/factor );
    TH2D *hAllowedRegion_LumVsAlpha_90Cut = new TH2D("lumvsalpha_90cut", ";#alpha;#varepsilon (10^{53} erg)", numalphabins, minalpha, maxalpha, numlumbins, minlum/factor, maxlum/factor );
    
    //set to zero
    tr_numGoodGridElementsAlphavsE0 = 0;
    tr_numGoodGridElementsLumvsE0 = 0;
    tr_numGoodGridElementsLumvsAlpha = 0;
    
    for(Double_t iAlpha = first_alpha; iAlpha <= last_alpha+step_alpha; iAlpha += step_alpha){
        
        //iAlpha tells us the alpha value we care about
        std::vector<double> chi2vals;
        
        for(iFlux = 0; iFlux < numfluxes; ++iFlux){
            //check to make sure we're at the right alpha value, and that alpha/E0 are physical
            
            Int_t diff = ( anue[iFlux] - iAlpha )*100.0;
            
            if(diff == 0.0 && anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0){
                
                chi2vals.emplace_back(chi2values[iFlux]);
            }
        }
        
        //std::cout << iAlpha << " " << chi2vals.size() << std::endl;
        
        //now we find the minimum chi2 from these values
        if(chi2vals.size() > 0){
            size_t iSmallest = std::distance(chi2vals.begin(), std::min_element( chi2vals.begin(), chi2vals.end() ));
            
            alphavals.emplace_back( iAlpha );
            
            alphachi2.emplace_back( chi2vals[iSmallest] );
            
            //std::cout << iAlpha << " " << chi2vals[iSmallest] << std::endl;
            
        }
        
        //2D plot
        for(Double_t iE0 = first_e0; iE0 <= last_e0; iE0 += step_e0){
            
            std::vector<double> chi2vals2d;
            
            for(iFlux = 0; iFlux < numfluxes; ++iFlux){
                //check to make sure we're at the right alpha value, and that alpha/E0 are physical
                
                Int_t diffalpha = ( anue[iFlux] - iAlpha )*100.0;
                Int_t diffe0 = ( bnue[iFlux] - iE0 )*100.0;
                
                if(diffalpha == 0.0 && diffe0 == 0.0 && anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0){
                    chi2vals2d.emplace_back(chi2values[iFlux]);
                }
                
            }
            
            //std::cout << iAlpha << " " << iE0 << " " << chi2vals2d.size() << std::endl;
            
            if(chi2vals2d.size() > 0){
                size_t iSmallest = std::distance(chi2vals2d.begin(), std::min_element( chi2vals2d.begin(), chi2vals2d.end() ));
                
                //std::cout << iAlpha << " " << iE0 << " " << chi2vals2d[iSmallest] << std::endl;
                Double_t toFill = chi2vals2d[iSmallest];
                if(toFill == 0.0) toFill = 1e-20;
                
                if(toFill <= chi2_90contour){ 
					++tr_numGoodGridElementsAlphavsE0;
					hAllowedRegion_AlphaVsE0_90Cut->Fill(iE0, iAlpha, toFill);
                }
                
                hAllowedRegion_AlphaVsE0->Fill(iE0, iAlpha, toFill);
                
            }
            
            chi2vals2d = std::vector<double>();
            
        }
        
        //2D plot:
        for(Double_t iLum = first_lum; iLum <= last_lum+step_lum; iLum += step_lum){
            std::vector<double> chi2vals2d;
            
            for(iFlux = 0; iFlux < numfluxes; ++iFlux){
                //check to make sure we're at the right alpha value, and that alpha/E0 are physical
                Int_t diffalpha = ( anue[iFlux] - iAlpha )*100.0;
                Int_t difflum = ( cnue[iFlux]/factor - iLum/factor )*100.0;
                
                if(diffalpha == 0.0 && difflum == 0.0 && anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0){
                    chi2vals2d.emplace_back(chi2values[iFlux]);
                }
                
            }
            
            
            if(chi2vals2d.size() > 0){
                size_t iSmallest = std::distance(chi2vals2d.begin(), std::min_element( chi2vals2d.begin(), chi2vals2d.end() ));
                
                
                Double_t toFill = chi2vals2d[iSmallest];
                if(toFill == 0.0) toFill = 1e-20;
                
                if(toFill <= chi2_90contour){ 
					++tr_numGoodGridElementsLumvsAlpha;
					hAllowedRegion_LumVsAlpha_90Cut->Fill(iAlpha, iLum/factor, toFill);
				}
                
                hAllowedRegion_LumVsAlpha->Fill(iAlpha, iLum/factor, toFill);
                
            }
            
            chi2vals2d = std::vector<double>();
        }
        
        chi2vals = std::vector<double>();
        
    }
    
    for(Double_t iE0 = first_e0; iE0 <= last_e0; iE0 += step_e0){
        
        //iE0 tells us the E0 value we care about
        std::vector<double> chi2vals;
        
        for(iFlux = 0; iFlux < numfluxes; ++iFlux){
            //check to make sure we're at the right alpha value, and that alpha/E0 are physical
            Int_t diffe0 = ( bnue[iFlux] - iE0 )*100.0;
            
            if(diffe0 == 0.0 && anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0){
                chi2vals.emplace_back(chi2values[iFlux]);
            }
        }
        
        //now we find the minimum chi2 from these values
        if(chi2vals.size() > 0){
            size_t iSmallest = std::distance(chi2vals.begin(), std::min_element( chi2vals.begin(), chi2vals.end() ));
            
            e0vals.emplace_back( iE0 );
            
            e0chi2.emplace_back( chi2vals[iSmallest] );
            
            //std::cout << iE0 << " " << chi2vals[iSmallest] << std::endl;
            
        }
        
        //2D plot
        for(Double_t iLum = first_lum; iLum <= last_lum+step_lum; iLum += step_lum){
            std::vector<double> chi2vals2d;
            
            for(iFlux = 0; iFlux < numfluxes; ++iFlux){
                //check to make sure we're at the right e0 value, and that alpha/E0 are physical
                Int_t diffe0 = ( bnue[iFlux] - iE0 )*100.0;
                Int_t difflum = ( cnue[iFlux]/factor - iLum/factor )*100.0;
                if(diffe0 == 0.0 && difflum == 0.0 && anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0){
                    chi2vals2d.emplace_back(chi2values[iFlux]);
                }
            }
            
            if(chi2vals2d.size() > 0){
                size_t iSmallest = std::distance(chi2vals2d.begin(), std::min_element( chi2vals2d.begin(), chi2vals2d.end() ));
                
                //std::cout << iE0 << " " << iLum << " " << chi2vals2d[iSmallest] << std::endl;
                
                Double_t toFill = chi2vals2d[iSmallest];
                if(toFill == 0.0) toFill = 1e-20;
                
                if(toFill <= chi2_90contour){ 
					++tr_numGoodGridElementsLumvsE0;
					hAllowedRegion_LumVsE0_90Cut->Fill(iE0, iLum/factor, toFill);
				}
                
                hAllowedRegion_LumVsE0->Fill(iE0, iLum/factor, toFill);
                
            }
            
            chi2vals2d = std::vector<double>();
        }
        
        chi2vals = std::vector<double>();
        
    }
    
    for(Double_t iLum = first_lum; iLum <= last_lum+step_lum; iLum += step_lum){
        
        //iLum tells us the luminosity value we care about
        std::vector<double> chi2vals;
        
        for(iFlux = 0; iFlux < numfluxes; ++iFlux){
            
            //check to make sure we're at the right lum value, and that alpha/E0 are physical
            Int_t difflum = ( cnue[iFlux]/factor - iLum/factor )*100.0;
            //if(cnue[iFlux] == iLum){
            if(difflum == 0.0){
                
                if(anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0) chi2vals.emplace_back(chi2values[iFlux]);
            }
        }
        
        //now we find the minimum chi2 from these values
        if(chi2vals.size() > 0){
            size_t iSmallest = std::distance(chi2vals.begin(), std::min_element( chi2vals.begin(), chi2vals.end() ));
            
            lumvals.emplace_back( iLum/factor );
            
            lumchi2.emplace_back( chi2vals[iSmallest] );
            
            //std::cout << iLum << " " << chi2vals[iSmallest] << std::endl;
            
        }
        
        chi2vals = std::vector<double>();
    }
    
    TGraph *g_Chi2VsAlpha = makeTGraphFromVectors(alphavals, alphachi2, titlealpha);
    TGraph *g_Chi2VsE0 = makeTGraphFromVectors(e0vals, e0chi2, titlee0);
    TGraph *g_Chi2VsLum = makeTGraphFromVectors(lumvals, lumchi2, titlelum);
    
    hAllowedRegion_AlphaVsE0->GetXaxis()->CenterTitle();
    hAllowedRegion_AlphaVsE0->GetXaxis()->SetTitleSize(0.05);
    hAllowedRegion_AlphaVsE0->GetYaxis()->CenterTitle();
    hAllowedRegion_AlphaVsE0->GetYaxis()->SetTitleSize(0.05);
    
    hAllowedRegion_LumVsE0->GetXaxis()->CenterTitle();
    hAllowedRegion_LumVsE0->GetXaxis()->SetTitleSize(0.04);
    hAllowedRegion_LumVsE0->GetYaxis()->CenterTitle();
    hAllowedRegion_LumVsE0->GetYaxis()->SetTitleSize(0.04);
    
    hAllowedRegion_LumVsAlpha->GetXaxis()->CenterTitle();
    hAllowedRegion_LumVsAlpha->GetXaxis()->SetTitleSize(0.05);
    hAllowedRegion_LumVsAlpha->GetYaxis()->CenterTitle();
    hAllowedRegion_LumVsAlpha->GetYaxis()->SetTitleSize(0.05);
    
    //change boundaries
    
    hAllowedRegion_AlphaVsE0->GetXaxis()->SetRangeUser(mine02, maxe02);
    hAllowedRegion_AlphaVsE0->GetYaxis()->SetRangeUser(minalpha2, maxalpha2);
    
    hAllowedRegion_LumVsE0->GetXaxis()->SetRangeUser(mine02, maxe02);
    hAllowedRegion_LumVsE0->GetYaxis()->SetRangeUser(minlum2/factor, maxlum2/factor);
    
    hAllowedRegion_LumVsAlpha->GetXaxis()->SetRangeUser(minalpha2, maxalpha2);
    hAllowedRegion_LumVsAlpha->GetYaxis()->SetRangeUser(minlum2/factor, maxlum2/factor);
    
    hAllowedRegion_AlphaVsE0_90Cut->GetXaxis()->SetRangeUser(mine02, maxe02);
    hAllowedRegion_AlphaVsE0_90Cut->GetYaxis()->SetRangeUser(minalpha2, maxalpha2);
    
    hAllowedRegion_LumVsE0_90Cut->GetXaxis()->SetRangeUser(mine02, maxe02);
    hAllowedRegion_LumVsE0_90Cut->GetYaxis()->SetRangeUser(minlum2/factor, maxlum2/factor);
    
    hAllowedRegion_LumVsAlpha_90Cut->GetXaxis()->SetRangeUser(minalpha2, maxalpha2);
    hAllowedRegion_LumVsAlpha_90Cut->GetYaxis()->SetRangeUser(minlum2/factor, maxlum2/factor);
    
    //limit 1D plots to ~5 sigma
    //if(g_Chi2VsAlpha->GetYaxis()->GetXmax() > 25.0) g_Chi2VsAlpha->SetMaximum(25.0);
    //if(g_Chi2VsE0->GetYaxis()->GetXmax() > 25.0) g_Chi2VsE0->SetMaximum(25.0);
    //if(g_Chi2VsLum->GetYaxis()->GetXmax() > 25.0) g_Chi2VsLum->SetMaximum(25.0);
    
    //set names
    g_Chi2VsAlpha->SetName(namealpha);
    g_Chi2VsE0->SetName(namee0);
    g_Chi2VsLum->SetName(namelum);
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // SAVE PLOTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
	//fill tree here
    tr_alpha = anue[ifbest];
    tr_e0 = bnue[ifbest];
    tr_lum = cnue[ifbest]/factor;
    tr_chi2 = chi2min;
    tr_dof = dofbest;
    tr_numGridElementsAlphavsE0 = numalphabins*nume0bins;
    tr_numGridElementsLumvsE0 = numlumbins*nume0bins;
    tr_numGridElementsLumvsAlpha = numalphabins*numlumbins;
    tr->Fill();
    
    TFile *fout = new TFile(outfile, "recreate");
    tr->Write();
    g_Chi2VsAlpha->Write();
    g_Chi2VsE0->Write();
    g_Chi2VsLum->Write();
    hAllowedRegion_AlphaVsE0->Write();
    hAllowedRegion_AlphaVsE0_90Cut->Write();
    hAllowedRegion_LumVsE0->Write();
    hAllowedRegion_LumVsE0_90Cut->Write();
    hAllowedRegion_LumVsAlpha->Write();
    hAllowedRegion_LumVsAlpha_90Cut->Write();
    fout->Close();
    
    //delete
    delete test_spectrum;
    delete hAllowedRegion_AlphaVsE0;
    delete hAllowedRegion_AlphaVsE0_90Cut;
    delete hAllowedRegion_LumVsE0;
    delete hAllowedRegion_LumVsE0_90Cut;
    delete hAllowedRegion_LumVsAlpha;
    delete hAllowedRegion_LumVsAlpha_90Cut;
	delete g_Chi2VsAlpha;
    delete g_Chi2VsE0;
    delete g_Chi2VsLum;
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
    chi2values = std::vector<Double_t>();
    alphavals = std::vector<double>();
    alphachi2 = std::vector<double>();
	e0vals = std::vector<double>();
	e0chi2 = std::vector<double>(); 
	lumvals = std::vector<double>(); 
	lumchi2 = std::vector<double>();
    
    //add line to separate the terminal output(s)
    std::cout << "--------------------------------" << std::endl;
    
}//end function

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// SaveBestFitElement: Produces chi2-minimized flux parameter measurements and
//                     plots for a supernova located distanceSN from Earth.
//                     Saves test spectrum and best-fit plot in a file.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void SaveBestFitElement(Double_t distanceSN, TString grid, TString grid_smear, TString grid_xscn, TString grid_effic, TString ts_smear, TString ts_xscn, TString ts_effic){
	//grid element we care about
    Double_t true_alpha = 2.5;
    Double_t true_e0 = 9.5;
    Double_t true_lum = 0.5;
    
    std::cout << "Producing best-fit plot for a " << distanceSN << "kpc supernova and the following input parameters:" << std::endl;
    std::cout << "	Grid: " << grid << " with smearing " << grid_smear << ", xscn " << grid_xscn << ", and effic " << grid_effic << std::endl;  
    std::cout << "	Test spectrum: " << ts_smear << ", xscn " << ts_xscn << " and effic " << ts_effic << std::endl;    
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DEFINE INPUTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //define string to hold grid parameters
    TString grid_pars = "smear_" + grid_smear + "_" + grid_xscn + "_" + grid_effic;
    TString ts_pars = ts_smear + "_" + ts_xscn + "_" + ts_effic;
	TString test_spect_filename = "input/test_spectra/pinched_test_smeared_sum_" + ts_pars + ".dat";
        
    TString smeareddir = "input/" + grid + "/" + grid_pars;    
    TString pinchedinfo = "input/pinched_info/pinched_info_" + grid + ".dat";
    TString gridinfo = "input/grid_info/grid_info_" + grid + ".dat";
    //define output file
    TString distanceStr; distanceStr.Form("%.2lf", distanceSN);
    
    //output file name
    TString outfilename = "asimov_bestfit_" + grid_pars + "_" + ts_pars + ".root";
    
    Double_t massfact = 1.0;//2.4; //see the full DUNE far detector response
    Double_t chi2_90contour = 4.61;  //pdg 2018
    
	//use same style settings that Kate used
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1,0);
    //set random seed
    gRandom->SetSeed(0);
    
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
    
    Double_t minalpha2 = first_alpha - step_alpha*2.0;
    Double_t maxalpha2 = last_alpha + step_alpha*2.0;
    Double_t mine02 = first_e0 - step_e0*2.0;
    Double_t maxe02 = last_e0 + step_e0*2.0;
    Double_t minlum2 = first_lum - step_lum*2.0;
    Double_t maxlum2 = last_lum + step_lum*2.0;
    
    Int_t numalphabins = int(last_alpha - first_alpha)/step_alpha+1;
    Int_t nume0bins = int(last_e0 - first_e0)/step_e0+1;
    Int_t numlumbins = (last_lum - first_lum)/(step_lum) + 1;
    
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
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // READ PINCHING PARAMETERS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
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
    
    std::cout << "This grid contains " << ip << " fluxes." << std::endl;
    
    std::vector<TH1D*> histos; 
    
    TH1D* test_spectrum = fill_spect_hist(test_spect_filename,"test_spectrum",distanceSN, massfact);
    
    //TString tslabel = "Test Spectrum: (" + TString::Format("%.1lf", true_alpha) + ", " + TString::Format("%.1lf", true_e0) + ", " + TString::Format("%.2lf", true_lum) + "e53)";
	histos.emplace_back(test_spectrum);
    
    //now we loop over templates in the grid
    Int_t iFlux; //iterator over the grid
    Int_t numfluxes = ip; //total number of fluxes we care about
    
    const Int_t maxhist = maxflux;
    
    Double_t chi2min = 100000000000.0;
    Int_t dofbest;
    Int_t ifbest;
    
    Double_t chi2;
    Int_t dof;
    
    for(iFlux = 0; iFlux < numfluxes; ++iFlux){
        
        // Do not go outside the limits
        if(anue[iFlux] <= maxgoodalpha && bnue[iFlux] <= maxgoode0 && anue[iFlux] >= mingoodalpha && bnue[iFlux] >= mingoode0 ){
            
            //define filename
            TString pargridfilename = smeareddir+"/pinched_"+TString::Format("%d",iFlux)+"_smeared_sum.dat";
            TString histname = "pinched_"+grid + "_" + grid_smear + "_" + grid_xscn + "_" + grid_effic + "_" + ts_smear + "_" + ts_xscn + "_" + ts_effic + "_" + distanceStr + "kpc_"+TString::Format("%d",iFlux);
            
            TH1D *pargridhist = fill_spect_hist(pargridfilename, histname, distanceSN, massfact);
            mychi2(test_spectrum,pargridhist,&chi2,&dof);
            
            
            if (chi2<chi2min) {
                chi2min = chi2;
                dofbest = dof;
                ifbest = iFlux;
            }
            
        }
    }

	std::cout << "Best-fit measurement: (" << anue[ifbest] << ", " << bnue[ifbest] << ", " << cnue[ifbest] << ") with reduced chi2 = " << chi2min/dofbest << std::endl;    
    
	TString bffilename = smeareddir+"/pinched_"+TString::Format("%d",ifbest)+"_smeared_sum.dat";
	TString bfname = "best_fit_" + TString::Format("%.1lf", anue[ifbest]) + "_" + TString::Format("%.1lf", bnue[ifbest]) + "_" + TString::Format("%.2lf", cnue[ifbest]/factor) + "e53";
	TH1D *bfhist = fill_spect_hist(bffilename, bfname, distanceSN, massfact);
	
	histos.emplace_back(bfhist);
	
    
    //save in ROOT file
    TFile *fout = new TFile(outfilename, "recreate");
    for(size_t i = 0; i < histos.size(); ++i) histos[i]->Write();
    fout->Close();
    
    //delete
    delete test_spectrum;
    delete bfhist;
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
    
    //add line to separate the terminal output(s)
    std::cout << "--------------------------------" << std::endl;
}

#endif /* Asimov_h */
