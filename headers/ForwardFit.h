//
//  ForwardFit.h
//  Erin Conley (erin.conley@duke.edu)
//

#ifndef ForwardFit_h
#define ForwardFit_h

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// getReferenceContourAreas: Open the file with reference contour areas and
//                           store the values in ref_alphavse0, ref_lumvse0, and
//                           ref_lumvsalpha
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void getReferenceContourAreas(Double_t distanceSN, TString grid, TString grid_xscn, TString grid_effic, TString ts_xscn, TString ts_effic, Double_t *ref_alphavse0, Double_t *ref_lumvse0, Double_t *ref_lumvsalpha){
    
    TString ref_filename = "input/reference_areas/referencearea_" + grid + "_smear0.00MARLEY_" + grid_xscn + "_" + grid_effic + "_spectra0.00MARLEY_" + ts_xscn + "_" + ts_effic + "_" + TString::Format("%.2lf", distanceSN) + "kpc.dat";
    
    ifstream refin;
    refin.open(ref_filename);
    while(1){
        
        refin >> *ref_alphavse0 >> *ref_lumvse0 >> *ref_lumvsalpha;
        
        if(!refin.good()) break;
        //std::cout << *ref_alphavse0 << " " << *ref_lumvse0 << " " << *ref_lumvsalpha << std::endl;
        
    }
    
    refin.close();
    
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// findCorrectGridForCombo: For a given list of grids, search for ROOT file of
//                          chi2 plots and contours and output the file name
//                          if the file exists
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TString findCorrectGridForCombo(TString directory, std::vector<TString> grids, TString grid_options, TString ts_options, TString distanceSN){
    
    //this needs to find the correct file name for a certain grid/ts combination
    TString actualfilename = "NULL";
    
    for(size_t i = 0; i < grids.size(); ++i){
        TString filename = directory + "chi2plots_" + grids[i] + "_smear" + grid_options + "_spectra" + ts_options + "_" + distanceSN + ".root";
        //check if this filename exists
        
        
        ifstream ifile(filename);
        if(ifile) return filename;
        
        
    }
    
    return actualfilename;
    
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// getMinOrMaxFromHists: For a given vector of histograms, find the minimum and
//                       and maximum x- and y-values.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void getMinOrMaxFromHists(std::vector<TH2D*> hists, double &xmin, double &xmax, double &ymin, double &ymax){
    
    std::vector<double> allxmins;
    std::vector<double> allxmaxes;
    std::vector<double> allymins;
    std::vector<double> allymaxes;
    
    for(size_t i = 0; i < hists.size(); ++i){
        
        //get the min and max of the
        allxmins.emplace_back( hists[i]->GetXaxis()->GetXmin() );
        allxmaxes.emplace_back( hists[i]->GetXaxis()->GetXmax() );
        allymins.emplace_back( hists[i]->GetYaxis()->GetXmin() );
        allymaxes.emplace_back( hists[i]->GetYaxis()->GetXmax() );
        
    }
    
    xmin = *std::min_element(allxmins.begin(),allxmins.end());
    xmax = *std::max_element(allxmaxes.begin(),allxmaxes.end());
    ymin = *std::min_element(allymins.begin(),allymins.end());
    ymax = *std::max_element(allymaxes.begin(),allymaxes.end());
    
    
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// makeTGraphFromVectors: Create a ROOT TGraph object from given vectors of
//                        values and a title.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TGraph* makeTGraphFromVectors(std::vector<double> xAxis, std::vector<double> yAxis, TString title){
	
	//first make the arrays
	Double_t* xArray = new Double_t[xAxis.size()];
	Double_t* yArray = new Double_t[yAxis.size()];
	
	for(size_t i = 0; i < xAxis.size(); ++i){
		xArray[i] = xAxis[i];
		yArray[i] = yAxis[i];
	}
	
	TGraph *g = new TGraph(xAxis.size(), xArray, yArray);
	g->SetTitle(title);
	g->SetMarkerStyle(20);
	g->SetMarkerSize(1.0);
	g->GetYaxis()->SetTitleOffset(1.2);
	g->SetLineColor(kWhite);
	
	return g;
	
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// getParameter: Read a TTree object for a given ROOT file defined by filename
//               and a TBranch name. Outputs a vector of values.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

std::vector<double> getParameter(TString filename, TString parameterName){
    
    std::vector<double> values;
    
    TFile *f = new TFile(filename);
    TTree *tr = (TTree*)gDirectory->Get("data;1");
    TTreeReader read(tr);
    TTreeReaderValue<double> val(read, parameterName);
    
    while(read.Next()) values.emplace_back(*val);
    
    f->Close();
    
    return values;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// getManyParameters: For a list of file names, outputs a 2D vector of values
//                    from TTree objects in the files.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

std::vector<std::vector<double>> getManyParameters(std::vector<TString> files, TString parameterName){
    
    std::vector<std::vector<double>> values;
    
    for(size_t i = 0; i < files.size(); ++i){
        std::vector<double> tmp_values = getParameter(files[i], parameterName);
        
        values.emplace_back(tmp_values);
    }
    
    return values;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// whoWins: Given a list of TH1D histograms, return the index corresponding to
//          the histogram with the highest maximum bin content.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

int whoWins(std::vector<TH1D*> hists){
	std::vector<int> bincontents;
	
	for(std::size_t i = 0; i < hists.size(); ++i){
		int bc = hists[i]->GetBinContent(hists[i]->GetMaximumBin());
		bincontents.emplace_back(bc);
	}
	
	return std::distance(bincontents.begin(), std::max_element(bincontents.begin(), bincontents.end()));
	
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// makeLegend: Creates a TLegend object given four values; x1,y1,x2,y2 are the
//             coordinates of the Legend in the current pad (in normalised
//             coordinates by default)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TLegend* makeLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2){
    //this function makes a legend
    TLegend *leg = new TLegend(x1, y1, x2, y2);
    leg->SetFillColor(10);
    leg->SetTextFont( (Font_t) 62);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    
    return leg;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// superimposeTH1D: For a given vector of TH1D histograms, legend labels, and a
//                  plot title, superimpose the histograms and add a legend.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void superimposeTH1D(std::vector<TH1D*> hists, std::vector<TString> names, TString title){
	
	TCanvas *c = new TCanvas(title, title, 1000, 700);
	gStyle->SetOptStat(0);
	
	for(size_t i = 0; i < hists.size(); ++i){
		
		//hists[i]->SetLineWidth(0);
		hists[i]->SetMarkerStyle(20);
		//hists[i]->SetMarkerSize(0.75);
		
		if(i+1==5){ 
			hists[i]->SetLineColor(i+2);
			hists[i]->SetMarkerColor(i+2);
		}
		else if(i+1 == 3){ 
			hists[i]->SetLineColor(kGreen+2);
			hists[i]->SetMarkerColor(kGreen+2);
		}
		else{ 
			hists[i]->SetLineColor(i+1);
			hists[i]->SetMarkerColor(i+1);
		}
	}
	
	TLegend *leg = makeLegend(0.3, 0.6, 0.7, 0.85);
	
	//change this later
	leg->SetTextSize(0.03);
	leg->SetBorderSize(0);
	
	for(size_t i = 0; i < hists.size(); ++i) leg->AddEntry(hists[i], names[i], "pe");
	
	int toDrawFirst = whoWins(hists);
	
	hists[toDrawFirst]->SetTitle(title);
	hists[toDrawFirst]->Draw();
	
	for(size_t i = 0; i < hists.size(); ++i){
		if(i != toDrawFirst) hists[i]->Draw("SAME");
	}
	leg->Draw("SAME");
	
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// superimposeTGraph: For a given vector of TGraph objects, legend labels, and
//                    plot title, superimpose the TGraphs and add a legend.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void superimposeTGraph(std::vector<TGraph*> graphs, std::vector<TString> names, TString title){
	
	TCanvas *c = new TCanvas(title, title, 1000, 700);
	gStyle->SetOptStat(0);
	
	for(size_t i = 0; i < graphs.size(); ++i){
		graphs[i]->SetLineColor(5);
		if(i+1==5){ 
			graphs[i]->SetLineColor(i+2);
			graphs[i]->SetMarkerColor(i+2);
		}
		else if(i+1 == 3){ 
			graphs[i]->SetLineColor(kGreen+2);
			graphs[i]->SetMarkerColor(kGreen+2);
		}
		else{ 
			graphs[i]->SetLineColor(i+1);
			graphs[i]->SetMarkerColor(i+1);
		}
		
	}
	
	TLegend *leg = makeLegend(0.3, 0.6, 0.7, 0.85);
	for(size_t i = 0; i < graphs.size(); ++i) leg->AddEntry(graphs[i], names[i], "p");
	leg->SetTextSize(0.03);
	leg->SetBorderSize(0);
	
	auto mg = new TMultiGraph();
	for(size_t i = 0; i < graphs.size(); ++i) mg->Add(graphs[i]);
	
	mg->Draw("AL");
	mg->SetTitle(title);
	mg->GetXaxis()->SetTitleSize(0.04);
	mg->GetXaxis()->SetLabelSize(0.04);
	mg->GetYaxis()->SetTitleSize(0.04);
	mg->GetYaxis()->SetTitleOffset(1.);
	mg->GetYaxis()->SetLabelSize(0.04);
	//mg->SetMinimum(0.0);
	mg->GetXaxis()->SetRangeUser(5.5, 45.0);
	//if(mg->GetYaxis()->GetXmax() > 25.0) mg->SetMaximum(25.0);	
	gPad->Modified();
	leg->Draw("SAME");
	
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// makeTGraphFromFilename: Reads in event rates contained in filename, and
//                         returns a TGraph named histname corresponding to SN
//                         distance and mass factor (to scale to specific
//                         detector volume)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TGraph* makeTGraphFromFilename(TString filename, TString histname, Double_t distance, Double_t massfactor) {
	//TGraph interpolation 
    // Read in the data from the file
    ifstream in;
    in.open(filename);
    
    std::vector<double> energyVector; std::vector<double> nueVector;
    
    while (1) {
		Double_t en, nue;
        in >> en>>nue;
        
        // Energy in file is in eV
        en *= 1000.;  // in MeV for the plots
        
        //now scale for the distance we want
        nue *= std::pow(10.0/distance, 2);
        //and scale for the mass we want
        nue *= massfactor;                
                        
        if (!in.good()) break;
		
		energyVector.emplace_back(en);
		nueVector.emplace_back(nue);
		        
    }
    
    TGraph *g = makeTGraphFromVectors(energyVector, nueVector, histname);
    
    return g;
    
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// scale_spect_hist: Reads in event rates contained in filename, and returns a
//                   histogram named histname corresponding to a given supernova
//                   distance and mass factor (to scale to specific detector
//                   volume). Change the scaling factor "scale" to shift the
//                   energy scale; typical scale range (-100, 100).
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TH1D* scale_spect_hist(TString filename, Double_t scale, TString histname, Double_t distance, Double_t massfactor) {
	//TGraph interpolation 
    // Read in the data from the file
    ifstream in;
    in.open(filename);
    
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
    
    TGraph *g = makeTGraphFromVectors(energyVector, nueVector, histname);
    
    TH1D* spectrum = new TH1D(histname,"Spectrum",200, 0.5, 100.);
    
    for(int i = 0; i < g->GetN(); ++i){
		Double_t en = g->GetX()[i];
				
		Double_t enShift = en*(1.0+ (scale/100.0) ); 
		Double_t evShift = g->Eval(enShift);
		
		//watch out for negative event rates!
		if(evShift < 0.0) evShift = 0.0;
		//also watch out for when the scaling is zero 
		//don't change the event rate!
		if(scale == 0.0) evShift = g->GetY()[i];
		
		// File has actually events per bin
        Double_t binfact=1;
        evShift /=  binfact;
        
        //now scale for the distance we want
        evShift *= std::pow(10.0/distance, 2);
        //and scale for the mass we want
        evShift *= massfactor;
		
		spectrum->Fill(en, evShift);
	}
    
    for(int j = 1; j <= spectrum->GetNbinsX(); ++j){
        spectrum->SetBinError(j, std::sqrt(spectrum->GetBinContent(j)) );
    }
    
    //Int_t numpoints=i;
    
    //  printf(" found %d points\n",numpoints);
    
    in.close();
    
    // Rebin-- Snowglobes output bin boundaries are on a slight "slant"
    //5 = merges five bins in one; so output will have 40 bins
    //"" = histogram is modified
    //0 = means we are not defining new bin low-edges
    spectrum->Rebin(5,"",0);
    
    return spectrum;
    
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// shift_spect_hist: Reads in event rates contained in filename, and returns a
//                   histogram named histname corresponding to a given supernova
//                   distance and mass factor (to scale to specific detector
//                   volume). Change the shift factor "shift" (MeV) to set the
//                   detection threshold of the histogram
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TH1D* shift_spect_hist(TString filename, TString histname, Double_t distance, Double_t massfactor, Double_t shift) {
    
    // Fill a spectrum histogram from a filename and return a pointer to it
    // returns a histogram
    
    TH1D* spectrum = new TH1D(histname,"Spectrum",200, 0.5, 100.);
    
    // Read in the data from the file
    ifstream in;
    in.open(filename);
    
    const Int_t maxpoints = 1000;
    //Double_t en[maxpoints],nue[maxpoints];
    std::vector<Double_t> en; en.resize(maxpoints);
    std::vector<Double_t> nue; nue.resize(maxpoints);
    
    Int_t i=0;
    while (1) {
        in >> en[i]>>nue[i];
        
        // Energy in file is in eV
        en[i] *= 1000.;  // in MeV for the plots
        
        // File has actually events per bin
        
        Double_t binfact=1;
        nue[i] /=  binfact;
        
        //now scale for the distance we want
        nue[i] *= std::pow(10.0/distance, 2);
        //and scale for the mass we want
        nue[i] *= massfactor;
        
        //    cout <<i<<" "<<en[i]<<" "<<nue[i]<<endl;
        
        //spectrum->Fill(en[i],nue[i]);
        if(en[i] >= shift) spectrum->Fill(en[i],nue[i]);
        else{ 
			
			//if(nue[i] > 0.0) std::cout << "Energy " << en[i] << " cut out." << std::endl;
			
			spectrum->Fill(en[i],0.0);
		}
        
        if (!in.good()) break;
        i++;
        
    }
    
	for(int j = 1; j <= spectrum->GetNbinsX(); ++j){
		spectrum->SetBinError(j, std::sqrt(spectrum->GetBinContent(j)) );
    }
    
    Int_t numpoints=i;
    
    //  printf(" found %d points\n",numpoints);
    
    in.close();
    
    // Rebin-- Snowglobes output bin boundaries are on a slight "slant"
    //5 = merges five bins in one; so output will have 40 bins
    //"" = histogram is modified
    //0 = means we are not defining new bin low-edges
    spectrum->Rebin(5,"",0);
    
    return spectrum;
    
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// fill_spect_hist: Reads in event rates contained in filename, and returns a
//                  histogram named histname corresponding to a given supernova
//                  distance and mass factor (to scale to specific detector
//                  volume).
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TH1D* fill_spect_hist(TString filename, TString histname, Double_t distance, Double_t massfactor) {
    
    // Fill a spectrum histogram from a filename and return a pointer to it
    // returns a histogram
    
    TH1D* spectrum = new TH1D(histname,"Spectrum",200, 0.5, 100.);
    
    // Read in the data from the file
    ifstream in;
    in.open(filename);
    
    const Int_t maxpoints = 1000;
    //Double_t en[maxpoints],nue[maxpoints];
    std::vector<Double_t> en; en.resize(maxpoints);
    std::vector<Double_t> nue; nue.resize(maxpoints);
    
    Int_t i=0;
    while (1) {
        in >> en[i]>>nue[i];
        
        // Energy in file is in eV
        en[i] *= 1000.;  // in MeV for the plots
        
        // File has actually events per bin
        
        Double_t binfact=1;
        nue[i] /=  binfact;
        
        //now scale for the distance we want
        nue[i] *= std::pow(10.0/distance, 2);
        //and scale for the mass we want
        nue[i] *= massfactor;
        
        //    cout <<i<<" "<<en[i]<<" "<<nue[i]<<endl;
        
        spectrum->Fill(en[i],nue[i]);
        
        if (!in.good()) break;
        i++;
        
    }
    
	for(int j = 1; j <= spectrum->GetNbinsX(); ++j){
		spectrum->SetBinError(j, std::sqrt(spectrum->GetBinContent(j)) );
    }
    
    Int_t numpoints=i;
    
    //  printf(" found %d points\n",numpoints);
    
    in.close();
    
    // Rebin-- Snowglobes output bin boundaries are on a slight "slant"
    //5 = merges five bins in one; so output will have 40 bins
    //"" = histogram is modified
    //0 = means we are not defining new bin low-edges
    spectrum->Rebin(5,"",0);
    
    return spectrum;
    
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// fill_spect_hist_NoRebin: Reads in event rates contained in filename, and
//                          returns a histogram named histname corresponding to
//                          a given supernova distance and mass factor (to scale
//                          to specific detector volume). No rebinning
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TH1D* fill_spect_hist_NoRebin(TString filename, TString histname, Double_t distance, Double_t massfactor) {
    
    // Fill a spectrum histogram from a filename and return a pointer to it
    // returns a histogram
    
    TH1D* spectrum = new TH1D(histname,"Spectrum",200, 0.5, 100.);
    
    // Read in the data from the file
    ifstream in;
    in.open(filename);
    
    const Int_t maxpoints = 1000;
    //Double_t en[maxpoints],nue[maxpoints];
    std::vector<Double_t> en; en.resize(maxpoints);
    std::vector<Double_t> nue; nue.resize(maxpoints);
    
    Int_t i=0;
    while (1) {
        in >> en[i]>>nue[i];
        
        // Energy in file is in eV
        en[i] *= 1000.;  // in MeV for the plots
        
        // File has actually events per bin
        
        Double_t binfact=1;
        nue[i] /=  binfact;
        
        //now scale for the distance we want
        nue[i] *= std::pow(10.0/distance, 2);
        //and scale for the mass we want
        nue[i] *= massfactor;
        
        //    cout <<i<<" "<<en[i]<<" "<<nue[i]<<endl;
        
        spectrum->Fill(en[i],nue[i]);
        
        if (!in.good()) break;
        i++;
        
    }
    
	for(int j = 1; j <= spectrum->GetNbinsX(); ++j){
		spectrum->SetBinError(j, std::sqrt(spectrum->GetBinContent(j)) );
    }
    
    Int_t numpoints=i;
    
    //  printf(" found %d points\n",numpoints);
    
    in.close();
    
    // Rebin-- Snowglobes output bin boundaries are on a slight "slant"
    //5 = merges five bins in one; so output will have 40 bins
    //"" = histogram is modified
    //0 = means we are not defining new bin low-edges
    //spectrum->Rebin(5,"",0);
    
    return spectrum;
    
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// mychi2: Calculates chi2 and dof between two TH1D histograms. Tests hist1
//         against hist2.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void mychi2(TH1D* hist1, TH1D* hist2, Double_t* chi2, Int_t* dof ) {
    
    // Test first histo against second
    // The histograms should have the same binning
    
    Double_t h1norm = hist1->Integral();
    Double_t h2norm = hist2->Integral();
    
      //cout << "h1norm, h2norm "<<h1norm<<" "<<h2norm<<endl;
    
    // Scale the second histogram to the norm of the first
    //hist2->Scale(h1norm/h2norm);
    
    //  cout << "New norm "<<hist2->Integral()<<endl;
    
    //setup to find chi2
    *chi2 = 0;
    Int_t numbin = hist1->GetNbinsX();
    Int_t i;
    Int_t numgoodbins = 0;
    
    //loop through the bins
    //for (i=0;i<numbin;i++) {
    for (i=1;i<=numbin;++i) {
        Double_t ni1 = hist1->GetBinContent(i);
        Double_t ni2 = hist2->GetBinContent(i);
        Double_t diff = ni1-ni2;
        if (ni1>0 && hist1->GetBinContent(i)>1) {
            //    if (ni1>0) {
            numgoodbins++;
            *chi2 += diff*diff/(hist1->GetBinContent(i));
            
        }
        
    }
    
    //Double_t chi2perdof = chi2/(numgoodbins-2);
    
    //set degrees of freedom
    *dof = numgoodbins-2;
    
}//end of mychi2

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// mychi2_distuncertainty: Calculates chi2 and dof between two TH1D histograms.
//                         Tests hist1 against hist2. Also incorporates
//                         uncertainty in supernova distance as a fraction
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void mychi2_distuncertainty(TH1D* hist1, TH1D* hist2, Double_t SNdistUncertainty, Double_t* chi2, Int_t* dof ) {
    
    // Test first histo against second
    // The histograms should have the same binning
    
    Double_t h1norm = hist1->Integral();
    Double_t h2norm = hist2->Integral();
    
      //cout << "h1norm, h2norm "<<h1norm<<" "<<h2norm<<endl;
    
    // Scale the second histogram to the norm of the first
    //hist2->Scale(h1norm/h2norm);
    
    //  cout << "New norm "<<hist2->Integral()<<endl;
    
    //setup to find chi2
    *chi2 = 0;
    Int_t numbin = hist1->GetNbinsX();
    Int_t i;
    Int_t numgoodbins = 0;
    
    //loop through the bins
    //for (i=0;i<numbin;i++) {
    for (i=1;i<=numbin;++i) {
        Double_t ni1 = hist1->GetBinContent(i);
        Double_t ni2 = hist2->GetBinContent(i);
        Double_t diff = ni1-ni2;
        if (ni1>0 && hist1->GetBinContent(i)>1) {
            //    if (ni1>0) {
            numgoodbins++;
            *chi2 += diff*diff/(hist1->GetBinContent(i)+SNdistUncertainty*SNdistUncertainty);
            
        }
        
    }
    
    //Double_t chi2perdof = chi2/(numgoodbins-2);
    
    //set degrees of freedom
    *dof = numgoodbins-2;
    
}//end of mychi2_distuncertainty

#endif /* ForwardFit_h */
