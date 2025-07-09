//**************************************
// Rosi Reed
// June 30, 2025
//
// This plot takes the CSV file downloaded from the Luminosity counting grafana
// GL1 & MBD_NS trigger rate panel
// Note: Go to inspect, data.  Then select Data Options and Series Joined by Time
// The time increment is based on the range you've chosen... so choose wisely
// 1 day = 1 min
// Probably should run the clean and repair first, though for a check it should do
// just dandy
//
// July 4th
// Adding functionality to put in the MBD ratios
// These come from Trigger-> MBD and ZDC
// At the moment, I'm going to ignore the issue of when we weren't running the |Vz|<10 trigger
// The macro assumes both a 1 minute time stamp AND that both files were cleaned and repaired
// It also lets one filter a single day

void PlotFromCSV(const char* filename = "GL1ZDC.csv", const char* filename2 = "Ratio.csv", const char* dayFilter = "") {
    
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    ifstream file(filename);
    if (!file.is_open()) {
        Error("Plot", "Could not open file: %s", filename);
        return;
    }
    
    ifstream file2(filename2);
    if (!file2.is_open()) {
        Error("Plot", "Could not open file: %s", filename2);
        return;
    }
    
    cout<<"Analyzing "<<filename<<" and "<<filename2<<endl;
    TString s(filename);
    s = gSystem->BaseName(s);
    s.ReplaceAll(".csv", "");
    cout<<s<<endl;
    
    vector<double> times;
    vector<double> rhicValues;
    vector<double> zdcValues;
    vector<double> mbdNarrowR;
    vector<double> mbdWideR;
    
    string line;
   // getline(file, line); // skip header --> I kill the header in cleaning because I am an ogre
    
    while (getline(file, line)) {
        stringstream ss(line);
        string timeStr, rhicStr, zdcStr;
        
        getline(ss, timeStr, ',');
        getline(ss, rhicStr, ',');
        getline(ss, zdcStr, ',');
        
        if (strlen(dayFilter) > 0) {
            if (strncmp(timeStr.c_str(), dayFilter, strlen(dayFilter)) != 0)
                continue; // skip this line if it doesn't match the date
        }
        
        int year, month, day, hour, min, sec;
        sscanf(timeStr.c_str(), "%d-%d-%d %d:%d:%d", &year, &month, &day, &hour, &min, &sec);
        TDatime t(year, month, day, hour, min, sec);
        times.push_back(t.Convert());
        
        
        
        rhicValues.push_back(atof(rhicStr.c_str()));
        zdcValues.push_back(atof(zdcStr.c_str()));
    }
    
    while (getline(file2, line)) {
        //sloppy code reuse, but essentially it's the same thing
        //don't do the time thing again, it makes the time lords unhappy
        stringstream ss(line);
        string timeStr, rhicStr, zdcStr;
        
        getline(ss, timeStr, ',');
        getline(ss, rhicStr, ',');
        getline(ss, zdcStr, ',');
        
        if (strlen(dayFilter) > 0) {
            if (strncmp(timeStr.c_str(), dayFilter, strlen(dayFilter)) != 0)
                continue; // skip this line if it doesn't match the date
        }
        
        double mb1 =atof(rhicStr.c_str());
        double mb2 =atof(zdcStr.c_str());
        if (std::isnan(mb1))
            mb1 = 0.0; //I set all nans to zero
        if (std::isnan(mb2))
            mb2 = 0.0; //I set all nans to zero
        
        mbdNarrowR.push_back(mb1);
        mbdWideR.push_back(mb2);
    }
    
    int n = times.size();
    
    // Create histograms with variable bin edges -> Just in case things get weird
    vector<double> binEdges;
    for (int i = 0; i < n; ++i) {
        binEdges.push_back(times[i]);
    }
    // Add one more bin edge to cover the last bin
    binEdges.push_back(times.back() + (times.back() - times[n-2]));

    //hRHIC = RHIC is running and ZDC rate
    TH1D* hRhic = new TH1D(Form("hRhic_%s",s.Data()), "RHIC Scalar vs Time", n, &binEdges[0]);
    //hZDC = sPHENIX is running and ZDC rate
    TH1D* hZdc  = new TH1D(Form("hZdc_%s",s.Data()),  "ZDC NS vs Time",     n, &binEdges[0]);
    TH1D* hMbdNarrowRate = (TH1D*)hZdc->Clone(Form("hMbdNarrowRate_%s", s.Data()));
    TH1D* hMbdWideRate = (TH1D*)hZdc->Clone(Form("hMbdWideRate_%s", s.Data()));
    
    //Ok, I've already built in the assumptions that everyone had 1 minute interval indentical axes...
    TH1D* hRhicCumulative = (TH1D*)hRhic->Clone(Form("hRhicCumulative_%s", s.Data()));
    TH1D* hZdcCumulative = (TH1D*)hZdc->Clone(Form("hZdcCumulative_%s", s.Data()));
    TH1D* hMbdNarrowCumulative = (TH1D*)hZdc->Clone(Form("hMbdNarrowCumulative_%s", s.Data()));
    TH1D* hMbdWideCumulative = (TH1D*)hZdc->Clone(Form("hMbdWideCumulative_%s", s.Data()));
    
    
    for (int i = 0; i < n; ++i) {
        hRhic->SetBinContent(i + 1, rhicValues[i]);
        hZdc->SetBinContent(i + 1, zdcValues[i]);
        hMbdNarrowRate->SetBinContent(i+1,zdcValues[i]*mbdNarrowR[i]); //turns the ratio into a rate
        hMbdWideRate->SetBinContent(i+1,zdcValues[i]*mbdWideR[i]);
        //if (rhicValues[i] > 0)
          //  cout<<"r "<<rhicValues[i]<<" sp "<<zdcValues[i]<<" mbn "<<rhicValues[i]*mbdNarrowR[i]<<" mbw "<<rhicValues[i]*mbdWideR[i]<<" ratios "<<mbdNarrowR[i]<<", "<<mbdWideR[i]<<endl;
    }

    //Theoretically I could do this in the above loop, but nah.
    double totalRHIC = 0; double totalsPHENIX = 0; double totalNarrow = 0; double totalWide = 0;
    double dtR   = hRhic->GetBinWidth(1);        // in seconds (time per bin)
    double dtZ   = hZdc->GetBinWidth(1);        // in seconds (time per bin)
    
    double total_time = 0;
    double nonzero_timeR = 0;
    double nonzero_timeZ = 0;
    
    for (int i = 1; i <= hRhic->GetNbinsX(); ++i) {
        double rateR = hRhic->GetBinContent(i);      // in kHz (events/sec)
        if (rateR < 200) //there are sometimes spikes
            totalRHIC += rateR * dtR*1000;
        double rateZ = hZdc->GetBinContent(i);      // in kHz (events/sec)
        totalsPHENIX += rateZ * dtZ*1000;
        totalNarrow += hMbdNarrowRate->GetBinContent(i)*1000*dtZ;
        totalWide += hMbdWideRate->GetBinContent(i)*1000*dtZ;
        total_time += dtR;
        hRhicCumulative->SetBinContent(i, totalRHIC);
        hZdcCumulative->SetBinContent(i, totalsPHENIX);
        hMbdNarrowCumulative->SetBinContent(i,totalNarrow);
        hMbdWideCumulative->SetBinContent(i,totalWide);
            
        if (rateR > 0)
                nonzero_timeR += dtR;
        if (rateZ > 0)
                nonzero_timeZ += dtR;
    }
    
    cout<<times[0]<<endl;
    time_t t = static_cast<time_t>(times[0]);
    struct tm *tm_info = localtime(&t);  // or use gmtime(&t) for UTC
    char buffer[80];
    strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", tm_info);

    cout << buffer << endl;
    cout<<"Total number of ZDC "<<totalRHIC<<" and total number of sPHENIX ZDC "<<totalsPHENIX<<endl;
    cout<<"Out of "<<total_time/3600<<" hours RHIC was on "<<nonzero_timeR/3600<<" hours and sPHENIX "<<nonzero_timeZ/3600<<" hours"<<endl;
    cout<<"Wall Time "<<nonzero_timeR/total_time<<" and sPHENIX Uptime "<<nonzero_timeZ/nonzero_timeR<<endl;

    // Style
    hRhic->SetLineColor(kBlack);
    hRhic->SetLineWidth(3);

    hZdc->SetLineColor(kBlue+2);
    hZdc->SetLineWidth(3);
    hZdc->SetFillColorAlpha(kBlue - 9,0.3);
    hZdc->SetFillStyle(1001); // solid
    
    hMbdNarrowRate->SetLineWidth(3);
    hMbdNarrowRate->SetLineColor(kRed);
    hMbdNarrowRate->SetFillStyle(1001);
    hMbdNarrowRate->SetFillColorAlpha(kRed,0.1);
    
    hMbdWideRate->SetLineWidth(3);
    hMbdWideRate->SetLineColor(kMagenta);
    hMbdWideRate->SetFillStyle(1001);
    hMbdWideRate->SetFillColorAlpha(kMagenta,0.1);

    // Canvas
    TCanvas* c = new TCanvas("c", "RHIC Scalar and ZDC NS vs Time", 1200, 600);
    c->SetLeftMargin(0.1);
    c->SetBottomMargin(0.15);
    gPad->SetLogy();
    hRhic->GetYaxis()->SetRangeUser(1, 1000);
    hRhic->GetXaxis()->SetTimeDisplay(1);
//    hRhic->GetXaxis()->SetTimeFormat("%m-%d %H:%M");
    hRhic->GetXaxis()->SetTimeFormat("%m/%d\n-%H:%M");
    hRhic->GetXaxis()->SetTimeOffset(0, "local");
    hRhic->GetXaxis()->SetNdivisions(508, kTRUE);
    hRhic->GetXaxis()->SetTitle("Time");
    hRhic->GetYaxis()->SetTitle("Rate (kHz)");
    
    hRhic->GetXaxis()->SetTitleSize(0.05);
    hRhic->GetYaxis()->SetTitleSize(0.05);
    hRhic->GetXaxis()->SetTitleOffset(1.2);
    hRhic->GetXaxis()->SetLabelSize(0.045);
    hRhic->GetYaxis()->SetLabelSize(0.045);
    
    hRhic->Draw("HIST");
    hZdc->Draw("HIST SAME");
    hMbdNarrowRate->Draw("HIST SAME");
    hMbdWideRate->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.3, 0.65, 0.85, 0.85);
    leg->SetFillStyle(1001); // solid background
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.05);
//    leg->AddEntry(hRhic, "RHIC", "l");
    leg->AddEntry(hRhic,Form("RHIC - Wall Time %.0f Percent",100*nonzero_timeR/total_time),"l");
    //leg->AddEntry(hZdc, "sPHENIX", "lf");
    leg->AddEntry(hZdc,Form("sPHENIX - Uptime Time %.0f Percent",100*nonzero_timeZ/nonzero_timeR),"l");
    leg->AddEntry(hMbdNarrowRate, "MBD Narrow", "l");
    leg->AddEntry(hMbdWideRate, "MBD Wide", "l");

    leg->Draw();
    
    c->SaveAs(Form("c%s.png",s.Data()));
    
    TCanvas* c2 = new TCanvas("c2", "Cumulative Event Count", 1200, 600);
    c2->SetLeftMargin(0.1);
    c2->SetBottomMargin(0.15);
    hRhicCumulative->GetYaxis()->SetTitle("Cumulative Events");
    hRhicCumulative->GetXaxis()->SetTitle("Time");
    hRhicCumulative->GetXaxis()->SetTimeDisplay(1);
    hRhicCumulative->GetXaxis()->SetTimeFormat("%m/%d\n-%H:%M");
    hRhicCumulative->GetXaxis()->SetTimeOffset(0, "local");
    hRhicCumulative->GetXaxis()->SetNdivisions(508, kTRUE);
    hRhicCumulative->GetXaxis()->SetTitleSize(0.05);
    hRhicCumulative->GetYaxis()->SetTitleSize(0.05);
    hRhicCumulative->GetXaxis()->SetTitleOffset(1.2);
    hRhicCumulative->GetXaxis()->SetLabelSize(0.045);
    hRhicCumulative->GetYaxis()->SetLabelSize(0.045);

    hRhicCumulative->SetLineColor(kBlack);
    hRhicCumulative->SetLineWidth(3);

    hZdcCumulative->SetLineColor(kBlue+2);
    hZdcCumulative->SetLineWidth(3);
    hZdcCumulative->SetFillColorAlpha(kBlue - 9, 0.4);
    hZdcCumulative->SetFillStyle(1001);
    
    hMbdNarrowCumulative->SetLineWidth(3);
    hMbdNarrowCumulative->SetLineColor(kRed);
    hMbdNarrowCumulative->SetFillStyle(1001);
    hMbdNarrowCumulative->SetFillColorAlpha(kRed,0.1);
    
    hMbdWideCumulative->SetLineWidth(3);
    hMbdWideCumulative->SetLineColor(kMagenta);
    hMbdWideCumulative->SetFillStyle(1001);
    hMbdWideCumulative->SetFillColorAlpha(kMagenta,0.1);

    hRhicCumulative->SetMaximum(hRhicCumulative->GetMaximum()*1.5);
    hRhicCumulative->Draw("HIST");
    hZdcCumulative->Draw("HIST SAME");
    hMbdWideCumulative->Draw("HIST SAME");
    hMbdNarrowCumulative->Draw("HIST SAME");

    TLegend* leg2 = new TLegend(0.15, 0.6, 0.4, 0.88);
    leg2->SetFillStyle(1001);
    leg2->SetFillColor(kWhite);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.05);
    leg2->AddEntry(hRhicCumulative, Form("Total ZDCNS Events %.1e",totalRHIC), "l");
    leg2->AddEntry(hZdcCumulative, Form("sPHENIX on, ZDCNS Events %.1e",totalsPHENIX), "lf");
    leg2->AddEntry(hMbdNarrowCumulative, Form("MBDNS>2 |z|< 10 cm %.1e",totalNarrow), "lf");
    leg2->AddEntry(hMbdWideCumulative, Form("MBDNS>2 |z|< 150 cm %.1e",totalWide), "lf");
    leg2->Draw();

    TPad *imagePad = new TPad("imagePad", "Image Pad", 0.15, 0.35, 0.4, 0.55);
    imagePad->SetFillColor(0);
    imagePad->SetBorderMode(0);
    imagePad->SetFrameBorderMode(0);
    imagePad->Draw();
    imagePad->cd();
    
    TImage *img = TImage::Open("sphenix-logo-white-bg.png");  // must be in working directory
    if (img) {
        cout<<"Drawing logo"<<endl;
        img->Draw(); // in NDC coordinates (xlow, ylow, xhigh, yhigh)
    }
    c2->SaveAs(Form("cumulative_%s.png", s.Data()));
    
    TFile *_file0 = TFile::Open(Form("r%s.root",s.Data()),"RECREATE");
    
    hRhic->Write();
    hZdc->Write();
    hMbdNarrowRate->Write();
    hMbdWideRate->Write();
    
    hRhicCumulative->Write();
    hZdcCumulative->Write();
    hMbdNarrowCumulative->Write();
    hMbdWideCumulative->Write();
    
    
}
