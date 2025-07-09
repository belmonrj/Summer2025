//**************************************
// Rosi Reed
// July 6th
// Same as the plot from the CSV macro, but now doesn't do any calculations
// just pulls in the root file

void PlotFromRoot(const char* filename = "rGL1_615-705.root") {
    
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    TString base = filename;
    base.ReplaceAll(".root", "");
    base.ReplaceAll("r", "");
    
    TFile *_file0 = TFile::Open(filename);
    
    TH1D* hRhic = (TH1D*)gROOT->FindObject(Form("hRhic_%s",base.Data()));
    TH1D* hZdc  = (TH1D*)gROOT->FindObject(Form("hZdc_%s",base.Data()));
    TH1D* hMbdNarrowRate = (TH1D*)gROOT->FindObject(Form("hMbdNarrowRate_%s", base.Data()));
    TH1D* hMbdWideRate = (TH1D*)gROOT->FindObject(Form("hMbdWideRate_%s", base.Data()));
    

    TH1D* hRhicCumulative = (TH1D*)gROOT->FindObject(Form("hRhicCumulative_%s", base.Data()));
    TH1D* hZdcCumulative = (TH1D*)gROOT->FindObject(Form("hZdcCumulative_%s", base.Data()));
    TH1D* hMbdNarrowCumulative = (TH1D*)gROOT->FindObject(Form("hMbdNarrowCumulative_%s", base.Data()));
    TH1D* hMbdWideCumulative = (TH1D*)gROOT->FindObject(Form("hMbdWideCumulative_%s", base.Data()));
    
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
    hRhic->GetYaxis()->SetRangeUser(1, 400);
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
    

    double nonzero_timeR = 0; double total_time = 0; double nonzero_timeZ = 0;
    
    for (int i = 0;i<hRhic->GetNbinsX();i++){
        total_time++;
        if (hRhic->GetBinContent(i+1) > 0)
            nonzero_timeR++;
        if (hZdc->GetBinContent(i+1) > 0)
            nonzero_timeZ++;
    }
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
    
    c->SaveAs(Form("c%s.png",base.Data()));
    
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
    hMbdNarrowCumulative->Draw("HIST SAME");
    hMbdWideCumulative->Draw("HIST SAME");

    double totalRHIC = hRhicCumulative->GetBinContent(hRhicCumulative->GetNbinsX());
    double totalsPHENIX = hZdcCumulative->GetBinContent(hZdcCumulative->GetNbinsX());
    double totalNarrow = hMbdNarrowCumulative->GetBinContent(hMbdNarrowCumulative->GetNbinsX());
    double totalWide = hMbdWideCumulative->GetBinContent(hMbdWideCumulative->GetNbinsX());
    
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
    if (img)
        img->Draw();
    c2->SaveAs(Form("cumulative_%s.png", base.Data()));

    
    
}
