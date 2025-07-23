#include "sepd_helper.h"

void take_runnumber(int);

void do_plots()
{
  //take_runnumber(54912);

  int run;
  ifstream fin;
  fin.open("runlist.txt");
  while ( fin >> run )
    {
      take_runnumber(run);
    }
  fin.close();

}

void take_runnumber(int runnumber)
{

  //int runnumber = 54912;
  char filename[200];
  sprintf(filename,"/sphenix/lustre01/sphnxpro/commissioning/online_monitoring/histograms/SEPDMON/Run_%05d-SEPDMON_0.root",runnumber);
  TFile* file = TFile::Open(filename);

  if ( file == NULL )
    {
      cout << "Missing file " << filename << ", skipping" << endl;
      return;
    }

  TH1F* h_event = (TH1F*)file->Get("h_event");
  TH1F* h_hits_all_channel = (TH1F*)file->Get("h_hits_all_channel");

  TCanvas* c1 = new TCanvas("c1","",1000,500);

  TPad* Pad[15];
  Pad[0] = new TPad("sepdpad0", "Left", 0., 0., 0.5, 1);
  Pad[1] = new TPad("sepdpad1", "Right", 0.5, 0., 1, 1);
  Pad[0]->Draw();
  Pad[1]->Draw();

  // ----------------------------
  // --- begin Rosi (mostly ) ---
  // ----------------------------

  TH2* polar_histS = new TH2F("polar_histS","polar_hist",
                               24, 0, 2*M_PI,
                               16, 0.15, 3.5);
  TH2* polar_histN = new TH2F("polar_histN","polar_hist",
                               24, 0, 2*M_PI,
                               16, 0.15, 3.5);

  //tile 0 is 2x the angular size of the rest of the tiles and needs a separate histogram
  TH2* polar_histS01 = new TH2F("polar_histS01","polar_hist",
                                 12, 0, 2*M_PI,
                                 16, 0.15, 3.5);
  TH2* polar_histN01 = new TH2F("polar_histN01","polar_hist",
                                 12, 0, 2*M_PI,
                                 16, 0.15, 3.5);

  // --- normalize
  //h_ADC_all_channel->Divide(h_hits_all_channel);
  int nevt = h_event->GetEntries();
  //h_ADC_all_channel->Scale(1.0/nevt);
  h_hits_all_channel->Scale(1.0/nevt);
  for ( int i = 0; i < 768; ++i )
    {
      int adc_channel = i;
      float adc_signal = h_hits_all_channel->GetBinContent(i+1);
      if ( adc_signal <= 0.0001 ) adc_signal = 0.0001;
      int tile = returnTile(i);
      int odd = (tile+1)%2;
      //int ring = returnRing(adc_channel);
      int sector = returnSector(adc_channel);
      int arm = returnArm(adc_channel);
      if ( arm == 0 )
        {
          if ( tile == 0 ) polar_histS01->SetBinContent(sector+1,1,adc_signal);
          else polar_histS->SetBinContent(sector*2+1+odd,(tile+1)/2+1,adc_signal);
        }
      if ( arm == 1 )
        {
          if ( tile == 0 ) polar_histN01->SetBinContent(sector+1,1,adc_signal);
          else polar_histN->SetBinContent(sector*2+1+odd,(tile+1)/2+1,adc_signal);
        }
    }

  // -------------------------
  // --- end Rosi (mostly) ---
  // -------------------------

  // --- may need to update these depending on whether there are "hot" tiles
  double zmin = 0.0;
  //double zmax = 0.1;
  double zmax = 1.0;
  //double zmax = 300;
  //double zmax = 1.1*h_ADC_all_channel->GetMaximum();

  TText trun;
  trun.SetNDC();
  trun.SetTextFont(42);
  trun.SetTextSize(0.05);

  trun.DrawText(0.1,0.9,Form("%05d",runnumber));

  TText tarm;
  tarm.SetNDC();
  tarm.SetTextFont(42);
  tarm.SetTextSize(0.05);

  gStyle->SetOptStat(0);
  // ---
  Pad[0]->cd();
  polar_histS->GetZaxis()->SetRangeUser(zmin,zmax);
  polar_histS01->GetZaxis()->SetRangeUser(zmin,zmax);
  // gPad->SetLeftMargin(0.2);
  // gPad->SetRightMargin(0.0);
  gPad->SetTicks(1,1);
  gPad->DrawFrame(-3.8, -3.8,3.8, 3.8);
  polar_histS->Draw("same col pol AH");
  polar_histS01->Draw("same col pol AH");
  tarm.DrawText(0.45,0.91,"South");
  gStyle->SetPalette(57);
  // ---
  Pad[1]->cd();
  polar_histN->GetZaxis()->SetRangeUser(zmin,zmax);
  polar_histN01->GetZaxis()->SetRangeUser(zmin,zmax);
  gPad->SetLeftMargin(0.05);
  gPad->SetRightMargin(0.15);
  gPad->SetTicks(1,1);
  gPad->DrawFrame(-3.8, -3.8,3.8, 3.8);
  polar_histN->Draw("same colz pol AH");
  polar_histN01->Draw("same col pol AH");
  tarm.DrawText(0.40,0.91,"North");
  //tarm.DrawText(0.35,0.91,"North");
  //tarm.DrawText(0.45,0.91,"North");
  gStyle->SetPalette(57);

  c1->Print(Form("AllPlots/sepd_wheel_%05d.png",runnumber));

  // -------------------------------------------------------

  if ( c1 ) delete c1;

  TCanvas* c2 = new TCanvas("c2","",1000,500);
  trun.DrawText(0.1,0.9,Form("%05d",runnumber));

  TH2* h_ADC_corr = (TH2*)file->Get("h_ADC_corr");
  TH2* h_hits_corr = (TH2*)file->Get( "h_hits_corr");

  Pad[4] = new TPad("sepdpad4", "Left", 0., 0., 0.5, 1);
  Pad[5] = new TPad("sepdpad5", "Right", 0.5, 0., 1, 1);
  Pad[4]->Draw();
  Pad[5]->Draw();

  Pad[4]->cd();
  h_ADC_corr->GetYaxis()->SetNdivisions(505);
  h_ADC_corr->GetXaxis()->SetNdivisions(505);
  h_ADC_corr->GetYaxis()->SetRangeUser(0,1.5e6);
  h_ADC_corr->GetXaxis()->SetRangeUser(0,1.5e6);
  h_ADC_corr->Draw("COLZ");
  // ---
  gPad->SetLogz();
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.2);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(57);
  gPad->SetTicky();
  gPad->SetTickx();
  // ---
  Pad[5]->cd();
  h_hits_corr->GetYaxis()->SetNdivisions(505);
  h_hits_corr->GetXaxis()->SetNdivisions(505);
  h_hits_corr->GetYaxis()->SetRangeUser(0,380);
  h_hits_corr->GetXaxis()->SetRangeUser(0,380);
  h_hits_corr->Draw("COLZ");
  // ---
  gPad->SetLogz();
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.2);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(57);
  gPad->SetTicky();
  gPad->SetTickx();
  // ---

  trun.DrawText(0.1,0.9,Form("%05d",runnumber));
  c2->Print(Form("AllPlots/sepd_north_south_%05d.png",runnumber));

  // -------------------------------------------------------------

  if ( c2 ) delete c2;

  TCanvas* c3 = new TCanvas("c3","",500,800);
  trun.DrawText(0.1,0.9,Form("%05d",runnumber));

  TH1 *h_waveform_time = (TH1*)file->Get("h1_waveform_time");
  TH1 *h_waveform_pedestal = (TH1*)file->Get("h1_waveform_pedestal");
  TH2 *h2_sepd_waveform = (TH2*)file->Get("h2_sepd_waveform");

  //TH1 *h_event = cl->getHisto("SEPDMON_0", "h_event");
  //int nevt = h_event->GetEntries();

  Pad[6] = new TPad("sepdpad6", "ADC vs sample #", 0.0, 0.6, 1.0, 0.95, 0);
  Pad[7] = new TPad("sepdpad7", "counts vs sample #", 0.0, 0.3, 1.0, 0.6, 0);
  Pad[8] = new TPad("sepdpad8", "pedestals", 0.0, 0.0, 1.0, 0.3, 0);
  Pad[6]->Draw();
  Pad[7]->Draw();
  Pad[8]->Draw();

  // Pad[6]->cd();
  // if (!h2_sepd_waveform || !h_waveform_time || !h_waveform_pedestal)
  // {
  //   //cout which one is not found
  //   if (!h2_sepd_waveform) std::cout << "h2_sepd_waveform not found" << std::endl;
  //   if (!h_waveform_time) std::cout << "h_waveform_time not found" << std::endl;
  //   if (!h_waveform_pedestal) std::cout << "h_waveform_pedestal not found" << std::endl;
  //   DrawDeadServer(transparent[canvasindex]);
  //   TC[canvasindex]->SetEditable(0);
  //   return -1;
  // }

  Pad[6]->cd();
  gStyle->SetTitleFontSize(0.03);
  float ymaxp = h2_sepd_waveform->ProfileX()->GetMaximum();
  float ymaxdraw = ymaxp * 10; // was originally 20, but that is too much
  h2_sepd_waveform->GetYaxis()->SetRangeUser(0,ymaxdraw);
  h2_sepd_waveform->GetXaxis()->SetRangeUser(0, 11);
  h2_sepd_waveform->Draw("colz");
  // --- add a profile on top
  TProfile* tp1f_sepd_waveform = h2_sepd_waveform->ProfileX();
  tp1f_sepd_waveform->SetLineColor(kBlack);
  tp1f_sepd_waveform->SetLineWidth(2);
  tp1f_sepd_waveform->Draw("same");
  // --- draw vertical lines where the waveform should be
  // x1 y1 x2 y2
  TLine* lineleft = new TLine(4.5,0,4.5,ymaxdraw);
  TLine* lineright = new TLine(7.5,0,7.5,ymaxdraw);
  lineleft->SetLineColor(kBlack);
  lineleft->SetLineWidth(2);
  lineleft->SetLineStyle(2);
  lineleft->Draw();
  lineright->SetLineColor(kBlack);
  lineright->SetLineWidth(2);
  lineright->SetLineStyle(2);
  lineright->Draw();

  float tsize = 0.09;
  //h2_sepd_waveform->GetXaxis()->SetNdivisions(510, kTRUE);
  h2_sepd_waveform->GetXaxis()->SetNdivisions(12);
  h2_sepd_waveform->GetXaxis()->SetTitle("Sample #");
  h2_sepd_waveform->GetYaxis()->SetTitle("Waveform [ADC]");
  h2_sepd_waveform->GetXaxis()->SetLabelSize(tsize/1.15);
  h2_sepd_waveform->GetYaxis()->SetLabelSize(tsize/1.15);
  h2_sepd_waveform->GetXaxis()->SetTitleSize(tsize/1.15);
  h2_sepd_waveform->GetYaxis()->SetTitleSize(tsize/1.15);
  h2_sepd_waveform->GetXaxis()->SetTitleOffset(1.0);
  h2_sepd_waveform->GetYaxis()->SetTitleOffset(1.3);
  gPad->SetLogz();
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.2);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(57);
  gPad->SetTicky();
  gPad->SetTickx();


  Pad[7]->cd();

  gStyle->SetTitleFontSize(0.06);

  h_waveform_time->GetXaxis()->SetRangeUser(0,11);
  //h_waveform_time->Scale(1.0/nevt);
  h_waveform_time->Draw("hist");
  // ---
  //h_waveform_time->GetXaxis()->SetNdivisions(510, kTRUE);
  h_waveform_time->GetXaxis()->SetNdivisions(12);
  h_waveform_time->GetXaxis()->SetTitle("Sample #");
  h_waveform_time->GetYaxis()->SetTitle("Normalized Counts");
  h_waveform_time->GetXaxis()->SetLabelSize(tsize);
  h_waveform_time->GetYaxis()->SetLabelSize(tsize);
  h_waveform_time->GetXaxis()->SetTitleSize(tsize);
  h_waveform_time->GetYaxis()->SetTitleSize(tsize);
  h_waveform_time->GetXaxis()->SetTitleOffset(1.0);
  h_waveform_time->GetYaxis()->SetTitleOffset(1.25);
  h_waveform_time->SetFillColorAlpha(kBlue, 0.1);
  if ( h_waveform_time->GetEntries() )
    {
      h_waveform_time->Scale(1.0/h_waveform_time->GetEntries());
    }
  gPad->Update();
  // draw two black lines for the okay timing range
  TLine line3(4.5, 0, 4.5, gPad->GetFrame()->GetY2());
  line3.SetLineColor(1);
  line3.SetLineWidth(3);
  line3.SetLineStyle(1);
  line3.DrawLine(4.5, 0, 4.5, gPad->GetFrame()->GetY2());
  // high line
  TLine line4(7.5, 0, 7.5, gPad->GetFrame()->GetY2());
  line4.SetLineColor(1);
  line4.SetLineWidth(3);
  line4.SetLineStyle(1);
  line4.DrawLine(7.5, 0, 7.5, gPad->GetFrame()->GetY2());
  // Draw a red line at mean x
  TLine line5(h_waveform_time->GetMean(), 0, h_waveform_time->GetMean(), gPad->GetFrame()->GetY2());
  line5.SetLineColor(2);
  line5.SetLineWidth(3);
  line5.SetLineStyle(1);
  line5.DrawLine(h_waveform_time->GetMean(), 0, h_waveform_time->GetMean(), gPad->GetFrame()->GetY2());
  // ---
  gPad->SetTopMargin(0.06);
  gPad->SetBottomMargin(0.18);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.2);
  gStyle->SetOptStat(0);
  gPad->SetTicky();
  gPad->SetTickx();

  Pad[8]->cd();

  gStyle->SetTitleFontSize(0.06);

  // x-axis range is set in SepdMon.cc, need to change there if want a wider range
  //h_waveform_pedestal->Scale(1.0/nevt);
  h_waveform_pedestal->Draw("hist");
  h_waveform_pedestal->GetXaxis()->SetNdivisions(505);
  h_waveform_pedestal->GetXaxis()->SetTitle("ADC Pedestal");
  h_waveform_pedestal->GetYaxis()->SetTitle("Normalized Counts");
  h_waveform_pedestal->GetXaxis()->SetLabelSize(tsize);
  h_waveform_pedestal->GetYaxis()->SetLabelSize(tsize);
  h_waveform_pedestal->GetXaxis()->SetTitleSize(tsize);
  h_waveform_pedestal->GetYaxis()->SetTitleSize(tsize);
  h_waveform_pedestal->GetXaxis()->SetTitleOffset(0.9);
  h_waveform_pedestal->GetYaxis()->SetTitleOffset(1.25);
  h_waveform_pedestal->SetFillColorAlpha(kBlue, 0.1);
  if ( h_waveform_pedestal->GetEntries() )
    {
      h_waveform_pedestal->Scale(1.0/h_waveform_pedestal->GetEntries());
    }
  gPad->Update();
  TLine line6(1000, 0, 1000, gPad->GetFrame()->GetY2());
  line6.SetLineColor(1);
  line6.SetLineWidth(3);
  line6.SetLineStyle(1);
  line6.DrawLine(1000, 0, 1000, gPad->GetFrame()->GetY2());
  TLine line7(2000, 0, 2000, gPad->GetFrame()->GetY2());
  line7.SetLineColor(1);
  line7.SetLineWidth(3);
  line7.SetLineStyle(1);
  line7.DrawLine(2000, 0, 2000, gPad->GetFrame()->GetY2());
  gPad->SetTopMargin(0.06);
  gPad->SetBottomMargin(0.18);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.2);
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  gPad->SetTicky();
  gPad->SetTickx();

  trun.DrawText(0.1,0.9,Form("%05d",runnumber));
  c3->Print(Form("AllPlots/sepd_waveform_%05d.png",runnumber));

  // ----------------------------------------------------------

  if ( c3 ) delete c3;

  TCanvas* c4 = new TCanvas("c4","",800,800);
  trun.DrawText(0.1,0.9,Form("%05d",runnumber));

  Pad[10] = new TPad("sepdpad10", "packet event check", 0.0, 0.6, 1.0 / 2, 0.95, 0);
  Pad[11] = new TPad("sepdpad11", "packet size", 0.0, 0.3, 1.0 / 2, 0.6, 0);
  Pad[12] = new TPad("sepdpad12", "packet channels", 0.0, 0.0, 1.0 / 2, 0.3, 0);
  Pad[13] = new TPad("sepdpad13", "event number offset", 0.5, 0.6, 1.0, 0.95, 0);

  Pad[10]->Draw();
  Pad[11]->Draw();
  Pad[12]->Draw();
  Pad[13]->Draw();

  TH1 *h1_packet_number = (TH1*)file->Get("h1_packet_number");
  TH1 *h1_packet_length = (TH1*)file->Get("h1_packet_length");
  TH1 *h1_packet_chans =  (TH1*)file->Get("h1_packet_chans");
  TH1 *h1_packet_event =  (TH1*)file->Get("h1_packet_event");

  if (!h1_packet_number || !h1_packet_length || !h1_packet_chans || !h1_packet_event)
  {
    // print out which is not found
    if (!h1_packet_number) std::cout << "h1_packet_number not found" << std::endl;
    if (!h1_packet_length) std::cout << "h1_packet_length not found" << std::endl;
    if (!h1_packet_chans) std::cout << "h1_packet_chans not found" << std::endl;
    if (!h1_packet_event) std::cout << "h1_packet_event not found" << std::endl;
  }

  // int maxbin = h1_packet_event->GetMaximumBin();
  int maxy = h1_packet_event->GetMaximum();
  // substract all none zero bin by maxy
  for (int i = 1; i <= h1_packet_event->GetNbinsX(); i++)
  {
    if (h1_packet_event->GetBinContent(i) != 0)
    {
      h1_packet_event->SetBinContent(i, h1_packet_event->GetBinContent(i) - maxy);
    }
  }

  // find the x range for h1_packet_number
  double xmin = h1_packet_number->GetXaxis()->GetXmin();
  double xmax = h1_packet_number->GetXaxis()->GetXmax();

  TLine *one = new TLine(xmin, 1, xmax, 1);
  one->SetLineStyle(7);

  // --- Martin says 1047 for NZS
  int PACKET_SIZE = 1047;
  TLine *goodSize = new TLine(xmin, PACKET_SIZE, xmax, PACKET_SIZE);
  goodSize->SetLineStyle(7);

  // --- 128 channels per packet
  int N_CHANNELS = 128;
  TLine *goodChans = new TLine(xmin, N_CHANNELS, xmax, N_CHANNELS);
  goodChans->SetLineStyle(7);

  float param = 0.95;

  //float tsize = 0.08;
  tsize = 0.08;
  TLegend *leg = new TLegend(0.3, 0.16, 0.95, 0.4);
  leg->SetTextFont(42);
  leg->SetTextSize(tsize);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  TLine *warnLineOne = new TLine(xmin, param * 1, xmax, param * 1);
  warnLineOne->SetLineStyle(7);
  warnLineOne->SetLineColor(2);
  leg->AddEntry(warnLineOne, "95% Threshold", "l");

  // --- packet size is 1047 NZS
  TLine *warnLineSize = new TLine(xmin, param * PACKET_SIZE, xmax, param * PACKET_SIZE);
  warnLineSize->SetLineStyle(7);
  warnLineSize->SetLineColor(2);

  // --- 128 channels per packet
  TLine *warnLineChans = new TLine(xmin, param * N_CHANNELS, xmax, param * N_CHANNELS);
  warnLineChans->SetLineStyle(7);
  warnLineChans->SetLineColor(2);


  // --- this one is okay
  Pad[10]->cd();
  //float tsize = 0.08;
  h1_packet_number->GetYaxis()->SetRangeUser(0.0, 1.3);
  h1_packet_number->Draw("hist");
  one->Draw("same");
  warnLineOne->Draw("same");
  leg->Draw("same");
  h1_packet_number->GetXaxis()->SetNdivisions(6);
  h1_packet_number->GetXaxis()->SetTitle("Packet #");
  h1_packet_number->GetYaxis()->SetTitle("% Of Events Present");
  // the sizing is funny on this pad...
  h1_packet_number->GetXaxis()->SetLabelSize(tsize/1.15);
  h1_packet_number->GetYaxis()->SetLabelSize(tsize/1.15);
  h1_packet_number->GetXaxis()->SetTitleSize(tsize/1.15);
  h1_packet_number->GetYaxis()->SetTitleSize(tsize/1.15);
  h1_packet_number->GetXaxis()->SetTitleOffset(1);
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gStyle->SetOptStat(0);
  gPad->SetTicky();
  gPad->SetTickx();

  // --- this one is okay (1047 for NZS, variable for ZS)
  Pad[11]->cd();
  h1_packet_length->Draw("hist");
  h1_packet_length->GetYaxis()->SetRangeUser(0, 1200);
  goodSize->Draw("same");
  warnLineSize->Draw("same");
  h1_packet_length->GetXaxis()->SetNdivisions(6);
  h1_packet_length->GetXaxis()->SetTitle("Packet #");
  h1_packet_length->GetYaxis()->SetTitle("Average Packet Size");
  h1_packet_length->GetXaxis()->SetLabelSize(tsize);
  h1_packet_length->GetYaxis()->SetLabelSize(tsize);
  h1_packet_length->GetXaxis()->SetTitleSize(tsize);
  h1_packet_length->GetYaxis()->SetTitleSize(tsize);
  h1_packet_length->GetXaxis()->SetTitleOffset(1);
  h1_packet_length->GetYaxis()->SetTitleOffset(0.8);
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gStyle->SetOptStat(0);
  gPad->SetTicky();
  gPad->SetTickx();

  // --- this one is okay
  Pad[12]->cd();
  h1_packet_chans->Draw("hist");
  h1_packet_chans->GetYaxis()->SetRangeUser(0, 150);
  goodChans->Draw("same");
  warnLineChans->Draw("same");
  h1_packet_chans->GetXaxis()->SetNdivisions(6);
  h1_packet_chans->GetXaxis()->SetTitle("Packet #");
  h1_packet_chans->GetYaxis()->SetTitle("Average # of Channels");
  h1_packet_chans->GetXaxis()->SetLabelSize(tsize);
  h1_packet_chans->GetYaxis()->SetLabelSize(tsize);
  h1_packet_chans->GetXaxis()->SetTitleSize(tsize);
  h1_packet_chans->GetYaxis()->SetTitleSize(tsize);
  h1_packet_chans->GetXaxis()->SetTitleOffset(0.8);
  h1_packet_chans->GetYaxis()->SetTitleOffset(0.8);
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gStyle->SetOptStat(0);
  gPad->SetTicky();
  gPad->SetTickx();

  Pad[13]->cd();
  //  h1_packet_event->Draw("hist");
  // h1_packet_event->SetLineColor(kWhite);;
  // h1_packet_event->Draw("AH");
  double ymax = h1_packet_event->GetMaximum();
  double ymin = h1_packet_event->GetMinimum();

  // --- this one seems okay
  h1_packet_event->GetXaxis()->SetNdivisions(6);
  h1_packet_event->GetYaxis()->SetRangeUser(ymin - 0.3 * (ymax - ymin + 30), ymax + 0.3 * (ymax - ymin + 30));
  // h1_packet_event->GetXaxis()->SetTitle("Packet #");
  // h1_packet_event->GetYaxis()->SetTitle("clock offset");
  h1_packet_event->GetXaxis()->SetLabelSize(tsize/1.2);
  h1_packet_event->GetYaxis()->SetLabelSize(tsize/1.2);
  h1_packet_event->GetXaxis()->SetTitleSize(tsize/1.2);
  h1_packet_event->GetYaxis()->SetTitleSize(tsize/1.2);
  h1_packet_event->GetXaxis()->SetTitleOffset(0.8);
  h1_packet_event->GetYaxis()->SetTitleOffset(1.2);
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.2);
  gStyle->SetOptStat(0);
  gPad->SetTicky();
  gPad->SetTickx();




  std::vector<int> badPackets;
  std::string list_of_bad_packets;
  bool packet_is_bad[7];
  for (int i = 1; i <= 6; i++)
  {
    packet_is_bad[i] = false;
    bool missing = false;
    bool badnumber = false;
    //bool badlength = false;
    bool badchans = false;
    if (h1_packet_number->GetBinContent(i) == 0)
    {
      missing = true;
    }
    if (h1_packet_number->GetBinContent(i) < param)
    {
      badnumber = true;
    }
    // if (h1_packet_length->GetBinContent(i) < param * PACKET_SIZE)
    // {
    //   badlength = true;
    // }
    if (h1_packet_chans->GetBinContent(i) < param * N_CHANNELS)
    {
      badchans = true;
    }
    //if (badnumber || badlength || badchans || missing)
    if (badnumber || badchans || missing)
    {
      int the_bad_packet = (int)h1_packet_number->GetBinCenter(i);
      packet_is_bad[i] = true;
      badPackets.push_back(the_bad_packet);
      list_of_bad_packets += std::to_string(the_bad_packet);
      list_of_bad_packets += ", ";
    }
  }
  // remove the final comma and space
  if (! list_of_bad_packets.empty())
    {
      list_of_bad_packets.resize(list_of_bad_packets.size() - 2);
    }
  // --- draw the packet information
  TText PacketWarn;
  PacketWarn.SetTextFont(42);
  //PacketWarn.SetTextSize(0.04);
  PacketWarn.SetTextSize(0.05);
  PacketWarn.SetTextColor(kBlack);
  PacketWarn.SetNDC();
  //PacketWarn.SetTextAlign(23);
  PacketWarn.DrawText(0.01, 0.75, "Packet Status:");
  for (int i = 1; i <= 6; i++)
    {
      if ( packet_is_bad[i] ) PacketWarn.SetTextColor(kRed);
      // PacketWarn.DrawText(0.01, 0.7 - 0.05 * i, Form("%d: %d%% events, %d size, %d channels, %d offset", i+9000,
      //                                               int(100*h1_packet_number->GetBinContent(i)+0.5),
      //                                               (int)h1_packet_length->GetBinContent(i),
      //                                               (int)h1_packet_chans->GetBinContent(i),
      //                                               (int)h1_packet_event->GetBinContent(i)) );
      PacketWarn.DrawText(0.01, 0.7 - 0.05 * i, Form("%d: %d%% events, %d size, %d channels", i+9000,
                                                     int(100*h1_packet_number->GetBinContent(i)+0.5),
                                                     (int)h1_packet_length->GetBinContent(i),
                                                     (int)h1_packet_chans->GetBinContent(i)) );
      PacketWarn.SetTextColor(kBlack);
    }
  if ( badPackets.size() == 0 )
    {
      PacketWarn.DrawText(0.01, 0.30, Form("No bad packets, everything okay"));
    }
  if ( badPackets.size() == 1 )
    {
      PacketWarn.SetTextColor(kRed);
      PacketWarn.DrawText(0.01, 0.30, Form("%d bad packet: %s",(int)badPackets.size(),list_of_bad_packets.c_str()));
    }
  if ( badPackets.size() > 1 )
    {
      PacketWarn.SetTextColor(kRed);
      PacketWarn.DrawText(0.01, 0.30, Form("%d bad packets: %s",(int)badPackets.size(),list_of_bad_packets.c_str()));
    }

  trun.DrawText(0.1,0.9,Form("%05d",runnumber));
  c4->Print(Form("AllPlots/sepd_packet_%05d.png",runnumber));

  if ( c4 ) delete c4;

  file->Close();

}


