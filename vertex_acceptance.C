void doit(double);

void vertex_acceptance()
{
  doit(3);
  doit(5);
  doit(8);
  doit(10);
  doit(15);
  doit(20);
}

void doit(double shift)
{

  TCanvas c1("c1","",800,600);

  TF1 fun1("fun1","gaus",-100,100);
  fun1.SetParameter(0,1.0);
  fun1.SetParameter(1,0.0);
  fun1.SetParameter(2,15.0);
  fun1.SetLineColor(kBlack);

  TF1 fun2("fun2","gaus",-100,100);
  fun2.SetParameter(0,1.0);
  fun2.SetParameter(1,shift);
  fun2.SetParameter(2,15.0);
  fun2.SetLineColor(kRed);

  TH2D hdummy("hdummy","",1,-30.0,30.0,1,0.0,1.5);
  hdummy.GetXaxis()->SetTitle("z-vertex (cm)");
  hdummy.Draw();
  fun1.Draw("same");
  fun2.Draw("same");

  double int1 = fun1.Integral(-10,10);
  double int2 = fun2.Integral(-10,10);

  double frac = int2/int1;

  TLatex txtwidth(-9.0,1.3,Form("Vertex dist width is %.1f cm",fun2.GetParameter(2)));
  txtwidth.Draw("same");
  TLatex txtshift(-9.0,1.2,Form("Vertex shifted by %.1f cm",fun2.GetParameter(1)));
  txtshift.Draw("same");
  TLatex txtfrac(-9.0,1.1,Form("Ratio of shifted/ideal = %.3f",frac));
  txtfrac.Draw("same");

  TBox tbmvtx(-10.0,0.0,10.0,1.0);
  tbmvtx.SetFillColorAlpha(kBlack,0.2);
  tbmvtx.Draw("f,same");

  c1.Print(Form("acc_zv_s%.1f.png",shift));

}
