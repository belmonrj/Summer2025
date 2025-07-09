const double pi = 3.14159265358979323;

void random_angle()
{

  TH1D th1d_phidiff("th1d_phidiff","",100,-pi,pi);
  TH1D th1d_cosdiff("th1d_cosdiff","",100,-1,1);

  for ( int i = 0; i < 1e4; ++i )
    {
      double phi1 = gRandom->Uniform(-pi/2,pi/2);
      double phi2 = gRandom->Uniform(-pi/2,pi/2);
      th1d_phidiff.Fill(phi1-phi2);
      th1d_cosdiff.Fill(cos(2*(phi1-phi2)));
    }

  TCanvas c1("c1","",800,600);
  th1d_phidiff.Draw();
  th1d_phidiff.SetLineWidth(2);
  c1.Print("phidiff.png");
  th1d_cosdiff.Draw();
  th1d_cosdiff.SetLineWidth(2);
  c1.Print("cosdiff.png");

}
