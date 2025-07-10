//#include "TComplex.h"

// 10 allows calculations of up to v2{10} and v3{6}
// 12 allows calculations of up to v2{12} and v3{8}
// 15 allows calculations of up to v2{14} and v3{10}
const int max_harmonic = 10;
// --- the max power dictates the larges number of particles that can be correlated
const int max_power = 10;

// for generic formulas
// will try to make this smarter, not dependent on global variable...
TComplex Qvector[max_harmonic][max_power];

TComplex Recursion(int, int*);
TComplex Recursion(int, int*, int, int);

TComplex get_flow_vector(const std::vector<double>& phi_angles, const int harmonic)
{
  TComplex Q(0.0,0.0);
  for ( auto it = phi_angles.begin(); it != phi_angles.end(); ++it )
    {
      double phi = *it;
      TComplex u(cos(harmonic*phi),sin(harmonic*phi));
      Q += u;
    }
  return Q;
}

TComplex get_weighted_flow_vector(const std::vector<std::pair<double,double>>& phi_weight, const int harmonic)
{
  TComplex Q(0.0,0.0);
  for ( auto it = phi_weight.begin(); it != phi_weight.end(); ++it )
    {
      double phi = it->first;
      double wgt = it->second;
      TComplex u(cos(harmonic*phi),sin(harmonic*phi));
      Q += wgt*u;
    }
  return Q;
}

std::array<TComplex,max_harmonic> get_flow_vectors(const std::vector<double>& phi_angles)
{
  std::array<TComplex,max_harmonic> allQ{};
  for ( int i = 0; i < max_harmonic; ++i )
    {
      allQ[i] = get_flow_vector(phi_angles,i);
    }
  return allQ;
}

// std::array<std::array<TComplex,max_harmonic>,max_power> get_weighted_flow_vectors(const std::vector<double>& phi_angles)
// {
//   std::array<TComplex,max_harmonic> allQ{};
//   for ( int i = 0; i < max_harmonic; ++i )
//     {
//       allQ[i] = get_flow_vector(phi_angles,i);
//     }
//   return allQ;
// }

// <cos(nphi)>
double calccosevent(const std::array<TComplex,max_harmonic>& allQ, int harmonic)
{
  double M = allQ[0].Re();
  if ( M < 1 ) return -9999;
  return allQ[harmonic].Re()/M;
}

// <sin(nphi)>
double calccsinevent(const std::array<TComplex,max_harmonic>& allQ, int harmonic)
{
  double M = allQ[0].Re();
  if ( M < 1 ) return -9999;
  return allQ[harmonic].Im()/M;
}

// <cos(n(phi1-phi2))>
double calc2event(const std::array<TComplex,max_harmonic>& allQ, int harmonic)
{
  double M = allQ[0].Re();
  if ( M < 2 ) return -9999;
  // ---
  TComplex Q = allQ[harmonic];
  // ---
  double numerator = Q.Rho2()-M;
  double denominator = M*(M-1);
  // ---
  return numerator/denominator;
}

// <cos(n(phi1a-phi2b))>
double calcSPevent(const std::array<TComplex,max_harmonic>& allQA, const std::array<TComplex,max_harmonic>& allQB, int harmonic)
{
  double MA = allQA[0].Re();
  double MB = allQA[0].Re();
  if ( MA < 1 || MB < 1 ) return -9999;
  // ---
  TComplex QA = allQA[harmonic];
  TComplex QBstar = TComplex::Conjugate(allQB[harmonic]);
  TComplex tc_numerator = QA*QBstar;
  // ---
  double numerator = tc_numerator.Re();
  double denominator = MA*MB;
  // ---
  return numerator/denominator;
}

// <cos(n(phi1+phi2))>
double calccossum2event(const std::array<TComplex,max_harmonic>& allQ, int harmonic)
{
  double M = allQ[0].Re();
  if ( M < 2 ) return -9999;
  // ---
  TComplex Q = allQ[harmonic];
  TComplex Q2 = allQ[2*harmonic];
  TComplex result = Q*Q - Q2;
  // ---
  double numerator = result.Re();
  double denominator = M*(M-1);
  // ---
  return numerator/denominator;
}

// <sin(n(phi1+phi2))>
double calcsinsum2event(const std::array<TComplex,max_harmonic>& allQ, int harmonic)
{
  double M = allQ[0].Re();
  if ( M < 2 ) return -9999;
  // ---
  TComplex Q = allQ[harmonic];
  TComplex Q2 = allQ[2*harmonic];
  TComplex result = Q*Q - Q2;
  // ---
  double numerator = result.Im();
  double denominator = M*(M-1);
  // ---
  return numerator/denominator;
}

// <cos(n(phi1-phi2-phi3))>
double calccos3event(const std::array<TComplex,max_harmonic>& allQ, int harmonic)
{
  double M = allQ[0].Re();
  if ( M < 3 ) return -9999;
  // ---
  TComplex Q = allQ[harmonic];
  TComplex Q2 = allQ[2*harmonic];
  TComplex Qstar = TComplex::Conjugate(allQ[harmonic]);
  TComplex Q2star = TComplex::Conjugate(allQ[2*harmonic]);
  TComplex result = Q*Qstar*Qstar - Q*Q2star;
  // ---
  double numerator = result.Re() - 2*(M-1)*Qstar.Re();
  double denominator = M*(M-1)*(M-2);
  // ---
  return numerator/denominator;
}

// <sin(n(phi1-phi2-phi3))>
double calcsin3event(const std::array<TComplex,max_harmonic>& allQ, int harmonic)
{
  double M = allQ[0].Re();
  if ( M < 3 ) return -9999;
  // ---
  TComplex Q = allQ[harmonic];
  TComplex Q2 = allQ[2*harmonic];
  TComplex Qstar = TComplex::Conjugate(allQ[harmonic]);
  TComplex Q2star = TComplex::Conjugate(allQ[2*harmonic]);
  TComplex result = Q*Qstar*Qstar - Q*Q2star;
  // ---
  double numerator = result.Im() - 2*(M-1)*Qstar.Im();
  double denominator = M*(M-1)*(M-2);
  // ---
  return numerator/denominator;
}

// <cos(n(phi1+phi2-phi3-phi4))>
double calc4event(const std::array<TComplex,max_harmonic>& allQ, int harmonic)
{
  double M = allQ[0].Re();
  if ( M < 4 ) return -9999;
  // ---
  TComplex Q = allQ[harmonic];
  TComplex Q2 = allQ[2*harmonic];
  TComplex Qstar = TComplex::Conjugate(allQ[harmonic]);
  TComplex Q2star = TComplex::Conjugate(allQ[2*harmonic]);
  TComplex tc_three = 2*Q2*Qstar*Qstar;
  // ---
  double one   = pow(Q.Rho2(),2);
  double two   = Q2.Rho2();
  double three = tc_three.Re();
  double four  = 2*(2*(M-2)*Q.Rho2());
  double five  = 2*(M*(M-3));
  cout << "one   is " << one << endl;
  cout << "two   is " << two << endl;
  cout << "three is " << three << endl;
  cout << "four  is " << four << endl;
  cout << "five  is " << five << endl;
  // ---
  double numerator = one + two - three - four + five;
  double denominator = M*(M-1)*(M-2)*(M-3);
  // ---
  return numerator/denominator;
}

// --- from generic forumulas ----------------------------------------------------

TComplex recQ(int n, int p)
{
  // Using the fact that Q{-n,p} = Q{n,p}^*.
  if(n>=0){return Qvector[n][p];}
  return TComplex::Conjugate(Qvector[-n][p]);
} // TComplex recQ(int n, int p)

TComplex Recursion(int n, int* harmonic)
{
  return Recursion(n,harmonic,1,0); // 1 and 0 are defaults from above
}

TComplex Recursion(int n, int* harmonic, int mult, int skip)
{
 // Calculate multi-particle correlators by using recursion (an improved faster version) originally developed by
 // Kristjan Gulbrandsen (gulbrand@nbi.dk).

  cout << "At the top, mult is " << mult << endl;

  int nm1 = n-1;
  TComplex c(recQ(harmonic[nm1], mult));
  if (nm1 == 0) return c;
  c *= Recursion(nm1, harmonic);
  if (nm1 == skip) return c;

  int multp1 = mult+1;
  int nm2 = n-2;
  int counter1 = 0;
  int hhold = harmonic[counter1];
  harmonic[counter1] = harmonic[nm2];
  harmonic[nm2] = hhold + harmonic[nm1];
  TComplex c2(Recursion(nm1, harmonic, multp1, nm2));
  int counter2 = n-3;
  while (counter2 >= skip) {
    harmonic[nm2] = harmonic[counter1];
    harmonic[counter1] = hhold;
    ++counter1;
    hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[nm2];
    harmonic[nm2] = hhold + harmonic[nm1];
    c2 += Recursion(nm1, harmonic, multp1, counter2);
    --counter2;
  }
  harmonic[nm2] = harmonic[counter1];
  harmonic[counter1] = hhold;

  cout << "At the bottom, mult is " << mult << endl;

  if (mult == 1) return c-c2;
  return c-double(mult)*c2;

}

// -------------------------------------------------------------------------------



void start_simple()
{

  std::vector<double> fake_angles = {0.00,0.05,3.00,3.05,4.00,4.05};
  std::vector<double> fake_weights = {1.01,1.02,1.03,1.04,1.05,1.06};

  std::vector<std::pair<double,double>> fake_angles_weights;
  // --- this is a little dumb but let's not worry about that too much for the time being
  TComplex check_fake_flow_2(0.0,0.0);
  TComplex check_fake_flow_3(0.0,0.0);
  TComplex check_weighted_fake_flow_2(0.0,0.0);
  TComplex check_weighted_fake_flow_3(0.0,0.0);
  for ( int i = 0; i < 6; ++i )
    {
      fake_angles_weights.push_back(std::make_pair(fake_angles[i],fake_weights[i]));
      check_fake_flow_2 += TComplex(cos(2*fake_angles[i]),sin(2*fake_angles[i]));
      check_fake_flow_3 += TComplex(cos(3*fake_angles[i]),sin(3*fake_angles[i]));
      check_weighted_fake_flow_2 += TComplex(cos(2*fake_angles[i]),sin(2*fake_angles[i]));
      check_weighted_fake_flow_3 += TComplex(cos(3*fake_angles[i]),sin(3*fake_angles[i]));
    }

  TComplex fake_flow_2 = get_flow_vector(fake_angles,2);
  TComplex fake_flow_3 = get_flow_vector(fake_angles,3);

  cout << "Fake flow flow vector for harmonic 2 is " << fake_flow_2 << endl;
  cout << "Fake flow flow vector for harmonic 3 is " << fake_flow_3 << endl;
  cout << "Fake flow check flow vector for harmonic 2 is " << check_fake_flow_2 << endl;
  cout << "Fake flow check flow vector for harmonic 3 is " << check_fake_flow_3 << endl;

  TComplex fake_weighted_flow_2 = get_weighted_flow_vector(fake_angles_weights,2);
  TComplex fake_weighted_flow_3 = get_weighted_flow_vector(fake_angles_weights,3);

  cout << "Fake flow weighted flow vector for harmonic 2 is " << fake_weighted_flow_2 << endl;
  cout << "Fake flow weighted flow vector for harmonic 3 is " << fake_weighted_flow_3 << endl;

  return;

  std::array<TComplex,max_harmonic> fake_flow_all = get_flow_vectors(fake_angles);

  for ( int i = 0; i < max_harmonic; ++i )
    {
      cout << "Fake flow flow vector for all harmonic " << i << " is " << fake_flow_all[i] << endl;
      for ( int j = 0; j < max_power; ++j )
        {
          Qvector[i][j] = fake_flow_all[i];
        }
    }

  // loop over the angles to do some direct calculations and compare
  int counter = 0;
  double dumb_cos2phi1phi2 = 0;
  double dumb_cos2phi1phi2_p = 0;
  // for ( auto it = fake_angles.begin(); it != fake_angles.end(); ++it )
  //   {
  //     double phi1 = *it;
  //     cout << "phi1 is " << phi1 << endl;
  //       for ( auto jt = fake_angles.begin(); jt != fake_angles.end(); ++jt )
  //         {
  //           double phi2 = *jt;
  //           if ( phi1 == phi2 ) continue;
  //           dumb_cos2phi1phi2 += cos(2*(phi1-phi2));
  //           ++counter;
  //         }
  //   }
  // --- dumb but quick and easy
  for ( int i = 0; i < fake_angles.size(); ++i )
    {
      double phi1 = fake_angles[i];
      cout << "phi1 is " << phi1 << endl;
      for ( int j = i+1; j < fake_angles.size(); ++j )
          {
            double phi2 = fake_angles[j];;
            dumb_cos2phi1phi2 += cos(2*(phi1-phi2));
            dumb_cos2phi1phi2_p += cos(2*(phi1+phi2));
            ++counter;
          }
    }
  dumb_cos2phi1phi2 /= counter;
  dumb_cos2phi1phi2_p /= counter;
  cout << "Direct calculation of cos(2(phi1-phi2)) is " << dumb_cos2phi1phi2 << endl;
  cout << "Direct calculation of cos(2(phi1+phi2)) is " << dumb_cos2phi1phi2_p << endl;

  double smart_cos2phi1phi2 = calc2event(fake_flow_all,2);
  double smart_cos2phi1phi2_p = calccossum2event(fake_flow_all,2);

  cout << "Flow vector based calculation of cos(2(phi1-phi2)) is " << smart_cos2phi1phi2 << endl;
  cout << "Flow vector based calculation of cos(2(phi1+phi2)) is " << smart_cos2phi1phi2_p << endl;


  int test2num[2]={2,-2};
  int test2den[2]={0,0};
  double super_smart_cos2phi1phi2 = Recursion(2,test2num).Re()/Recursion(2,test2den).Re();

  int test2pnum[2]={2,2};
  int test2pden[2]={0,0};
  double super_smart_cos2phi1phi2_p = Recursion(2,test2pnum).Re()/Recursion(2,test2pden).Re();

  cout << "Recursion based calculation of cos(2(phi1-phi2)) is " << super_smart_cos2phi1phi2 << endl;
  cout << "Recursion based calculation of cos(2(phi1+phi2)) is " << super_smart_cos2phi1phi2_p << endl;

  int test4num[4]={2,2,-2,-2};
  int test4den[4]={0,0,0,0};
  double super_smart_calc4 = Recursion(4,test4num).Re()/Recursion(4,test4den).Re();
  double smart_calc4 = calc4event(fake_flow_all,2);

  cout << "Flow vector based calculation of <4> (harmonic=2) is " << smart_calc4 << endl;
  cout << "Recursion based calculation of <4> (harmonic=2) is " << super_smart_calc4 << endl;

  int test6num[6]={2,2,2,-2,-2,-2};
  int test6den[6]={0,0,0,0,0,0};
  double super_smart_calc6 = Recursion(6,test6num).Re()/Recursion(6,test6den).Re();
  //double smart_calc6 = calc6event(fake_flow_all,2);

  //cout << "Flow vector based calculation of <6> (harmonic=2) is " << smart_calc6 << endl;
  cout << "Recursion based calculation of <6> (harmonic=2) is " << super_smart_calc6 << endl;

}
