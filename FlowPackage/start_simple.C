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

// --- wondering if I want to change this...
std::array<TComplex,max_harmonic> get_flow_vectors(const std::vector<double>& phi_angles)
{
  std::array<TComplex,max_harmonic> allQ{};
  for ( int i = 0; i < max_harmonic; ++i )
    {
      allQ[i] = get_flow_vector(phi_angles,i);
    }
  return allQ;
}

// --- because I think I want to do this differently
std::array<std::array<TComplex,max_harmonic>,max_power> get_weighted_flow_vectors(const std::vector<std::pair<double,double>>& phi_weight)
{
  std::vector<std::pair<double,double>> new_phi_weight = phi_weight;
  std::array<std::array<TComplex,max_harmonic>,max_power> allQ{{}};
  for ( int i = 0; i < max_harmonic; ++i )
    {
      for ( int j = 0; j < max_power; ++j )
        {
          // weight needs to be raised to power j
          auto it_old = phi_weight.begin();
          auto it_new = new_phi_weight.begin();
          for ( ; it_old != phi_weight.end() && it_new != new_phi_weight.end(); ++it_old, ++it_new )
            {
              it_new->second = pow(it_old->second,j);
            }
          // use the new weight to get the weighted flow vector
          allQ[i][j] = get_weighted_flow_vector(new_phi_weight,i);
        }
    }
  return allQ;
}

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
  // ---
  double numerator = one + two - three - four + five;
  double denominator = M*(M-1)*(M-2)*(M-3);
  // ---
  return numerator/denominator;
}

// <cos(n(phi1+phi2+phi3-phi4-phi5-phi6))>
double calc6_event_jamie(TComplex& qn, TComplex& q2n, TComplex& q3n, float M)
{

  if ( M < 6 ) return -9999;

  // TComplex qn, q2n, q3n;
  // qn = TComplex(Q2x,Q2y);
  // q2n = TComplex(Q4x,Q4y);
  // q3n = TComplex(Q6x,Q6y);

  TComplex temp1;

  // first term
  // |Qn|^6 + 9*|Q2n|^2|Qn|^2 - 6 x Re[Q2n x Qn x Qn* x Qn* x Qn*] / (Mx(M-1)x(M-2)x(M-3)x(M-4)x(M-5)
  double term1a = TMath::Power((qn*TComplex::Conjugate(qn)),3);
  double term1b = 9.0 * q2n*TComplex::Conjugate(q2n) * qn*TComplex::Conjugate(qn);
  temp1 = q2n * qn * TComplex::Conjugate(qn) * TComplex::Conjugate(qn) * TComplex::Conjugate(qn);
  double term1c = -6.0 * temp1.Re();
  double term1 = (term1a+term1b+term1c)/(M*(M-1)*(M-2)*(M-3)*(M-4)*(M-5));

  // second term
  // 4 * [Re[Q3nQn*Qn*Qn*] - 3 Re[Q3nQ2n*Qn*]] / (M(M-1)(M-2)(M-3)(M-4)(M-5)
  temp1 = q3n * TComplex::Conjugate(qn) * TComplex::Conjugate(qn) * TComplex::Conjugate(qn);
  double term2a = temp1.Re();
  temp1 = q3n * TComplex::Conjugate(q2n) * TComplex::Conjugate(qn);
  double term2b = -3.0 * temp1.Re();
  double term2 = 4.0 * (term2a+term2b)/(M*(M-1)*(M-2)*(M-3)*(M-4)*(M-5));

  // third term
  // +2 * (9*(M-4)*Re[Q2nQn*qn*] + 2 |Q3n|^2) / ((M(M-1)(M-2)(M-3)(M-4)(M-5))
  temp1 = q2n*TComplex::Conjugate(qn)*TComplex::Conjugate(qn);
  double term3a = 9.0*(M-4)*temp1.Re();
  double term3b = 2.0*q3n*TComplex::Conjugate(q3n);
  double term3 = 2.0 * (term3a + term3b) / (M*(M-1)*(M-2)*(M-3)*(M-4)*(M-5));

  // fourth term
  //double term4 = -9.0 * (TMath::Power(qn*TComplex::Conjugate(qn),2)+q2n*TComplex::Conjugate(q2n)) / (M*(M-1)*(M-2)*(M-3)*(M-5));
  double term4 = -9.0 * (TMath::Power(qn*TComplex::Conjugate(qn),2)+q2n*TComplex::Conjugate(q2n)) ;
  term4 /= (M*(M-1)*(M-2)*(M-3)*(M-5));

  // fifth term
  //double term5 = 18.0 * qn*TComplex::Conjugate(qn) / (M*(M-1)*(M-3)*(M-4));
  double term5 = 18.0 * qn*TComplex::Conjugate(qn) ;
  term5 /=  (M*(M-1)*(M-3)*(M-4));

  // sixth term
  double term6 = -6.0/((M-1)*(M-2)*(M-3));

  // cos(n(phi1+phi2+phi3-phi4-phi5-phi6))
  double six = term1 + term2 + term3 + term4 + term5 + term6;

  return six;

}

// <cos(n(phi1+phi2+phi3-phi4-phi5-phi6))>
double calc6event(const std::array<TComplex,max_harmonic>& allQ, int harmonic)
{

  double M = allQ[0].Re();
  if ( M < 6 ) return -9999;
  // ---
  TComplex Q = allQ[harmonic];
  TComplex Q2 = allQ[2*harmonic];
  TComplex Q3 = allQ[3*harmonic];
  TComplex Qstar = TComplex::Conjugate(allQ[harmonic]);
  TComplex Q2star = TComplex::Conjugate(allQ[2*harmonic]);
  TComplex Q3star = TComplex::Conjugate(allQ[3*harmonic]);

  TComplex temp1;

  // first term
  // |Qn|^6 + 9*|Q2n|^2|Qn|^2 - 6 x Re[Q2n x Qn x Qn* x Qn* x Qn*] / (Mx(M-1)x(M-2)x(M-3)x(M-4)x(M-5)
  double term1a = TMath::Power((Q*Qstar),3);
  double term1b = 9.0 * Q2*Q2star * Q*Qstar;
  temp1 = Q2 * Q * Qstar * Qstar * Qstar;
  double term1c = -6.0 * temp1.Re();
  double term1 = (term1a+term1b+term1c)/(M*(M-1)*(M-2)*(M-3)*(M-4)*(M-5));

  // second term
  // 4 * [Re[Q3nQn*Qn*Qn*] - 3 Re[Q3nQ2n*Qn*]] / (M(M-1)(M-2)(M-3)(M-4)(M-5)
  temp1 = Q3 * Qstar * Qstar * Qstar;
  double term2a = temp1.Re();
  temp1 = Q3 * Q2star * Qstar;
  double term2b = -3.0 * temp1.Re();
  double term2 = 4.0 * (term2a+term2b)/(M*(M-1)*(M-2)*(M-3)*(M-4)*(M-5));

  // third term
  // +2 * (9*(M-4)*Re[Q2nQn*qn*] + 2 |Q3n|^2) / ((M(M-1)(M-2)(M-3)(M-4)(M-5))
  temp1 = Q2*Qstar*Qstar;
  double term3a = 9.0*(M-4)*temp1.Re();
  double term3b = 2.0*Q3*Q3star;
  double term3 = 2.0 * (term3a + term3b) / (M*(M-1)*(M-2)*(M-3)*(M-4)*(M-5));

  // fourth term
  //double term4 = -9.0 * (TMath::Power(Q*Qstar,2)+Q2*Q2star) / (M*(M-1)*(M-2)*(M-3)*(M-5));
  double term4 = -9.0 * (TMath::Power(Q*Qstar,2)+Q2*Q2star) ;
  term4 /= (M*(M-1)*(M-2)*(M-3)*(M-5));

  // fifth term
  //double term5 = 18.0 * Q*Qstar / (M*(M-1)*(M-3)*(M-4));
  double term5 = 18.0 * Q*Qstar ;
  term5 /=  (M*(M-1)*(M-3)*(M-4));

  // sixth term
  double term6 = -6.0/((M-1)*(M-2)*(M-3));

  // cos(n(phi1+phi2+phi3-phi4-phi5-phi6))
  double six = term1 + term2 + term3 + term4 + term5 + term6;

  return six;

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

  //std::cout << "At the top, mult is " << mult << std::endl;

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

  //std::cout << "At the bottom, mult is " << mult << std::endl;

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
      check_weighted_fake_flow_2 += fake_weights[i]*TComplex(cos(2*fake_angles[i]),sin(2*fake_angles[i]));
      check_weighted_fake_flow_3 += fake_weights[i]*TComplex(cos(3*fake_angles[i]),sin(3*fake_angles[i]));
    }

  TComplex fake_flow_2 = get_flow_vector(fake_angles,2);
  TComplex fake_flow_3 = get_flow_vector(fake_angles,3);

  std::cout << "Fake flow flow vector for harmonic 2 is " << fake_flow_2 << std::endl;
  std::cout << "Fake flow flow vector for harmonic 3 is " << fake_flow_3 << std::endl;
  std::cout << "Fake flow check flow vector for harmonic 2 is " << check_fake_flow_2 << std::endl;
  std::cout << "Fake flow check flow vector for harmonic 3 is " << check_fake_flow_3 << std::endl;

  TComplex fake_weighted_flow_2 = get_weighted_flow_vector(fake_angles_weights,2);
  TComplex fake_weighted_flow_3 = get_weighted_flow_vector(fake_angles_weights,3);

  std::cout << "Fake flow weighted flow vector for harmonic 2 is " << fake_weighted_flow_2 << std::endl;
  std::cout << "Fake flow weighted flow vector for harmonic 3 is " << fake_weighted_flow_3 << std::endl;
  std::cout << "Fake flow check weighted flow vector for harmonic 2 is " << check_weighted_fake_flow_2 << std::endl;
  std::cout << "Fake flow check weighted flow vector for harmonic 3 is " << check_weighted_fake_flow_3 << std::endl;

  //return;

  std::array<TComplex,max_harmonic> fake_flow_all = get_flow_vectors(fake_angles);

  std::array<std::array<TComplex,max_harmonic>,max_power> fake_weighted_flow_all = get_weighted_flow_vectors(fake_angles_weights);

  std::cout << "Without weights: " << std::endl;
  std::cout << fake_flow_2 << std::endl;
  std::cout << fake_flow_all[2] << std::endl;
  std::cout << fake_weighted_flow_all[2][0] << std::endl;

  std::cout << "With linear weights: " << std::endl;
  std::cout << fake_weighted_flow_2 << std::endl;
  // --- this is giving the wrong result at the moment, for all powers...
  // --- second index = 1 should give linear weight
  std::cout << fake_weighted_flow_all[2][1] << std::endl;
  std::cout << "With higher order weights: " << std::endl;
  std::cout << fake_weighted_flow_all[2][2] << std::endl;
  std::cout << fake_weighted_flow_all[2][3] << std::endl;
  std::cout << fake_weighted_flow_all[2][4] << std::endl;

  //return;

  for ( int i = 0; i < max_harmonic; ++i )
    {
      std::cout << "Fake flow flow vector for all harmonic " << i << " is " << fake_flow_all[i] << std::endl;
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
  //     std::cout << "phi1 is " << phi1 << std::endl;
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
      std::cout << "phi1 is " << phi1 << std::endl;
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
  std::cout << "Direct calculation of cos(2(phi1-phi2)) is " << dumb_cos2phi1phi2 << std::endl;
  std::cout << "Direct calculation of cos(2(phi1+phi2)) is " << dumb_cos2phi1phi2_p << std::endl;

  double smart_cos2phi1phi2 = calc2event(fake_flow_all,2);
  double smart_cos2phi1phi2_p = calccossum2event(fake_flow_all,2);

  std::cout << "Flow vector based calculation of cos(2(phi1-phi2)) is " << smart_cos2phi1phi2 << std::endl;
  std::cout << "Flow vector based calculation of cos(2(phi1+phi2)) is " << smart_cos2phi1phi2_p << std::endl;


  int test2num[2]={2,-2};
  int test2den[2]={0,0};
  double super_smart_cos2phi1phi2 = Recursion(2,test2num).Re()/Recursion(2,test2den).Re();

  int test2pnum[2]={2,2};
  int test2pden[2]={0,0};
  double super_smart_cos2phi1phi2_p = Recursion(2,test2pnum).Re()/Recursion(2,test2pden).Re();

  std::cout << "Recursion based calculation of cos(2(phi1-phi2)) is " << super_smart_cos2phi1phi2 << std::endl;
  std::cout << "Recursion based calculation of cos(2(phi1+phi2)) is " << super_smart_cos2phi1phi2_p << std::endl;

  int test4num[4]={2,2,-2,-2};
  int test4den[4]={0,0,0,0};
  double super_smart_calc4 = Recursion(4,test4num).Re()/Recursion(4,test4den).Re();
  double smart_calc4 = calc4event(fake_flow_all,2);

  std::cout << "Flow vector based calculation of <4> (harmonic=2) is " << smart_calc4 << std::endl;
  std::cout << "Recursion based calculation of <4> (harmonic=2) is " << super_smart_calc4 << std::endl;

  int test6num[6]={2,2,2,-2,-2,-2};
  int test6den[6]={0,0,0,0,0,0};
  std::cout << "Starting recursion calculation for 6..." << std::endl;
  double numerator6 = Recursion(6,test6num).Re();
  std::cout << "Now need denominator..." << std::endl;
  double denominator6 = Recursion(6,test6den).Re();
  double super_smart_calc6 = numerator6/denominator6;
  double jamie_smart_calc6 = calc6_event_jamie(fake_flow_all[2],fake_flow_all[4],fake_flow_all[6],fake_flow_all[0]);
  double smart_calc6 = calc6event(fake_flow_all,2);

  std::cout << "Jamie's flow vector based calculation of <6> (harmonic=2) is " << jamie_smart_calc6 << std::endl;
  std::cout << "Flow vector based calculation of <6> (harmonic=2) is " << smart_calc6 << std::endl;
  std::cout << "Recursion based calculation of <6> (harmonic=2) is " << super_smart_calc6 << std::endl;

}
