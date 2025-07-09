//#include "TComplex.h"

// 10 allows calculations of up to v2{10} and v3{6}
// 12 allows calculations of up to v2{12} and v3{8}
const int max_harmonic = 10;

// some of the generic formulas make use of weights to various powers
// we're not going to worry about that for the time being
const int max_power = 1;

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

std::array<TComplex, max_harmonic> get_flow_vectors(const std::vector<double>& phi_angles)
{
  std::array<TComplex, max_harmonic> allQ{};
  for ( int i = 0; i < max_harmonic; ++i )
    {
      allQ[i] = get_flow_vector(phi_angles,i);
    }
  return allQ;
}

double calc2event(const std::array<TComplex, max_harmonic>& allQ, int harmonic)
{
  TComplex Q = allQ[harmonic];
  TComplex Qstar = TComplex::Conjugate(allQ[harmonic]);
  TComplex tc_numerator = Q*Qstar;
  double numerator = tc_numerator.Re();
  double M = allQ[0].Re();
  double denominator = M*M(-1);
  return numerator/denominator;
}

double calcSPevent(const std::array<TComplex, max_harmonic>& allQA, std::array<TComplex, max_harmonic>& allQB, int harmonic)
{
  TComplex QA = allQA[harmonic];
  TComplex QBstar = TComplex::Conjugate(allQB[harmonic]);
  TComplex tc_numerator = QA*QBstar;
  double MA = allQA[0].Re();
  double MB = allQA[0].Re();
  double numerator = tc_numerator.Re();
  double denominator = MA*MB;
  return numerator/denominator;
}

double calc4event(const std::array<TComplex, max_harmonic>& allQ, int harmonic)
{
  TComplex Q = allQ[harmonic];
  TComplex Q2 = allQ[2*harmonic];
  TComplex Qstar = TComplex::Conjugate(allQ[harmonic]);
  TComplex Q2star = TComplex::Conjugate(allQ[2*harmonic]);

  TComplex tc_three = Q2*Qstar*Qstar;

  double one   = pow(Q.Rho2(),2);
  double two   = Q2.Rho2();
  double three = tc_three.Re();
  double four  = 2*(2*(M-2)*Q.Rho2());
  double five  = 2*(M*(M-3));

  double M = allQ[0].Re();

  double numerator = one + two - three - four + five;
  double denominator = M*(M-1)*(M-2)*(M-3);
  return numerator/denominator;
}



void start_simple()
{

}
