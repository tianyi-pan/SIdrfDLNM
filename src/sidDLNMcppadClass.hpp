#ifndef SIDDLNMCPPADCLASS_HPP
#define SIDDLNMCPPADCLASS_HPP


// include CppAD
#include <cppad/cppad.hpp>

// autodiff include
// from https://github.com/awstringer1/varcomptest/blob/main/src/reml-ad.cpp
// #include "autodiff/common/meta.hpp"
// #include "autodiff/common/numbertraits.hpp"
// #include "autodiff/common/binomialcoefficient.hpp"
// #include "autodiff/common/vectortraits.hpp"
// #include "autodiff/forward/dual/dual.hpp"
// #include "autodiff/common/eigen.hpp"
// #include "autodiff/common/classtraits.hpp"
// #include "autodiff/forward/utils/derivative.hpp"
// #include "autodiff/forward/real/real.hpp"
// #include "autodiff/forward/utils/taylorseries.hpp"
// #include "autodiff/forward/real.hpp"
// #include "autodiff/forward/real/eigen.hpp"
// #include "autodiff/forward/utils/gradient.hpp"
// using namespace autodiff;



#include "defheader.h"


#include <cmath>

#include <iostream>
using namespace std;



// Eigen::VectorXd convertToDouble(const Vec input) {
//   int n = input.size();
//   Eigen::VectorXd output(n);

//   for (int i = 0; i < n; ++i) {
//       output(i) = CppAD::Value(input(i));
//   }
//   return output;
// }


// Eigen::MatrixXd convertToDouble(const Mat input) {
//   int nrow = input.rows();
//   int ncol = input.cols();
//   Eigen::MatrixXd output(nrow, ncol);

//   for (int i = 0; i < nrow; ++i) {
//       for (int j = 0; j < ncol; ++j) {
//           output(i, j) = CppAD::Value(input(i, j));
//       }
//   }

//   return output;
// }


// ************** PART 1: Define some functions **********************


// TODO: the bspline function evaluated at the points outside the boundaries are incorret!
// Bspline(l=0) = 0,0,0,0... It should be a linear function of l, not always equal to 0.
int knotindex(Scalar x,const Vec t) {
  int q = t.size();
  int k=0;
  if (x < t(0)) return -1;
  while(x>=t(k)){
    k++;
    if (k >= q) break;
  }

  return k-1;
}

Scalar weight(Scalar x, const Vec& t,int i,int k) {
  if (t(i+k-1) != t(i-1))
    return((x - t(i-1))/(t(i+k-1)-t(i-1)));
  return 0.;
}


Scalar Bspline(Scalar x, int j, const Vec& t,int p) {
  // Evaluate the jth B-spline
  // B_p(x) of order p (degree p-1) at x
  if (p==1)
    return(x>=t(j-1) && x<t(j+1-1));

  Scalar w1 = weight(x,t,j,p-1);
  Scalar w2 = weight(x,t,j+1,p-1);
  Scalar b1 = Bspline(x,j,t,p-1);
  Scalar b2 = Bspline(x,j+1,t,p-1);

  return w1*b1 + (1.-w2)*b2;
}

Scalar Bspline1st(Scalar x, int j, const Vec& t,int p) {
  // 1st derivative of the jth B-spline
  // https://stats.stackexchange.com/questions/258786/what-is-the-second-derivative-of-a-b-spline
  Scalar bb1 = Bspline(x,j+1,t,p-1);
  Scalar bb2 = Bspline(x,j,t,p-1);

  Scalar ww1 = 0.0;
  if (t(j+p-1) != t(j))
    ww1 = -1.0/(t(j+p-1)-t(j));
  Scalar ww2 = 0.0;
  if (t(j+p-2) != t(j-1))
    ww2 = 1.0/(t(j+p-2)-t(j-1));

  return (p-1.0) * (ww1 * bb1 + ww2 * bb2);
}

Scalar Bspline2nd(Scalar x, int j, const Vec& t,int p) {
  // 1st derivative of the jth B-spline
  // https://stats.stackexchange.com/questions/258786/what-is-the-second-derivative-of-a-b-spline
  Scalar bb1 = Bspline1st(x,j+1,t,p-1);
  Scalar bb2 = Bspline1st(x,j,t,p-1);

  Scalar ww1 = 0.0;
  if (t(j+p-1) != t(j))
    ww1 = -1.0/(t(j+p-1)-t(j));
  Scalar ww2 = 0.0;
  if (t(j+p-2) != t(j-1))
    ww2 = 1.0/(t(j+p-2)-t(j-1));

  return (p-1.0) * (ww1 * bb1 + ww2 * bb2);
}

Vec Bsplinevec(Scalar x, const Vec& t,int p) {
  int m = t.size() - p;
  Vec b(m);
  b.setZero();
  // int k = knotindex(x,t);
  // for (int i=(k-(p-1));i<k+1;i++)
  for (int i=0;i<m;i++)
    b(i) = Bspline(x,i+1,t,p);
  return b;
}

Vec BsplinevecCon(Scalar x, const Vec& t, int p, Mat Z) {
  int m = t.size() - p;
  Vec b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline(x,i+1,t,p);
  return Z.transpose()*b;
  // return b.transpose()*Z;
}

Vec Bsplinevec1st(Scalar x, const Vec& t,int p) {
  int m = t.size() - p;
  Vec b(m);
  b.setZero();
  // int k = knotindex(x,t);
  // for (int i=(k-(p-1));i<k+1;i++)
  for (int i=0;i<m;i++)
    b(i) = Bspline1st(x,i+1,t,p);
  return b;
}

Vec BsplinevecCon1st(Scalar x, const Vec& t, int p, Mat Z) {
  int m = t.size() - p;
  Vec b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline1st(x,i+1,t,p);
  return Z.transpose()*b;
  // return b.transpose()*Z;
}

Vec Bsplinevec2nd(Scalar x, const Vec& t,int p) {
  int m = t.size() - p;
  Vec b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline2nd(x,i+1,t,p);
  return b;
}


Vec BsplinevecCon2nd(Scalar x, const Vec& t,int p, Mat Z) {
  int m = t.size() - p;
  Vec b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline2nd(x,i+1,t,p);
  return Z.transpose()*b;
  // return b.transpose()*Z;
}



// Lanczos approximation
// Source code from https://github.com/brianmartens/BetaFunction/blob/master/BetaFunction/bmath.h
Scalar lanczos_lgamma(Scalar z) {
    const Scalar LG_g = 7.0;
    const int LG_N = 9;

    const Scalar ln_sqrt_2_pi = 0.91893853320467274178;

    Vec lct(LG_N+1);
    lct << 0.9999999999998099322768470047347,
    676.520368121885098567009190444019,
   -1259.13921672240287047156078755283,
    771.3234287776530788486528258894,
    -176.61502916214059906584551354,
     12.507343278686904814458936853,
    -0.13857109526572011689554707,
    9.984369578019570859563e-6,
    1.50563273514931155834e-7;

    Scalar sum;
    Scalar base;

    // To avoid if condition z < 0.5, we calculate gamma(z+1) which is equal to z*gamma(z).
    // WAS:
    // z = z - 1.0;
    // if (z < 0.5) {
    //   // Use Euler's reflection formula:
    //   // Gamma(z) = Pi / [Sin[Pi*z] * Gamma[1-z]];
    //   out = log(g_pi / sin(g_pi * z)) - lanczos_lgamma(1.0 - z);
    //   return out;
    // }
    // gamma(z) ...

    // New: indeed gamma(z+1)
    base = z + LG_g + 0.5;  // Base of the Lanczos exponential
    sum = 0;
    // We start with the terms that have the smallest coefficients and largest
    // denominator.
    for(int i=LG_N; i>=1; i--) {
      sum += lct[i] / (z + ((double) i));
    }
    sum += lct[0];
    Scalar gammazplus1 = ((ln_sqrt_2_pi + log(sum)) - base) + log(base)*(z+0.5);
    Scalar out = gammazplus1 - log(z);

    return out;
}

// 1st derivative of log gamma function
// https://math.stackexchange.com/questions/481253/differentiate-log-gamma-function
// VERY SLOW!!!
// Scalar lgamma1st (Scalar z) {
//   // Euler's constant https://en.wikipedia.org/wiki/Euler%27s_constant#
//   const Scalar Euler = 0.57721566490153286060651209008240243104215933593992;
//   Scalar out = -1.0 * Euler;
//   int K = 1e6;
//   for (int i = 1; i < K; i++) {
//     out += 1.0/i - 1.0/(i+z-1.0);
//   }
//   std::cout << "z" << CppAD::Value(z) << std::endl;
//   std::cout << "lgamma1st" << CppAD::Value(out) << std::endl;
//   return out;
// }



// https://github.com/tminka/lightspeed/blob/master/digamma.m
Scalar lgamma1st (Scalar x) {
  const Scalar pi = 3.141592653589793238462643383279;
  const Scalar large = 9.5;
  const Scalar d1 = -0.5772156649015328606065121;
  const Scalar d2 = pi*pi/6.0;
  const Scalar small = 1e-6;
  const Scalar s3 = 1.0/12.0;
  const Scalar s4 = 1.0/120.0;
  const Scalar s5 = 1.0/252.0;
  const Scalar s6 = 1.0/240.0;
  const Scalar s7 = 1.0/132.0;
  const Scalar s8 = 691.0/32760.0;
  const Scalar s9 = 1.0/12.0;
  const Scalar s10 = 3617.0/8160.0;

  // Use de Moivre's expansion if x >= large = 9.5
  // calculate lgamma1st(x+10)
  Scalar xplus10 = x + 10.0;
  Scalar y = 0.0;
  Scalar r = 1.0 / xplus10;
  y += log(xplus10) - 0.5 * r;
  r = r * r;
  y = y - r * ( s3 - r * ( s4 - r * (s5 - r * (s6 - r * s7))));

  // lgamma1st(x+10) = (1/x + 1/(x+1) + ... + 1/(x+9)) + lgamma1st(x)
  y = y - 1.0/x - 1.0/(x+1.0) - 1.0/(x+2.0) - 1.0/(x+3.0) - 1.0/(x+4.0) - 1.0/(x+5.0) - 1.0/(x+6.0) - 1.0/(x+7.0) - 1.0/(x+8.0) - 1/(x+9);

  return y;
}



// https://github.com/tminka/lightspeed/blob/master/trigamma.m

Scalar lgamma2nd (Scalar x) {
  const Scalar pi = 3.141592653589793238462643383279;
  const Scalar c = pi*pi/6;
  const Scalar c1 = -2.404113806319188570799476;
  const Scalar b2 =  1.0/6.0;
  const Scalar b4 = -1.0/30.0;
  const Scalar b6 =  1.0/42.0;
  const Scalar b8 = -1.0/30.0;
  const Scalar b10 = 5.0/66.0;

  // TO DO: % Reduce to trigamma(x+n) where ( X + N ) >= large.

  Scalar z = 1./(x*x);
  Scalar y = 0.5*z + (1.0 + z*(b2 + z*(b4 + z*(b6 + z*(b8 + z*b10))))) / x;

  // std::cout << "x" << CppAD::Value(x) << std::endl;
  // std::cout << "trigamma" << CppAD::Value(y) << std::endl;
  return y;
}

// **************** PART 2: g(mu) = DL term + linear term + smooth term *************************
class Modelcppad {
  // The sidDLNM model

private:
  // DATA
  const Vec y; // Response
  const Mat Sw; // penalty matrix for marginal lag
  const Mat SwR;
  const Mat Sf; // penalty matrix for marginal x
  const Mat SfR;
  const std::vector<Mat> B_inner_list;
  const Vec knots_f; // knots
  const Vec knots_w; // knots

  const Mat Xfix; // fixed effects
  const Mat Xrand; // random effects
  const Vec r; // rank of each smooth

  const Mat Zf;

  const Vec Xoffset; // offset


public:

  int n;
  int kE;
  int kl; // knots for marginal lag
  int kx; // knots for marginal x
  int L; // max Lag
  int kbetaR;
  int kbetaF;
  int mE;
  int p; // number of smooth terms in Xrand

  Mat Blag; // lag basis
  Vec Lseq; // sequence of lags

  const Mat Bindex; 

  // PARAMETERS
  Vec alpha_f;
  Vec con_index_par; // unconstrained index parameter
  Scalar log_theta;
  Scalar log_smoothing_f;
  Scalar log_smoothing_w;

  Vec betaF; // parameters for fixed effects
  Vec betaR; // parameters for random effects
  Vec logsmoothing; // log smoothing parameters for random effects

  // Components generated
  Mat B_inner; // B_inner = sum_m gamma_m B_inner_list.at(m)
  Scalar theta;
  Scalar smoothing_f;
  Scalar smoothing_w;
  Vec smoothing;
  Vec con_index_par_long;
  Scalar index_par_denominator;
  Vec index_par; // index parameter = con_index_par_long/index_par_denominator
  Mat Bf_matrix;
  Vec eta;
  Vec eta_remaining; // remaining terms = Xfix * betaF + Xrand * betaR
  Vec mu; // log(mu) = eta + eta_remaining + Xoffset
  Mat BtBindex;

  // Components for derivatives
  Mat dlogmu_dindex_par_mat;
  Vec dlogdensity_dmu_vec;
  Mat dmu_dindex_par_mat;
  Mat dindex_par_dcon_index_par_mat;

  Vec Bf;
  Vec Bf_tmp;
  Vec d2logdensity_dmudmu_vec;
  std::vector<Mat> d2index_par_dcon_index_pardcon_index_par_list;
  Mat he_index_par_mat;
  Mat he_alpha_f_index_par_mat;


  // gradient and hessian for updating alpha_f, betaR and betaF
  Mat he_inner_mat;



  // full hessian
  Mat he_alpha_f_mat;
  Mat he_betaR_mat;
  Mat he_betaF_mat;
  Mat he_con_index_par_mat;
  Mat he_alpha_f_con_index_par_mat;
  Mat he_alpha_f_betaF_mat;
  Mat he_alpha_f_betaR_mat;
  Mat he_betaR_betaF_mat;
  Mat he_con_index_par_betaF_mat;
  Mat he_con_index_par_betaR_mat;

  Mat he_s_u_mat;


  // for AIC
  Vec he_alpha_f_log_theta_vec;
  Vec he_betaR_log_theta_vec;
  Vec he_con_index_par_log_theta_vec;
  Vec he_betaF_log_theta_vec;
  Vec d2logdensity_dmudtheta_vec;
  Scalar he_log_theta_scalar;
  Mat I_mat;
  Mat I_alpha_f_mat;
  Mat I_betaR_mat;


  Mat K_alpha_f_mat;
  Mat K_betaR_mat;
  Mat K_betaF_mat;
  Mat K_index_par_mat;
  Mat K_con_index_par_mat;
  Vec K_log_theta_vec;

  Mat Kleft;
  Mat Khat; // Khat = Kleft * Kleft.transpose()


  // AD tape for LAML
  // CppAD::ADFun<double> gr;
  bool ifhastape = false;

  // Constructor
  Modelcppad(const Vec& y_,
            const std::vector<Mat>& B_inner_list_,
            const Vec& knots_f_,
            const Vec& knots_w_,
            const Mat& Sw_,
            const Mat& SwR_,
            const Mat& Sf_,
            const Mat& SfR_,
            const Mat& Xrand_,
            const Mat& Xfix_,
            const Mat& Zf_,
            const Vec& Xoffset_,
            const Vec& r_,
            const Mat& Bindex_,
            Vec& alpha_f_,
            Vec& con_index_par_,
            Scalar log_theta_,
            Scalar log_smoothing_f_,
            Scalar log_smoothing_w_,
            Vec& betaR_,
            Vec& betaF_,
            Vec& logsmoothing_) :
    y(y_), B_inner_list(B_inner_list_), knots_f(knots_f_), knots_w(knots_w_), Sw(Sw_), SwR(SwR_), Sf(Sf_), SfR(SfR_), Xrand(Xrand_), Xfix(Xfix_), Zf(Zf_), Xoffset(Xoffset_), r(r_), Bindex(Bindex_), 
    alpha_f(alpha_f_), con_index_par(con_index_par_), log_theta(log_theta_), log_smoothing_f(log_smoothing_f_), log_smoothing_w(log_smoothing_w_), betaR(betaR_), betaF(betaF_), logsmoothing(logsmoothing_) {

      n = y.size(); // sample size
      kE = alpha_f.size();
      kx = knots_f.size() - 4 - 1; // dim for marginal x. one for identifiability constraint (intercept).
      kl = knots_w.size() - 4; // dim for marginal lag
      L = B_inner_list.at(0).cols()-1; // max lag. ncol = L + 1
      kbetaR = betaR.size();
      kbetaF = betaF.size();
      mE = B_inner_list.size(); // number of exposures, which is the length of
      p = r.size();


      B_inner.resize(n, L+1);
      B_inner.setZero();

      Bf.resize(kE);
      Bf.setZero();
      Bf_tmp.resize(kE);
      Bf_tmp.setZero();

      Bf_matrix.resize(n, kE);

      theta = exp(log_theta);
      smoothing_f = exp(log_smoothing_f);
      smoothing_w = exp(log_smoothing_w);

      smoothing.resize(p);
      for (int i = 0; i < p; i++) smoothing(i) = exp(logsmoothing(i));

      BtBindex = Bindex.transpose() * Bindex;

      con_index_par_long.resize(mE); // con_index_par_long = c(1,con_index_par)
      con_index_par_long(0) = 1.0;
      for (int i = 0; i < (mE - 1); i++) {
        con_index_par_long(i + 1) = con_index_par(i);
      }
      index_par_denominator = sqrt(con_index_par_long.dot(BtBindex * con_index_par_long));
      index_par = Bindex * con_index_par_long/index_par_denominator;

      for (int i = 0; i < mE; i++) {
        B_inner += B_inner_list.at(i) * index_par(i);
      }



      Lseq.resize(L+1); // Set the length of Lseq to L+1
      for (int l = 0; l < (L+1); l++) {
        Lseq(l) = Scalar(l);
      }

      Blag.resize(kl, L+1); // basis for lag
      for (int j = 0; j < (L+1); j++) {
        Blag.col(j) = Bsplinevec(Lseq(j), knots_w, 4);
      }

      eta.resize(n);
      eta_remaining.resize(n);
      mu.resize(n);
      // START Bf_matrix and eta
      // B_inner.resize(n, L+1);
      for (int i = 0; i < n; i++) {
        Bf.setZero();
        Bf_tmp.setZero();
        for (int j = 0; j < (L+1); j++) {
          Vec bx = BsplinevecCon(B_inner(i, j), knots_f, 4, Zf);
          for (int ii = 0; ii < kx; ii++) {
            // TODO: optimize the computation of Bf_tmp
            Bf_tmp.segment(ii*kl, kl) = bx(ii) * Blag.col(j);
          }
          Bf += Bf_tmp;
        }
        Bf_matrix.row(i) = Bf;
        eta(i) = Bf.dot(alpha_f);
        eta_remaining(i) = Xfix.row(i).dot(betaF) + Xrand.row(i).dot(betaR);
        mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
      }
      // END Bf_matrix and eta





      // Initialize the derivative components and NegativeLogLikelihood

      dindex_par_dcon_index_par_mat = dindex_par_dcon_index_par();
      d2index_par_dcon_index_pardcon_index_par_list = d2index_par_dcon_index_pardcon_index_par();


      he_s_u_mat.resize(kE+mE-1+kbetaR+kbetaF, kE+mE-1+kbetaR+kbetaF);
      I_mat.resize(kE+mE-1+kbetaR+kbetaF + 1, kE+mE-1+kbetaR+kbetaF + 1);

      derivative_coef();
      derivative_he();
    }

  // Functions to set parameters
  void setAlphaF(const Vec alpha_f_) {
    alpha_f = alpha_f_;
  }

  void setCon_Index_Par(const Vec con_index_par_){
    con_index_par = con_index_par_;
  }
  void setBetaF(const Vec betaF_) {
    betaF = betaF_;
  }
  void setBetaR(const Vec betaR_) {
    betaR = betaR_;
  }
  void setLogTheta(const Scalar log_theta_) {
    log_theta = log_theta_;
    theta = exp(log_theta);
  }

  void setLogSmoothingF(const Scalar log_smoothing_f_) {
    log_smoothing_f = log_smoothing_f_;
    smoothing_f = exp(log_smoothing_f);
  }

  void setLogSmoothingW(const Scalar log_smoothing_w_) {
    log_smoothing_w = log_smoothing_w_;
    smoothing_w = exp(log_smoothing_w);
  }
  void setLogsmoothing(const Vec logsmoothing_) { // log smoothing parameters for remaining terms
    logsmoothing = logsmoothing_;
    for (int i = 0; i < p; i++) smoothing(i) = exp(logsmoothing(i));
  }


  void derivative_coef() {

    for (int i = 0; i < (mE - 1); i++) {
      con_index_par_long(i + 1) = con_index_par(i);
    }
    index_par_denominator = sqrt(con_index_par_long.dot(BtBindex * con_index_par_long));
    index_par = Bindex * con_index_par_long/index_par_denominator;
    B_inner.setZero();
    for (int i = 0; i < mE; i++) {
      B_inner += B_inner_list.at(i) * index_par(i);
    }


    // START Bf_matrix and eta
    // B_inner.resize(n, L+1);
    for (int i = 0; i < n; i++) {
      Bf.setZero();
      Bf_tmp.setZero();
      for (int j = 0; j < (L+1); j++) {
        Vec bx = BsplinevecCon(B_inner(i, j), knots_f, 4, Zf);
        for (int ii = 0; ii < kx; ii++) {
          // TODO: optimize the computation of Bf_tmp
          Bf_tmp.segment(ii*kl, kl) = bx(ii) * Blag.col(j);
        }
        Bf += Bf_tmp;
      }
      Bf_matrix.row(i) = Bf;
      eta(i) = Bf.dot(alpha_f);
      eta_remaining(i) = Xfix.row(i).dot(betaF) + Xrand.row(i).dot(betaR);
      mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
    }
    // END Bf_matrix and eta

    dindex_par_dcon_index_par_mat = dindex_par_dcon_index_par();
    d2index_par_dcon_index_pardcon_index_par_list = d2index_par_dcon_index_pardcon_index_par();


    dlogmu_dindex_par_mat = dlogmu_dindex_par();
    dlogdensity_dmu_vec = dlogdensity_dmu();
    dmu_dindex_par_mat = dmu_dindex_par();
 
  }

  // update full gradient and hessian of alpha_f, phi and betaR and betaF
  void derivative_he () {
    // obtain hessian
    d2logdensity_dmudmu_vec = d2logdensity_dmudmu();
    he_index_par_mat = he_index_par();
    he_alpha_f_index_par_mat = he_alpha_f_index_par();

    he_alpha_f_mat = he_alpha_f();
    he_betaR_mat = he_betaR();
    he_betaF_mat = he_betaF();
    he_con_index_par_mat = he_con_index_par();

    he_alpha_f_con_index_par_mat = he_alpha_f_con_index_par();
    he_alpha_f_betaF_mat = he_alpha_f_betaF();
    he_alpha_f_betaR_mat = he_alpha_f_betaR();
    he_betaR_betaF_mat = he_betaR_betaF();
    he_con_index_par_betaF_mat = he_con_index_par_betaF();
    he_con_index_par_betaR_mat = he_con_index_par_betaR();

    he_s_u_mat.setZero();


    he_s_u_mat.block(0, 0, kE, kE)  = he_alpha_f_mat;

    he_s_u_mat.block(0, kE, kE, mE-1) = he_alpha_f_con_index_par_mat;
    he_s_u_mat.block(kE, 0, mE-1, kE) = he_alpha_f_con_index_par_mat.transpose();
    he_s_u_mat.block(kE, kE, mE-1, mE-1) = he_con_index_par_mat;

    he_s_u_mat.block(kE+mE-1, kE+mE-1, kbetaR, kbetaR) = he_betaR_mat;
    he_s_u_mat.block(kE+kbetaR+mE-1, kE+kbetaR+mE-1, kbetaF, kbetaF) = he_betaF_mat;

    he_s_u_mat.block(0,kE+mE-1,kE,kbetaR) = he_alpha_f_betaR_mat;
    he_s_u_mat.block(kE,kE+mE-1,mE-1,kbetaR) = he_con_index_par_betaR_mat;
    he_s_u_mat.block(0,kE+mE-1+kbetaR,kE,kbetaF) = he_alpha_f_betaF_mat;
    he_s_u_mat.block(kE,kE+mE-1+kbetaR,mE-1,kbetaF) = he_con_index_par_betaF_mat;
    he_s_u_mat.block(kE+mE-1, kE+mE-1+kbetaR, kbetaR, kbetaF) = he_betaR_betaF_mat;

    he_s_u_mat.block(kE+mE-1,0,kbetaR,kE) = he_alpha_f_betaR_mat.transpose();
    he_s_u_mat.block(kE+mE-1,kE,kbetaR,mE-1) = he_con_index_par_betaR_mat.transpose();
    he_s_u_mat.block(kE+mE-1+kbetaR,0,kbetaF,kE) = he_alpha_f_betaF_mat.transpose();
    he_s_u_mat.block(kE+mE-1+kbetaR,kE,kbetaF,mE-1) = he_con_index_par_betaF_mat.transpose();
    he_s_u_mat.block(kE+mE-1+kbetaR, kE+mE-1, kbetaF, kbetaR) = he_betaR_betaF_mat.transpose();

  }

  void prepare_AIC () {

    I_mat.setZero();
    // obtain hessian
    d2logdensity_dmudmu_vec = d2logdensity_dmudmu();
    he_index_par_mat = he_index_par();
    he_alpha_f_index_par_mat = he_alpha_f_index_par();

    I_alpha_f_mat = I_alpha_f();
    I_betaR_mat = I_betaR();
    he_betaF_mat = he_betaF();
    he_con_index_par_mat = he_con_index_par();

    he_alpha_f_con_index_par_mat = he_alpha_f_con_index_par();
    he_alpha_f_betaF_mat = he_alpha_f_betaF();
    he_alpha_f_betaR_mat = he_alpha_f_betaR();
    he_betaR_betaF_mat = he_betaR_betaF();
    he_con_index_par_betaF_mat = he_con_index_par_betaF();
    he_con_index_par_betaR_mat = he_con_index_par_betaR();

   


    I_mat.block(0, 0, kE, kE)  = I_alpha_f_mat;

    I_mat.block(0, kE, kE, mE-1) = he_alpha_f_con_index_par_mat;
    I_mat.block(kE, 0, mE-1, kE) = he_alpha_f_con_index_par_mat.transpose();
    I_mat.block(kE, kE, mE-1, mE-1) = he_con_index_par_mat;

    I_mat.block(kE+mE-1, kE+mE-1, kbetaR, kbetaR) = I_betaR_mat;
    I_mat.block(kE+kbetaR+mE-1, kE+kbetaR+mE-1, kbetaF, kbetaF) = he_betaF_mat;

    I_mat.block(0,kE+mE-1,kE,kbetaR) = he_alpha_f_betaR_mat;
    I_mat.block(kE,kE+mE-1,mE-1,kbetaR) = he_con_index_par_betaR_mat;
    I_mat.block(0,kE+mE-1+kbetaR,kE,kbetaF) = he_alpha_f_betaF_mat;
    I_mat.block(kE,kE+mE-1+kbetaR,mE-1,kbetaF) = he_con_index_par_betaF_mat;
    I_mat.block(kE+mE-1, kE+mE-1+kbetaR, kbetaR, kbetaF) = he_betaR_betaF_mat;

    I_mat.block(kE+mE-1,0,kbetaR,kE) = he_alpha_f_betaR_mat.transpose();
    I_mat.block(kE+mE-1,kE,kbetaR,mE-1) = he_con_index_par_betaR_mat.transpose();
    I_mat.block(kE+mE-1+kbetaR,0,kbetaF,kE) = he_alpha_f_betaF_mat.transpose();
    I_mat.block(kE+mE-1+kbetaR,kE,kbetaF,mE-1) = he_con_index_par_betaF_mat.transpose();
    I_mat.block(kE+mE-1+kbetaR, kE+mE-1, kbetaF, kbetaR) = he_betaR_betaF_mat.transpose();

    d2logdensity_dmudtheta_vec = d2logdensity_dmudtheta();
    he_alpha_f_log_theta_vec = he_alpha_f_log_theta();
    he_betaR_log_theta_vec = he_betaR_log_theta();
    he_con_index_par_log_theta_vec = he_con_index_par_log_theta();
    he_betaF_log_theta_vec = he_betaF_log_theta();
    he_log_theta_scalar = he_log_theta();      


    I_mat.row(kE+mE-1+kbetaR+kbetaF) << he_alpha_f_log_theta_vec.transpose(), he_con_index_par_log_theta_vec.transpose(), he_betaR_log_theta_vec.transpose(), he_betaF_log_theta_vec.transpose(), he_log_theta_scalar;
    I_mat.col(kE+mE-1+kbetaR+kbetaF) = I_mat.row(kE+mE-1+kbetaR+kbetaF).transpose();
  }



  void prepare_AIC_proposed () {
    // Matrix K
    K_alpha_f_mat = K_alpha_f();
    K_betaR_mat = K_betaR();
    K_betaF_mat = K_betaF();
    K_index_par_mat = K_index_par();
    K_con_index_par_mat = K_con_index_par();
    K_log_theta_vec = K_log_theta();

    Kleft.resize(kE+mE-1+kbetaR+kbetaF+1, n);
    Kleft.setZero();

    Khat.resize(kE+mE-1+kbetaR+kbetaF+1, kE+mE-1+kbetaR+kbetaF+1);

    Kleft.block(0, 0, kE, n) = K_alpha_f_mat;
    Kleft.block(kE,0, mE-1,n) = K_con_index_par_mat;
    Kleft.block(kE+mE-1,0, kbetaR,n) = K_betaR_mat;
    Kleft.block(kE+mE-1+kbetaR,0, kbetaF,n) = K_betaF_mat;
    Kleft.row(kE+mE-1+kbetaR+kbetaF) = K_log_theta_vec.transpose();

    Khat = Kleft * Kleft.transpose();
  }



  // ********* Derivatives *************

  // FUNCTIONS
  // 1. density function
  // d log(exponential family density) / d mu
  Vec dlogdensity_dmu () {
    Vec out(n);
    for (int i = 0; i < n; i++) {
      out(i) = y(i) / mu(i) - (theta + y(i)) / (theta + mu(i));
    }
    return out;
  }
  // d^2 log(exponential family density) / d mu^2
  Vec d2logdensity_dmudmu () {
    Vec out(n);
    for (int i = 0; i < n; i++) {
      out(i) = - y(i) / pow(mu(i), 2) + (theta + y(i)) / pow(theta + mu(i), 2);
    }
    return out;
  }




  // 2. mean model
  // d log mu / d index_par
  Mat dlogmu_dindex_par () {
    Mat out(n, mE);
    Mat out_tmp(mE, kE);
    Vec Bf1st_tmp(kE);
    Vec bx(kx);

    for (int i = 0; i < n; i++) {
      out_tmp.setZero();
      Bf1st_tmp.setZero();
      for (int j = 0; j < (L+1); j++) {
        bx = BsplinevecCon1st(B_inner(i, j), knots_f, 4, Zf);
        for (int ii = 0; ii < kx; ii++) {
          // TODO: optimize the computation of Bf1st_tmp
          Bf1st_tmp.segment(ii*kl, kl) = bx(ii) * Blag.col(j);
        }
        for (int mm = 0; mm < mE; mm++) {
          out_tmp.row(mm) += B_inner_list.at(mm)(i, j) * Bf1st_tmp.transpose();
        }
      }
      out.row(i) = (out_tmp * alpha_f).transpose();
    }
    return out;
  }
  // d mu / d index_par
  Mat dmu_dindex_par () {
    Mat out(n, mE);
    for (int i = 0; i < n; i++) {
      out.row(i) = dlogmu_dindex_par_mat.row(i) * mu(i);
    }
    return out;
  }

  // 3. Re-parameterization
  // d index_par / d con_index_par
  Mat dindex_par_dcon_index_par () {

    Eigen::DiagonalMatrix<Scalar, Eigen::Dynamic> D(mE);
    D.diagonal().setConstant(1.0/index_par_denominator);
    Mat Ddense = D.toDenseMatrix();

    Mat deriv_g2 = (con_index_par_long * (BtBindex * con_index_par_long).transpose()) * pow(index_par_denominator, -3);

    Mat deriv_g = Ddense - deriv_g2;
    // Remove the first column
    Mat out = Bindex * (deriv_g.block(0, 1, mE, mE - 1));
    return out;
  }



  std::vector<Mat> d2index_par_dcon_index_pardcon_index_par () {
    std::vector<Mat> out;
    std::vector<Mat> outreal;
    
    Scalar tmp1 = pow(index_par_denominator, -3); // pow(tmp, -1.5)
    Scalar tmp2 = pow(index_par_denominator, -5); // pow(tmp, -2.5)
    Vec BtBindexcon_index_par_long = BtBindex * con_index_par_long;
    Mat outlarge(mE, mE);

    
    for (int s = 0; s < mE; s++) {
      if (s == 0) {
        outlarge = -1.0 * tmp1 * BtBindex + 3.0 * tmp2 * BtBindexcon_index_par_long * BtBindexcon_index_par_long.transpose();
      } else {
        Mat m1(mE, mE);
        m1.setZero();
        m1.row(s) = BtBindexcon_index_par_long.transpose()*tmp1;
        m1.col(s) = m1.col(s) + BtBindexcon_index_par_long*tmp1;
        Mat m2 = -1.0 * tmp1* BtBindex + 3.0 * tmp2 * BtBindexcon_index_par_long * BtBindexcon_index_par_long.transpose();
        outlarge = -1.0 * m1 + m2;
      }

      out.push_back(outlarge.block(1, 1, mE-1, mE-1));
    }

    Mat outrealtmp(mE-1, mE-1);
    for (int s = 0; s < mE; s++) {
      outrealtmp.setZero();
      for (int i = 0; i < mE; i++) {
        outrealtmp += Bindex(s, i) * out.at(i);
      }
      outreal.push_back(outrealtmp);
    }

    return outreal;
  }

  // components for K
  Mat K_alpha_f () {
    // dmu_df_mat: n * kE
    // dlogdensity_dmu_vec: n * 1
    // out: kE * n
    Mat out(kE, n);
    for (int i = 0; i < n; i++) {
      out.col(i) = mu(i) * Bf_matrix.row(i).transpose() * dlogdensity_dmu_vec(i);
    }
    return out;
  }

  
  Mat K_betaR () {
    // dmu_dbetaR_mat: n * kbetaR
    // dlogdensity_dmu_vec: n * 1
    // out: kbetaR * n
    Mat out(kbetaR, n);
    for (int i = 0; i < n; i++) {
      out.col(i) = mu(i) * Xrand.row(i).transpose() * dlogdensity_dmu_vec(i);
    }
    return out;
  }

  Mat K_betaF () {
    // dmu_dbetaF_mat: n * kbetaF
    // dlogdensity_dmu_vec: n * 1
    // out: kbetaF * n
    Mat out(kbetaF, n);
    for (int i = 0; i < n; i++) {
      out.col(i) = mu(i) * Xfix.row(i).transpose() * dlogdensity_dmu_vec(i);
    }
    return out;
  }


  Mat K_index_par () {
    // dmu_dindex_par_mat: n * mE
    // dlogdensity_dmu_vec: n * 1
    // out: mE * n
    Mat out(mE, n);
    for (int i = 0; i < n; i++) {
      out.col(i) = mu(i) * dlogmu_dindex_par_mat.row(i).transpose() * dlogdensity_dmu_vec(i);
    }
    return out;
  }

  Mat K_con_index_par () {
    // K_index_par_mat: mE * n
    // dindex_par_dcon_index_par_mat: mE * (mE-1)
    // out: (mE-1) * n
    return dindex_par_dcon_index_par_mat.transpose() * K_index_par_mat;
  }


  Vec K_log_theta () {
    // out: n * 1. DO NOT FORGET TO TRANSPOSE WHEN USING
    Vec out(n);
    for (int i = 0; i < n; i++) {
      out(i) = -1.0 * theta * (log_theta - log(theta + mu(i)) + (mu(i) - y(i))/(theta+mu(i)) + lgamma1st(theta+y(i)) - lgamma1st(theta));
    }
    return out;
  }




  // *** Hessian ***
  Mat he_alpha_f () {
    Mat out1(kE, kE);
    Mat out2(kE, kE);
    Mat muibftbf(kE, kE); 
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      muibftbf = mu(i) * Bf_matrix.row(i).transpose() * Bf_matrix.row(i);
      out1 += d2logdensity_dmudmu_vec(i) * mu(i) * muibftbf;
      out2 += dlogdensity_dmu_vec(i) * muibftbf;
    }
    return - out1 - out2 + smoothing_f*Sf + smoothing_w*Sw;
  }

  Mat I_alpha_f () {
    Mat out1(kE, kE);
    Mat out2(kE, kE);
    Mat muibftbf(kE, kE); 
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      muibftbf = mu(i) * Bf_matrix.row(i).transpose() * Bf_matrix.row(i);
      out1 += d2logdensity_dmudmu_vec(i) * mu(i) * muibftbf;
      out2 += dlogdensity_dmu_vec(i) * muibftbf;
    }
    return - out1 - out2;
  }



  Mat he_betaR () {
    Mat out1(kbetaR, kbetaR);
    Mat out2(kbetaR, kbetaR);
    Mat Ones(kbetaR, kbetaR); // identity matrix with diagonal smoothing
    Mat muixrtxr(kbetaR, kbetaR);
    out1.setZero();
    out2.setZero();
    Ones.setZero();
    for (int i = 0; i < n; i++) {
      muixrtxr = (mu(i) * Xrand.row(i).transpose() * Xrand.row(i));
      out1 += d2logdensity_dmudmu_vec(i) * mu(i) * muixrtxr;
      out2 += dlogdensity_dmu_vec(i) * muixrtxr;
    }
    int begin = 0;
    for (int i = 0; i < p; i++) {
      // Smooth Penalty
      double ri = CppAD::Value(r(i));
      int ki = static_cast<int>(ri);
      for (int j = 0; j < ki; j++) Ones(begin + j, begin + j) = smoothing(i);
      begin += ki;
    }
    return - out1 - out2 + Ones;
  }



  Mat I_betaR () {
    Mat out1(kbetaR, kbetaR);
    Mat out2(kbetaR, kbetaR);
    Mat muixrtxr(kbetaR, kbetaR);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      muixrtxr = (mu(i) * Xrand.row(i).transpose() * Xrand.row(i));
      out1 += d2logdensity_dmudmu_vec(i) * mu(i) * muixrtxr;
      out2 += dlogdensity_dmu_vec(i) * muixrtxr;
    }
    return - out1 - out2;
  }


  Mat he_betaF () {
    Mat out1(kbetaF, kbetaF);
    Mat out2(kbetaF, kbetaF);
    Mat muixftxf(kbetaF, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      muixftxf = mu(i) * Xfix.row(i).transpose() * Xfix.row(i);
      out1 += d2logdensity_dmudmu_vec(i) * mu(i) * muixftxf;
      out2 += dlogdensity_dmu_vec(i) * muixftxf;
    }
    return - out1 - out2;
  }


  Mat he_index_par () {
    Mat out1(mE, mE);
    Mat out2(mE, mE);
    Mat out_tmp(mE, mE);
    Mat xltxl_tmp(mE, mE);
    Vec bx;
    out1.setZero();
    out2.setZero();

    Vec Bf2nd_tmp(kE);
    double tmp;
    for (int i = 0; i < n; i++) {
      out_tmp.setZero();
      Bf2nd_tmp.setZero();
      for (int j = 0; j < (L+1); j++) {
        bx = BsplinevecCon2nd(B_inner(i, j), knots_f, 4, Zf);
        for (int ii = 0; ii < kx; ii++) {
          // TODO: optimize the computation of Bf2nd_tmp
          Bf2nd_tmp.segment(ii*kl, kl) = bx(ii) * Blag.col(j);
        }
        for (int mi = 0; mi < mE; mi++) {
          for (int mj = 0; mj < mE; mj++) {
            xltxl_tmp(mi,mj) = B_inner_list.at(mi)(i, j) * B_inner_list.at(mj)(i, j);
          }
        }

        out_tmp += xltxl_tmp * (Bf2nd_tmp.dot(alpha_f));
      }


      out1 += d2logdensity_dmudmu_vec(i) * dmu_dindex_par_mat.row(i).transpose() * dmu_dindex_par_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dindex_par_mat.row(i).transpose() * dlogmu_dindex_par_mat.row(i) + mu(i) * out_tmp);
    }

    return - out1 - out2;
  }


  Mat he_con_index_par () {
    Mat out1 = dindex_par_dcon_index_par_mat.transpose() * he_index_par_mat * dindex_par_dcon_index_par_mat;
    Mat out2(mE-1, mE-1);
    out2.setZero();
    Vec tmp = - dmu_dindex_par_mat.transpose() * dlogdensity_dmu_vec;
    for (int s = 0; s < mE; s++) {
      out2 = out2 + tmp(s) * d2index_par_dcon_index_pardcon_index_par_list.at(s);
    }
    return out1 + out2;
  }


  Mat he_alpha_f_index_par () {
    Mat out1(kE, mE);
    Mat out2(kE, mE);
    Mat out_mat(kE, mE);
    Mat out_tmp(mE, kE); // out_tmp = out_mat.transpose()
    Vec Bf1st_tmp(kE);
    Vec bx;
    out1.setZero();
    out2.setZero();
    Vec Bf1st;
    for (int i = 0; i < n; i++) {
      out_tmp.setZero();
      Bf1st_tmp.setZero();
      for (int j = 0; j < (L+1); j++) {
        bx = BsplinevecCon1st(B_inner(i, j), knots_f, 4, Zf);
        for (int ii = 0; ii < kx; ii++) {
          // TODO: optimize the computation of Bf1st_tmp
          Bf1st_tmp.segment(ii*kl, kl) = bx(ii) * Blag.col(j);
        }
        for (int mm = 0; mm < mE; mm++) {
          out_tmp.row(mm) += B_inner_list.at(mm)(i, j) * Bf1st_tmp.transpose();
        }
      }
      out_mat = out_tmp.transpose();

      out1 += d2logdensity_dmudmu_vec(i) * Bf_matrix.row(i).transpose() * dmu_dindex_par_mat.row(i) * mu(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * Bf_matrix.row(i).transpose() * dlogmu_dindex_par_mat.row(i) + mu(i)*out_mat);
    }
    return - out1 - out2;
  }


  Mat he_alpha_f_con_index_par () {
    Mat out = he_alpha_f_index_par_mat * dindex_par_dcon_index_par_mat;
    return out;
  }

  Mat he_alpha_f_betaF () {
    Mat out1(kE, kbetaF);
    Mat out2(kE, kbetaF);
    Mat muibgtxf(kE, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      muibgtxf = Bf_matrix.row(i).transpose() * (Xfix.row(i) * mu(i));
      out1 += d2logdensity_dmudmu_vec(i) * muibgtxf * mu(i);
      out2 += dlogdensity_dmu_vec(i) * muibgtxf;
    }
    return - out1 - out2;
  }
  Mat he_alpha_f_betaR () {
    Mat out1(kE, kbetaR);
    Mat out2(kE, kbetaR);
    Mat muibftxr(kE, kbetaR);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      muibftxr = Bf_matrix.row(i).transpose() * (Xrand.row(i) * mu(i));
      out1 += d2logdensity_dmudmu_vec(i) * muibftxr * mu(i);
      out2 += dlogdensity_dmu_vec(i) * muibftxr;
    }
    return - out1 - out2;
  }

  Mat he_betaR_betaF () {
    Mat out1(kbetaR, kbetaF);
    Mat out2(kbetaR, kbetaF);
    Mat muixrtxf(kbetaR, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      muixrtxf = Xrand.row(i).transpose() * (Xfix.row(i) * mu(i));
      out1 += d2logdensity_dmudmu_vec(i) * muixrtxf * mu(i);
      out2 += dlogdensity_dmu_vec(i) * muixrtxf;
    }
    return - out1 - out2;
  }


  Mat he_con_index_par_betaF () {
    Mat out1(mE-1, kbetaF);
    Mat out2(mE-1, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += dindex_par_dcon_index_par_mat.transpose() * (d2logdensity_dmudmu_vec(i) * dmu_dindex_par_mat.row(i).transpose() * (Xfix.row(i) * mu(i)));
      out2 += dindex_par_dcon_index_par_mat.transpose() * (dlogdensity_dmu_vec(i) * dlogmu_dindex_par_mat.row(i).transpose() * (Xfix.row(i) * mu(i)));
    }
    return - out1 - out2;
  }



  Mat he_con_index_par_betaR () {
    Mat out1(mE-1, kbetaR);
    Mat out2(mE-1, kbetaR);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += dindex_par_dcon_index_par_mat.transpose() * (d2logdensity_dmudmu_vec(i) * dmu_dindex_par_mat.row(i).transpose() * (Xrand.row(i) * mu(i)));
      out2 += dindex_par_dcon_index_par_mat.transpose() * (dlogdensity_dmu_vec(i) * dlogmu_dindex_par_mat.row(i).transpose() * (Xrand.row(i) * mu(i)));
    }
    return - out1 - out2;
  }

  // For AIC
  Vec d2logdensity_dmudtheta () {
    Vec out(n);
    for (int i = 0; i < n; i++) {
      out(i) = (y(i) - mu(i)) / pow(theta+mu(i), 2);
    }
    return out;
  }

  Vec he_alpha_f_log_theta () {
    Mat dmu_df_mat(n, kE);
    for (int i = 0; i < n; i++) {
      dmu_df_mat.row(i) = Bf_matrix.row(i) * mu(i);
    }

    return -1.0*dmu_df_mat.transpose() * d2logdensity_dmudtheta_vec * theta;
  }

  Vec he_betaR_log_theta () {
    Mat dmu_dbetaR_mat(n, kbetaR);
    for (int i = 0; i < n; i++) {
      dmu_dbetaR_mat.row(i) = Xrand.row(i) * mu(i);
    }

    return -1.0*dmu_dbetaR_mat.transpose() * d2logdensity_dmudtheta_vec * theta;
  }


  Vec he_con_index_par_log_theta () {
    return -1.0*theta * dindex_par_dcon_index_par_mat.transpose() * ( dmu_dindex_par_mat.transpose() * d2logdensity_dmudtheta_vec );
  }
  
  Vec he_betaF_log_theta () {
    Mat dmu_dbetaF_mat(n, kbetaF);
    for (int i = 0; i < n; i++) {
      dmu_dbetaF_mat.row(i) = Xfix.row(i) * mu(i);
    }
    return -1.0*dmu_dbetaF_mat.transpose() * d2logdensity_dmudtheta_vec * theta;
  }

  Scalar he_log_theta () {
    Scalar dlogdensity_dtheta_scalar = 0.0;
    for (int i = 0; i < n; i++) {
      dlogdensity_dtheta_scalar += log_theta - log(theta + mu(i)) + (mu(i) - y(i))/(theta+mu(i)) + lgamma1st(theta+y(i)) - lgamma1st(theta);
    }
    Scalar d2logdensity_dthetadtheta_scalar = 0.0;
    for (int i = 0; i < n; i++) {
      d2logdensity_dthetadtheta_scalar += 1/theta - 1/(theta + mu(i)) - (mu(i) - y(i)) / ((theta + mu(i))*(theta + mu(i))) + lgamma2nd(y(i) + theta) - lgamma2nd(theta);
    }

    return -1.0 * theta * dlogdensity_dtheta_scalar - theta * theta * d2logdensity_dthetadtheta_scalar;
  }


  // *********** LAML ***********
  Scalar logdetH05() {
    Scalar logdetH05 = 0.0;
    // Scalar eps = 1e-8; // small value to avoid log(0)
    Eigen::FullPivLU<Mat> lu(he_s_u_mat);
    Mat LU = lu.matrixLU();
    // Scalar c = lu.permutationP().determinant(); // -1 or 1
    Scalar lii;
    for (int i = 0; i < (mE-1+kE+kbetaR+kbetaF); i++) {
      lii = LU(i,i);
      // std::cout << "lii : " << CppAD::Value(lii) << std::endl;
      // std::cout << "c : " << CppAD::Value(c) << std::endl;
      // if (lii < 0.0) c *= -1;

      logdetH05 += log(CppAD::abs(lii));
      // logdetH05 += log(CppAD::abs(lii) + eps);  // add eps to avoid log(0)
      // logdetH05 += log(std::max(lii, eps));
    }
    // logdetH05 += log(c);



    return logdetH05/2.0;

  }
};





#endif
