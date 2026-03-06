#ifndef SIDDLNMEIGENCLASS_HPP
#define SIDDLNMEIGENCLASS_HPP




#include <lambda_lanczos.hpp>
using lambda_lanczos::LambdaLanczos;


#include <cmath>

#include <iostream>
using namespace std;



// ************** PART 1: Define some functions **********************
void choleskyAD(Eigen::MatrixXd& L) {
  // Eigen::MatrixXd will be overwritten; lower triangle will be its Cholesky. Only the lower triangle is computed/stored
  int s = L.cols();
  for (int k = 0; k < s; k++) {
    // (a) define pivot
    L(k,k) = sqrt(L(k,k));

    // (b) adjust lead column
    for (int j = k+1; j < s; j++) L(j,k) /= L(k,k);
    for (int j = k+1; j < s; j++)
      for (int i = j; i < s; i++) L(i,j) -= L(i,k) * L(j,k);
  }
}

Eigen::MatrixXd invertL(Eigen::MatrixXd &L) {
  // inverse of a lower triangular matrix
  int n = L.cols();
  Eigen::MatrixXd M(n, n);
  M.setZero();
  for (int i = 0; i < n; i++)
  {
    M(i,i) = 1.0 / L(i,i);
    for (int j = 0; j < i; j++)
    {
      for (int k = j; k < i; k++) M(i,j) += L(i,k) * M(k,j);
      M(i,j) = -M(i,j) / L(i,i);
    }
  }
  return M;
}

// check whether there is nan in the input vector
bool hasNaN(Eigen::VectorXd vec) {
    for (int i = 0; i < vec.size(); i++) {
      if( std::isnan(vec(i))) return true; // has nan
    }
    return false; // No nan
}

// TODO: the bspline function evaluated at the points outside the boundaries are incorret!
// Bspline(l=0) = 0,0,0,0... It should be a linear function of l, not always equal to 0.
int knotindexEigen(double x,const Eigen::VectorXd t) {
  int q = t.size();
  int k=0;
  if (x < t(0)) return -1;
  while(x>=t(k)){
    k++;
    if (k >= q) break;
  }

  return k-1;
}

double weight(double x, const Eigen::VectorXd& t,int i,int k) {
  if (t(i+k-1) != t(i-1))
    return((x - t(i-1))/(t(i+k-1)-t(i-1)));
  return 0.;
}


double Bspline(double x, int j, const Eigen::VectorXd& t,int p) {
  // Evaluate the jth B-spline
  // B_p(x) of order p (degree p-1) at x
  if (p==1)
    return(x>=t(j-1) && x<t(j+1-1));

  double w1 = weight(x,t,j,p-1);
  double w2 = weight(x,t,j+1,p-1);
  double b1 = Bspline(x,j,t,p-1);
  double b2 = Bspline(x,j+1,t,p-1);

  return w1*b1 + (1.-w2)*b2;
}

double Bspline1st(double x, int j, const Eigen::VectorXd& t,int p) {
  // 1st derivative of the jth B-spline
  // https://stats.stackexchange.com/questions/258786/what-is-the-second-derivative-of-a-b-spline
  double bb1 = Bspline(x,j+1,t,p-1);
  double bb2 = Bspline(x,j,t,p-1);

  double ww1 = 0.0;
  if (t(j+p-1) != t(j))
    ww1 = -1.0/(t(j+p-1)-t(j));
  double ww2 = 0.0;
  if (t(j+p-2) != t(j-1))
    ww2 = 1.0/(t(j+p-2)-t(j-1));

  return (p-1.0) * (ww1 * bb1 + ww2 * bb2);
}

double Bspline2nd(double x, int j, const Eigen::VectorXd& t,int p) {
  // 1st derivative of the jth B-spline
  // https://stats.stackexchange.com/questions/258786/what-is-the-second-derivative-of-a-b-spline
  double bb1 = Bspline1st(x,j+1,t,p-1);
  double bb2 = Bspline1st(x,j,t,p-1);

  double ww1 = 0.0;
  if (t(j+p-1) != t(j))
    ww1 = -1.0/(t(j+p-1)-t(j));
  double ww2 = 0.0;
  if (t(j+p-2) != t(j-1))
    ww2 = 1.0/(t(j+p-2)-t(j-1));

  return (p-1.0) * (ww1 * bb1 + ww2 * bb2);
}

Eigen::VectorXd Bsplinevec(double x, const Eigen::VectorXd& t,int p) {
  int m = t.size() - p;
  Eigen::VectorXd b(m);
  b.setZero();
  // int k = knotindexEigen(x,t);
  // for (int i=(k-(p-1));i<k+1;i++)
  for (int i=0;i<m;i++)
    b(i) = Bspline(x,i+1,t,p);
  return b;
}

Eigen::VectorXd BsplinevecCon(double x, const Eigen::VectorXd& t, int p, Eigen::MatrixXd Z) {
  int m = t.size() - p;
  Eigen::VectorXd b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline(x,i+1,t,p);
  return Z.transpose()*b;
  // return b.transpose()*Z;
}

Eigen::VectorXd Bsplinevec1st(double x, const Eigen::VectorXd& t,int p) {
  int m = t.size() - p;
  Eigen::VectorXd b(m);
  b.setZero();
  // int k = knotindexEigen(x,t);
  // for (int i=(k-(p-1));i<k+1;i++)
  for (int i=0;i<m;i++)
    b(i) = Bspline1st(x,i+1,t,p);
  return b;
}

// [[Rcpp::export]]
Eigen::VectorXd BsplinevecCon1st(double x, const Eigen::VectorXd& t, int p, Eigen::MatrixXd Z) {
  int m = t.size() - p;
  Eigen::VectorXd b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline1st(x,i+1,t,p);
  return Z.transpose()*b;
  // return b.transpose()*Z;
}

Eigen::VectorXd Bsplinevec2nd(double x, const Eigen::VectorXd& t,int p) {
  int m = t.size() - p;
  Eigen::VectorXd b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline2nd(x,i+1,t,p);
  return b;
}

// [[Rcpp::export]]
Eigen::VectorXd BsplinevecCon2nd(double x, const Eigen::VectorXd& t,int p, Eigen::MatrixXd Z) {
  int m = t.size() - p;
  Eigen::VectorXd b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline2nd(x,i+1,t,p);
  return Z.transpose()*b;
  // return b.transpose()*Z;
}





// https://github.com/tminka/lightspeed/blob/master/digamma.m
double lgamma1st (double x) {
  const double pi = 3.141592653589793238462643383279;
  const double large = 9.5;
  const double d1 = -0.5772156649015328606065121;
  const double d2 = pi*pi/6.0;
  const double small = 1e-6;
  const double s3 = 1.0/12.0;
  const double s4 = 1.0/120.0;
  const double s5 = 1.0/252.0;
  const double s6 = 1.0/240.0;
  const double s7 = 1.0/132.0;
  const double s8 = 691.0/32760.0;
  const double s9 = 1.0/12.0;
  const double s10 = 3617.0/8160.0;

  // Use de Moivre's expansion if x >= large = 9.5
  // calculate lgamma1st(x+10)
  double xplus10 = x + 10.0;
  double y = 0.0;
  double r = 1.0 / xplus10;
  y += log(xplus10) - 0.5 * r;
  r = r * r;
  y = y - r * ( s3 - r * ( s4 - r * (s5 - r * (s6 - r * s7))));

  // lgamma1st(x+10) = (1/x + 1/(x+1) + ... + 1/(x+9)) + lgamma1st(x)
  y = y - 1.0/x - 1.0/(x+1.0) - 1.0/(x+2.0) - 1.0/(x+3.0) - 1.0/(x+4.0) - 1.0/(x+5.0) - 1.0/(x+6.0) - 1.0/(x+7.0) - 1.0/(x+8.0) - 1/(x+9);

  return y;
}

int containsNaN(const Eigen::MatrixXd& mat) {
    int out = 0;
    for (int i = 0; i < mat.rows(); ++i) {
        for (int j = 0; j < mat.cols(); ++j) {
            if (std::isnan(mat(i, j))) {
              out++;
              // return true;
            }
        }
    }
    return out;
}



// https://github.com/tminka/lightspeed/blob/master/trigamma.m

double lgamma2nd (double x) {
  const double pi = 3.141592653589793238462643383279;
  const double c = pi*pi/6;
  const double c1 = -2.404113806319188570799476;
  const double b2 =  1.0/6.0;
  const double b4 = -1.0/30.0;
  const double b6 =  1.0/42.0;
  const double b8 = -1.0/30.0;
  const double b10 = 5.0/66.0;

  // TO DO: % Reduce to trigamma(x+n) where ( X + N ) >= large.

  double z = 1./(x*x);
  double y = 0.5*z + (1.0 + z*(b2 + z*(b4 + z*(b6 + z*(b8 + z*b10))))) / x;

  // std::cout << "x" << CppAD::Value(x) << std::endl;
  // std::cout << "trigamma" << CppAD::Value(y) << std::endl;
  return y;
}

// **************** PART 2: g(mu) = DL term + linear term + smooth term *************************
class Model {
  // The DLNM model

private:
  // DATA
  const Eigen::VectorXd y; // Response
  const Eigen::MatrixXd Sf; // penalty matrix for f(E)
  // const Eigen::MatrixXd SfR;
  const Eigen::MatrixXd Sw;
  // const Eigen::MatrixXd SwR;
  const std::vector<Eigen::MatrixXd> B_inner_list;
  const Eigen::VectorXd knots_f;
  const Eigen::VectorXd knots_w;


  const Eigen::MatrixXd Xfix; // fixed effects
  const Eigen::MatrixXd Xrand; // random effects
  const Eigen::MatrixXd Zf;

  const Eigen::VectorXd Xoffset; // offset


public:

  // DATA
  int n;
  int kE;
  int kl; // knots for marginal lag
  int kx; // knots for marginal x
  int L; // max Lag
  int kbetaR;
  int kbetaF;
  int mE;
  int p; // number of smooth terms in Xrand

  Eigen::MatrixXd SfR;
  Eigen::MatrixXd SwR;

  Eigen::MatrixXd Blag; // lag basis
  Eigen::VectorXd Lseq; // sequence of lags

  const Eigen::MatrixXd Bindex;

  const Eigen::VectorXd r; // rank of each smooth

  // PARAMETERS
  Eigen::VectorXd alpha_f;
  double log_theta;
  double log_smoothing_f;
  double log_smoothing_w;

  Eigen::VectorXd betaF; // parameters for fixed effects
  Eigen::VectorXd betaR; // parameters for random effects
  Eigen::VectorXd con_index_par; // unconstrained index parameter
  Eigen::VectorXd logsmoothing; // log smoothing parameters for random effects

  // Components generated
  Eigen::MatrixXd B_inner; // B_inner = sum_m gamma_m B_inner_list.at(m)
  double theta;
  double smoothing_f;
  double smoothing_w;
  Eigen::VectorXd smoothing;
  Eigen::VectorXd con_index_par_long;
  double index_par_denominator;
  Eigen::VectorXd index_par; // index parameter = con_index_par_long/index_par_denominator
  Eigen::MatrixXd Bf_matrix;
  Eigen::VectorXd Bf;
  Eigen::VectorXd Bf_tmp;
  Eigen::VectorXd eta;
  Eigen::VectorXd eta_remaining; // remaining terms = Xfix * betaF + Xrand * betaR
  Eigen::VectorXd mu; // log(mu) = eta + eta_remaining + Xoffset
  double NegLogL; // NegativeLogLikelihood value
  Eigen::MatrixXd BtBindex; // BtBindex = Bindex.transpose() * Bindex


  // Components for derivatives

  Eigen::MatrixXd dlogmu_df_mat;
  Eigen::MatrixXd dlogmu_dbetaR_mat;
  Eigen::MatrixXd dlogmu_dbetaF_mat;
  Eigen::MatrixXd dlogmu_dindex_par_mat;
  Eigen::VectorXd dlogdensity_dmu_vec;
  Eigen::MatrixXd dmu_df_mat;
  Eigen::MatrixXd dmu_dbetaR_mat;
  Eigen::MatrixXd dmu_dbetaF_mat;
  Eigen::MatrixXd dmu_dindex_par_mat;
  Eigen::MatrixXd dindex_par_dcon_index_par_mat;
  Eigen::VectorXd gr_index_par_vec;
  Eigen::VectorXd d2logdensity_dmudmu_vec;
  // std::vector<Eigen::MatrixXd> d2mu_dfdf_list;
  // std::vector<Eigen::MatrixXd> d2mu_dbetaRdbetaR_list;
  // std::vector<Eigen::MatrixXd> d2mu_dbetaFdbetaF_list;
  std::vector<Eigen::MatrixXd> d2index_par_dcon_index_pardcon_index_par_list;
  Eigen::MatrixXd he_index_par_mat;
  Eigen::MatrixXd he_alpha_f_index_par_mat;
  double dlogdensity_dtheta_scalar;
  double d2logdensity_dthetadtheta_scalar;
  Eigen::VectorXd d2logdensity_dmudtheta_vec;


  // gradient and hessian for updating alpha_f, betaR and betaF
  Eigen::VectorXd gr_inner_vec;
  Eigen::MatrixXd he_inner_mat;

  // full gradient
  Eigen::VectorXd gr_alpha_f_vec;
  Eigen::VectorXd gr_betaR_vec;
  Eigen::VectorXd gr_betaF_vec;
  Eigen::VectorXd gr_con_index_par_vec;
  double gr_log_smoothing_f_scalar;
  double gr_log_smoothing_w_scalar;
  double gr_log_theta_scalar;
  Eigen::VectorXd gr_logsmoothing_vec;

  Eigen::VectorXd gr_s_u_vec;
  Eigen::VectorXd gr_s_par_vec;

  // full hessian
  Eigen::MatrixXd he_alpha_f_mat;
  Eigen::MatrixXd he_betaR_mat;
  Eigen::MatrixXd he_betaF_mat;
  Eigen::MatrixXd he_con_index_par_mat;
  Eigen::MatrixXd he_alpha_f_con_index_par_mat;
  Eigen::MatrixXd he_alpha_f_betaF_mat;
  Eigen::MatrixXd he_alpha_f_betaR_mat;
  Eigen::MatrixXd he_betaR_betaF_mat;
  Eigen::MatrixXd he_con_index_par_betaF_mat;
  Eigen::MatrixXd he_con_index_par_betaR_mat;
  double he_log_smoothing_f_scalar;
  double he_log_smoothing_w_scalar;
  Eigen::MatrixXd he_logsmoothing_mat;
  double he_log_theta_scalar;
  Eigen::VectorXd he_alpha_f_log_smoothing_f_vec;
  Eigen::VectorXd he_alpha_f_log_smoothing_w_vec;
  Eigen::MatrixXd he_betaR_logsmoothing_mat;
  Eigen::VectorXd he_alpha_f_log_theta_vec;
  Eigen::VectorXd he_betaR_log_theta_vec;
  Eigen::VectorXd he_betaF_log_theta_vec;
  Eigen::VectorXd he_con_index_par_log_theta_vec;

  Eigen::MatrixXd he_s_u_mat;
  Eigen::MatrixXd he_s_par_u_mat;

  // To compute AIC
  double NegLogL_l; // NegativeLogLikelihood without penalty
  // matrix for I (hessian of log likelihood without penalty)
  Eigen::MatrixXd I_alpha_f_mat;
  Eigen::MatrixXd I_betaR_mat;
  Eigen::MatrixXd I_mat;
  Eigen::MatrixXd IS_mat; // I_mat with smoothing penalty
  Eigen::MatrixXd he_s_par_u_log_theta_mat;
  Eigen::MatrixXd I_log_smoothing_mat;

  Eigen::MatrixXd K_alpha_f_mat;
  Eigen::MatrixXd K_betaR_mat;
  Eigen::MatrixXd K_betaF_mat;
  Eigen::MatrixXd K_index_par_mat;
  Eigen::MatrixXd K_con_index_par_mat;
  Eigen::VectorXd K_log_theta_vec;

  Eigen::MatrixXd Kleft;
  Eigen::MatrixXd Khat; // Khat = Kleft * Kleft.transpose()


  // NCV
  Eigen::MatrixXd I_alpha_f_i_mat;
  Eigen::MatrixXd I_betaR_i_mat;
  Eigen::MatrixXd he_betaF_i_mat;
  Eigen::MatrixXd he_index_par_i_mat;
  Eigen::MatrixXd he_con_index_par_i_mat;
  Eigen::MatrixXd he_alpha_f_index_par_i_mat;
  Eigen::MatrixXd he_alpha_f_con_index_par_i_mat;
  Eigen::MatrixXd he_alpha_f_betaF_i_mat;
  Eigen::MatrixXd he_alpha_f_betaR_i_mat;
  Eigen::MatrixXd he_con_index_par_betaF_i_mat;
  Eigen::MatrixXd he_betaR_betaF_i_mat;
  Eigen::MatrixXd he_con_index_par_betaR_i_mat;

  Eigen::VectorXd he_alpha_f_log_theta_i_vec;
  Eigen::VectorXd he_betaR_log_theta_i_vec;
  Eigen::VectorXd he_betaF_log_theta_i_vec;
  Eigen::VectorXd he_con_index_par_log_theta_i_vec;


  // results for profile likelihood
  Eigen::VectorXd PL_gradient;
  Eigen::MatrixXd PL_hessian;
  int converge; // 0: converge. 99: not converge

  double PLg; // modelobj.PLg = g.maxCoeff() in inner()

  double logdetSplus; // log of product of positive eigenvalues of penalty matrix.
  Eigen::MatrixXd Sinv; // inverse of penalty matrix
  Eigen::VectorXd eigvalS; // eigenvalues of penalty matrix
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigS;

  // Constructor
  Model(const Eigen::VectorXd& y_,
        const std::vector<Eigen::MatrixXd>& B_inner_list_,
        const Eigen::VectorXd& knots_f_,
        const Eigen::VectorXd& knots_w_,
        const Eigen::MatrixXd& Sw_,
        const Eigen::MatrixXd& SwR_,
        const Eigen::MatrixXd& Sf_,
        const Eigen::MatrixXd& SfR_,
        const Eigen::MatrixXd& Xrand_,
        const Eigen::MatrixXd& Xfix_,
        const Eigen::MatrixXd& Zf_,
        const Eigen::VectorXd& Xoffset_,
        const Eigen::VectorXd& r_,
        const Eigen::MatrixXd& Bindex_,
        Eigen::VectorXd& alpha_f_,
        Eigen::VectorXd& con_index_par_,
        double log_theta_,
        double log_smoothing_f_,
        double log_smoothing_w_,
        Eigen::VectorXd& betaR_,
        Eigen::VectorXd& betaF_,
        Eigen::VectorXd& logsmoothing_) :
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
      index_par = Bindex * con_index_par_long /index_par_denominator;


      for (int i = 0; i < mE; i++) {
        B_inner += B_inner_list.at(i) * index_par(i);
      }




      Lseq.resize(L+1); // Set the length of Lseq to L+1
      for (int l = 0; l < (L+1); l++) {
        Lseq(l) = (double) l;
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
          Eigen::VectorXd bx = BsplinevecCon(B_inner(i, j), knots_f, 4, Zf);
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


      gr_s_u_vec.resize(kE+mE-1+kbetaR+kbetaF);
      he_s_u_mat.resize(kE+mE-1+kbetaR+kbetaF, kE+mE-1+kbetaR+kbetaF);
      I_mat.resize(kE+mE-1+kbetaR+kbetaF + 1, kE+mE-1+kbetaR+kbetaF + 1);
      IS_mat.resize(kE+mE-1+kbetaR+kbetaF + 1, kE+mE-1+kbetaR+kbetaF + 1);
      I_log_smoothing_mat.resize(2+p, 2+p);
      gr_s_par_vec.resize(3+p);
      he_s_par_u_mat.resize(3+p, kE+mE-1+kbetaR+kbetaF);
      he_s_par_u_log_theta_mat.resize(3+p-1, kE+mE-1+kbetaR+kbetaF + 1);

      gr_inner_vec.resize(kE+kbetaR+kbetaF);
      he_inner_mat.resize(kE+kbetaR+kbetaF, kE+kbetaR+kbetaF);





      // calculate logdetSplus
      // eigS.compute(smoothing_f*SfR + smoothing_w*SwR, true);
      eigS.compute(smoothing_f*SfR + smoothing_w*SwR, Eigen::EigenvaluesOnly);
      // std::cout << "smoothing_f: " << smoothing_f << std::endl;
      // std::cout << "smoothing_w: " << smoothing_w << std::endl;
      logdetSplus = 0.0;
      eigvalS = eigS.eigenvalues();
      // std::cout << "eigvalS: " << eigvalS.transpose() << std::endl;
      for (int i = 0; i < eigvalS.size(); i++) {
        if (eigvalS(i) > 1e-8) logdetSplus += log(eigvalS(i));
      }
      // std::cout << "logdetSplus: " << logdetSplus << std::endl;
      // Sinv = eigS.eigenvectors() * eigvalS.cwiseInverse().asDiagonal() * (eigS.eigenvectors().transpose());
      // std::cout << "eigS.eigenvectors()*eigS.eigenvectors().transpose(): " << std::endl;
      // std::cout << eigS.eigenvectors().transpose() * eigS.eigenvectors() << std::endl;
      // Sinv = (smoothing_f*SfR + smoothing_w*SwR).inverse();
      Sinv = (smoothing_f*SfR + smoothing_w*SwR).ldlt().solve(Eigen::MatrixXd::Identity(SfR.rows(), SfR.cols()));

      derivative_coef();
      derivative_he();
      derivative_full();

      NegativeLogLikelihood();

      // Initialize PL
      PL_gradient.resize(mE-1);
      PL_hessian.resize(mE-1, mE-1);

    }


  // Functions to set parameters
  void setAlphaF(const Eigen::VectorXd alpha_f_) {
    alpha_f = alpha_f_;

    for (int i = 0; i < n; i++) {
      eta(i) = Bf_matrix.row(i).dot(alpha_f);
      mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
    }
  }

  void setAlphaFBetaRBetaFCon_Index_Par(const Eigen::VectorXd alpha_f_, const Eigen::VectorXd betaR_, const Eigen::VectorXd betaF_, const Eigen::VectorXd con_index_par_, const double log_theta_) {
    log_theta = log_theta_;
    theta = exp(log_theta);
    
    con_index_par = con_index_par_;
    for (int i = 0; i < (mE - 1); i++) {
      con_index_par_long(i + 1) = con_index_par(i);
    }
    index_par_denominator = sqrt(con_index_par_long.dot(BtBindex * con_index_par_long));
    index_par = Bindex * con_index_par_long /index_par_denominator;
    B_inner.setZero();
    for (int i = 0; i < mE; i++) {
      B_inner += B_inner_list.at(i) * index_par(i);
    }

    for (int i = 0; i < n; i++) {
      Bf.setZero();
      Bf_tmp.setZero();
      for (int j = 0; j < (L+1); j++) {
        Eigen::VectorXd bx = BsplinevecCon(B_inner(i, j), knots_f, 4, Zf);
        for (int ii = 0; ii < kx; ii++) {
          Bf_tmp.segment(ii*kl, kl) = bx(ii) * Blag.col(j);
        }
        Bf += Bf_tmp;
      }
      Bf_matrix.row(i) = Bf;
      eta(i) = Bf.dot(alpha_f);
      eta_remaining(i) = Xfix.row(i).dot(betaF) + Xrand.row(i).dot(betaR);
      mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
    }

    dindex_par_dcon_index_par_mat = dindex_par_dcon_index_par();
    d2index_par_dcon_index_pardcon_index_par_list = d2index_par_dcon_index_pardcon_index_par();
  }


  void setCon_Index_Par(const Eigen::VectorXd con_index_par_){
    con_index_par = con_index_par_;
    for (int i = 0; i < (mE - 1); i++) {
      con_index_par_long(i + 1) = con_index_par(i);
    }
    index_par_denominator = sqrt(con_index_par_long.dot(BtBindex * con_index_par_long));
    index_par = Bindex * con_index_par_long /index_par_denominator;
    B_inner.setZero();
    for (int i = 0; i < mE; i++) {
      B_inner += B_inner_list.at(i) * index_par(i);
    }

    for (int i = 0; i < n; i++) {
      Bf.setZero();
      Bf_tmp.setZero();
      for (int j = 0; j < (L+1); j++) {
        Eigen::VectorXd bx = BsplinevecCon(B_inner(i, j), knots_f, 4, Zf);
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

    dindex_par_dcon_index_par_mat = dindex_par_dcon_index_par();
    d2index_par_dcon_index_pardcon_index_par_list = d2index_par_dcon_index_pardcon_index_par();

  }

  void setBetaF(const Eigen::VectorXd betaF_) {
    betaF = betaF_;
    for (int i = 0; i < n; i++) {
      eta_remaining(i) = Xfix.row(i).dot(betaF) + Xrand.row(i).dot(betaR);
      mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
    }
  }
  void setBetaR(const Eigen::VectorXd betaR_) {
    betaR = betaR_;
    for (int i = 0; i < n; i++) {
      eta_remaining(i) = Xfix.row(i).dot(betaF) + Xrand.row(i).dot(betaR);
      mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
    }
  }
  void setLogTheta(const double log_theta_) {
    log_theta = log_theta_;
    theta = exp(log_theta);
  }

  void setLogSmoothingF(const double log_smoothing_f_) {
    log_smoothing_f = log_smoothing_f_;
    smoothing_f = exp(log_smoothing_f);
    // calculate logdetSplus
    eigS.compute(smoothing_f*SfR + smoothing_w*SwR, Eigen::EigenvaluesOnly);
    logdetSplus = 0.0;
    eigvalS = eigS.eigenvalues();
    for (int i = 0; i < eigvalS.size(); i++) {
      if (eigvalS(i) > 1e-8) logdetSplus += log(eigvalS(i));
    }
    // Sinv = eigS.eigenvectors() * eigvalS.asDiagonal().inverse() * (eigS.eigenvectors().transpose());
    Sinv = (smoothing_f*SfR + smoothing_w*SwR).ldlt().solve(Eigen::MatrixXd::Identity(SfR.rows(), SfR.cols()));
    // std::cout << "setf: logdetSplus: " << logdetSplus << std::endl;
  }

  void setLogSmoothingW(const double log_smoothing_w_) {
    log_smoothing_w = log_smoothing_w_;
    smoothing_w = exp(log_smoothing_w);
    // calculate logdetSplus
    eigS.compute(smoothing_f*SfR + smoothing_w*SwR, Eigen::EigenvaluesOnly);
    logdetSplus = 0.0;
    eigvalS = eigS.eigenvalues();
    for (int i = 0; i < eigvalS.size(); i++) {
      if (eigvalS(i) > 1e-8) logdetSplus += log(eigvalS(i));
    }
    // Sinv = eigS.eigenvectors() * eigvalS.asDiagonal().inverse() * (eigS.eigenvectors().transpose());
    Sinv = (smoothing_f*SfR + smoothing_w*SwR).ldlt().solve(Eigen::MatrixXd::Identity(SfR.rows(), SfR.cols()));
    // std::cout << "setw: logdetSplus: " << logdetSplus << std::endl;
  }
  void setLogsmoothing(const Eigen::VectorXd logsmoothing_) { // log smoothing parameters for remaining terms
    logsmoothing = logsmoothing_;
    for (int i = 0; i < p; i++) smoothing(i) = exp(logsmoothing(i));
  }


  double get_p_i (Eigen::VectorXd alpha_f_i, Eigen::VectorXd con_index_par_i,
                  Eigen::VectorXd betaR_i, Eigen::VectorXd betaF_i, double log_theta_i,
                  int i) {
      Eigen::VectorXd con_index_par_long_i(mE); // con_index_par_long = c(1,con_index_par)
      con_index_par_long_i(0) = 1.0;
      for (int i = 0; i < (mE - 1); i++) {
        con_index_par_long_i(i + 1) = con_index_par_i(i);
      }

      Eigen::VectorXd index_par_i = Bindex * con_index_par_long_i / sqrt(con_index_par_long_i.dot(BtBindex * con_index_par_long_i));

      Eigen::VectorXd B_inner_i(L+1); // B_inner.resize(n, L+1);
      B_inner_i.setZero();
      for (int j = 0; j < mE; j++) {
        B_inner_i += B_inner_list.at(j).row(i).transpose() * index_par_i(j);
      }



      Eigen::VectorXd Bf_i(kE);
      Eigen::VectorXd Bf_tmp_i(kE);
      Bf_i.setZero();
      Bf_tmp_i.setZero();
      for (int j = 0; j < (L+1); j++) {
        Eigen::VectorXd bx = BsplinevecCon(B_inner_i(j), knots_f, 4, Zf);
        for (int ii = 0; ii < kx; ii++) {
          Bf_tmp_i.segment(ii*kl, kl) = bx(ii) * Blag.col(j);
        }
        Bf_i += Bf_tmp_i;
      }
        // Bf_matrix.row(i) = Bf;
      double eta_i = Bf_i.dot(alpha_f_i);
      double eta_remaining_i = Xfix.row(i).dot(betaF_i) + Xrand.row(i).dot(betaR_i);
      double mu_i = exp(eta_i + eta_remaining_i + Xoffset(i));

      double theta_i = exp(log_theta_i);
      double log_p_i = lgamma(y(i) + theta_i) - lgamma(theta_i) - lgamma(y(i) + 1) -
                                    theta_i * log(1 + mu_i/theta_i) +
                                    y(i)*( eta_i + eta_remaining_i + Xoffset(i) - log_theta_i - log(1 + mu_i/theta_i) );
      return exp(log_p_i);

  }


  double get_D_i (Eigen::VectorXd alpha_f_i, Eigen::VectorXd con_index_par_i,
                  Eigen::VectorXd betaR_i, Eigen::VectorXd betaF_i, double log_theta_i,
                  int i) {
      Eigen::VectorXd con_index_par_long_i(mE); // con_index_par_long = c(1,con_index_par)
      con_index_par_long_i(0) = 1.0;
      for (int i = 0; i < (mE - 1); i++) {
        con_index_par_long_i(i + 1) = con_index_par_i(i);
      }

      Eigen::VectorXd index_par_i = Bindex * con_index_par_long_i / sqrt(con_index_par_long_i.dot(BtBindex * con_index_par_long_i));

      Eigen::VectorXd B_inner_i(L+1); // B_inner.resize(n, L+1);
      B_inner_i.setZero();
      for (int j = 0; j < mE; j++) {
        B_inner_i += B_inner_list.at(j).row(i).transpose() * index_par_i(j);
      }



      Eigen::VectorXd Bf_i(kE);
      Eigen::VectorXd Bf_tmp_i(kE);
      Bf_i.setZero();
      Bf_tmp_i.setZero();
      for (int j = 0; j < (L+1); j++) {
        Eigen::VectorXd bx = BsplinevecCon(B_inner_i(j), knots_f, 4, Zf);
        for (int ii = 0; ii < kx; ii++) {
          Bf_tmp_i.segment(ii*kl, kl) = bx(ii) * Blag.col(j);
        }
        Bf_i += Bf_tmp_i;
      }
        // Bf_matrix.row(i) = Bf;
      double eta_i = Bf_i.dot(alpha_f_i);
      double eta_remaining_i = Xfix.row(i).dot(betaF_i) + Xrand.row(i).dot(betaR_i);
      double mu_i = exp(eta_i + eta_remaining_i + Xoffset(i));


      double theta_i = exp(log_theta_i);
      double log_p_i = lgamma(y(i) + theta_i) - lgamma(theta_i) - lgamma(y(i) + 1) -
                                      theta_i * log(1 + mu_i/theta_i) +
                                      y(i)*( eta_i + eta_remaining_i + Xoffset(i) - log_theta_i - log(1 + mu_i/theta_i) );
      return -1.0*log_p_i; // loss function negative log likelihood
  }

  // get private members
  std::vector<Eigen::MatrixXd> getB_inner_list () {
    return B_inner_list;
  }
  Eigen::MatrixXd getB_inner () {
    return B_inner;
  }


  Eigen::VectorXd getknots_f () {
    return knots_f;
  }
  Eigen::VectorXd getknots_w () {
    return knots_w;
  }
  Eigen::MatrixXd getZf () {
    return Zf;
  }
  Eigen::MatrixXd getXfix () {
    return Xfix;
  }
  Eigen::MatrixXd getXrand () {
    return Xrand;
  }
  Eigen::VectorXd getXoffset() {
    return Xoffset;
  }
  Eigen::MatrixXd getSf () {
    return Sf;
  }
  Eigen::MatrixXd getSw () {
    return Sw;
  }
  // Function to update derivatives.
  // RUN the function derivative_coef(), derivative_he() and derivative_full() after update parameters.
  // update derivatives related to spline coefficients alpha_f, and betaR and betaF
  void derivative_coef() {
    dlogmu_dindex_par_mat = dlogmu_dindex_par();
    dlogmu_df_mat = dlogmu_df();
    dlogmu_dbetaR_mat = dlogmu_dbetaR();
    dlogmu_dbetaF_mat = dlogmu_dbetaF();
    dlogdensity_dmu_vec = dlogdensity_dmu();
    dmu_df_mat = dmu_df();
    dmu_dindex_par_mat = dmu_dindex_par();

    dmu_dbetaR_mat = dmu_dbetaR();
    dmu_dbetaF_mat = dmu_dbetaF();
    d2logdensity_dmudmu_vec = d2logdensity_dmudmu();


    he_index_par_mat = he_index_par();

    he_alpha_f_index_par_mat = he_alpha_f_index_par();

    dlogdensity_dtheta_scalar = dlogdensity_dtheta();
    d2logdensity_dmudtheta_vec = d2logdensity_dmudtheta();

    // obtain gradient
    gr_alpha_f_vec = gr_alpha_f();
    gr_index_par_vec = gr_index_par();
    gr_betaR_vec = gr_betaR();
    gr_betaF_vec = gr_betaF();
    gr_con_index_par_vec = gr_con_index_par();
    // obtain hessian
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
  }

  // update full gradient and hessian of alpha_f, con_index_par and betaR and betaF
  void derivative_he () {
    gr_s_u_vec << gr_alpha_f_vec, gr_con_index_par_vec, gr_betaR_vec, gr_betaF_vec;


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

    // make it symmetric. Comment out ...
    // he_s_u_mat = (he_s_u_mat + he_s_u_mat.transpose())/2.0;
  }

  // update derivatives related to overdispersion and smoothing parameters
  // Full derivative for LAML
  void derivative_full () {
    // obtain full gradient

    gr_log_smoothing_f_scalar = gr_log_smoothing_f();
    gr_log_smoothing_w_scalar = gr_log_smoothing_w();
    gr_log_theta_scalar = gr_log_theta();
    gr_logsmoothing_vec = gr_logsmoothing();



    // u represents spline coefficient alpha_f, and betaR and betaF
    // par represents overdispersion and smoothing parameters

    gr_s_par_vec << gr_log_theta_scalar, gr_log_smoothing_f_scalar, gr_log_smoothing_w_scalar, gr_logsmoothing_vec;


    // obtain full hessian
    he_alpha_f_log_smoothing_f_vec = he_alpha_f_log_smoothing_f();
    he_alpha_f_log_smoothing_w_vec = he_alpha_f_log_smoothing_w();
    he_betaR_logsmoothing_mat = he_betaR_logsmoothing();
    he_alpha_f_log_theta_vec = he_alpha_f_log_theta();
    he_betaR_log_theta_vec = he_betaR_log_theta();
    he_con_index_par_log_theta_vec = he_con_index_par_log_theta();
    he_betaF_log_theta_vec = he_betaF_log_theta();


    he_s_par_u_mat.setZero();

    he_s_par_u_mat.row(0) << he_alpha_f_log_theta_vec.transpose(), he_con_index_par_log_theta_vec.transpose(), he_betaR_log_theta_vec.transpose(), he_betaF_log_theta_vec.transpose();
    he_s_par_u_mat.block(1, 0, 1, kE) = he_alpha_f_log_smoothing_f_vec.transpose();
    he_s_par_u_mat.block(2, 0, 1, kE) = he_alpha_f_log_smoothing_w_vec.transpose();
    he_s_par_u_mat.block(3, kE+mE-1, p, kbetaR) = he_betaR_logsmoothing_mat.transpose();
  }

  // update variables related to alpha_f, betaR and betaF.
  // Used only in updating alpha_f, betaR and betaF. .
  void derivative_f () {
    dlogmu_df_mat = dlogmu_df();
    dlogmu_dbetaR_mat = dlogmu_dbetaR();
    dlogmu_dbetaF_mat = dlogmu_dbetaF();
    dlogdensity_dmu_vec = dlogdensity_dmu();
    dmu_df_mat = dmu_df();
    dmu_dbetaR_mat = dmu_dbetaR();
    dmu_dbetaF_mat = dmu_dbetaF();
    d2logdensity_dmudmu_vec = d2logdensity_dmudmu();

    gr_alpha_f_vec = gr_alpha_f();
    gr_betaR_vec = gr_betaR();
    gr_betaF_vec = gr_betaF();

    he_alpha_f_mat = he_alpha_f();
    he_betaR_mat = he_betaR();
    he_betaF_mat = he_betaF();
    he_alpha_f_betaF_mat = he_alpha_f_betaF();
    he_alpha_f_betaR_mat = he_alpha_f_betaR();
    he_betaR_betaF_mat = he_betaR_betaF();

    gr_inner_vec << gr_alpha_f_vec, gr_betaR_vec, gr_betaF_vec;

    he_inner_mat.setZero();
    he_inner_mat.block(0,0,kE,kE) = he_alpha_f_mat;
    he_inner_mat.block(kE,kE,kbetaR,kbetaR) = he_betaR_mat;
    he_inner_mat.block(kE+kbetaR,kE+kbetaR,kbetaF,kbetaF) = he_betaF_mat;

    he_inner_mat.block(0,kE,kE,kbetaR) = he_alpha_f_betaR_mat;
    he_inner_mat.block(0,kE+kbetaR,kE,kbetaF) = he_alpha_f_betaF_mat;
    he_inner_mat.block(kE,kE+kbetaR,kbetaR,kbetaF) = he_betaR_betaF_mat;

    he_inner_mat.block(kE,0,kbetaR,kE) = he_alpha_f_betaR_mat.transpose();
    he_inner_mat.block(kE+kbetaR,0,kbetaF,kE) = he_alpha_f_betaF_mat.transpose();
    he_inner_mat.block(kE+kbetaR,kE,kbetaF,kbetaR) = he_betaR_betaF_mat.transpose();
  }


  // functions for NegativeLogLikelihood
  void NegativeLogLikelihood() {

    double loglik = 0;
    for (int i = 0; i < n; i++) {
      loglik += lgamma(y(i) + theta) - lgamma(theta) - lgamma(y(i) + 1) -
                                    theta * log(1 + mu(i)/theta) +
                                    y(i)*( eta(i) + eta_remaining(i) + Xoffset(i) - log_theta - log(1 + mu(i)/theta) );
    }
    // part 1: DLNM
    // Smooth Penalty
    loglik += -0.5 * smoothing_w * alpha_f.dot(Sw * alpha_f) - 0.5 * smoothing_f * alpha_f.dot(Sf * alpha_f);
    // Scale
    // loglik += ((double)(kx*(kl-2)) / 2.0) * log_smoothing_w + ((double)((kx-1)*kl) / 2.0) * log_smoothing_f;
    loglik += logdetSplus/2.0;

    // part 2: Remaining smooth terms
    int begin = 0;
    for (int i = 0; i < p; i++) {
      // Smooth Penalty
      int ki = static_cast<int>(r(i));
      Eigen::VectorXd betaRi(ki);
      for (int j = 0; j < ki; j++) betaRi(j) = betaR(begin + j);
      loglik += -0.5 * smoothing(i) * betaRi.dot(betaRi); // smooth penalty
      loglik += ki/2.0 * logsmoothing(i); // scale

      begin += ki;
    }


    NegLogL = -1.0 * loglik; // NEGATIVE log-likelihood
  }

  // functions for NegativeLogLikelihood WITHOUT penalty for AIC
  void NegativeLogLikelihood_l() {

    double loglik = 0;
    for (int i = 0; i < n; i++) {
      loglik += lgamma(y(i) + theta) - lgamma(theta) - lgamma(y(i) + 1) -
                                    theta * log(1 + mu(i)/theta) +
                                    y(i)*( eta(i) + eta_remaining(i) + Xoffset(i) - log_theta - log(1 + mu(i)/theta) );
    }

    NegLogL_l = -1.0 * loglik; // NEGATIVE log-likelihood
  }

  void prepare_AIC () {
    NegativeLogLikelihood_l();
    // hessian of log likelihood without penalty
    I_alpha_f_mat = I_alpha_f();
    I_betaR_mat = I_betaR();

    I_mat.block(0,0,kE+mE-1+kbetaR+kbetaF, kE+mE-1+kbetaR+kbetaF) = he_s_u_mat;
    I_mat.block(0, 0, kE, kE)  = I_alpha_f_mat;
    I_mat.block(kE+mE-1, kE+mE-1, kbetaR, kbetaR) = I_betaR_mat;


    d2logdensity_dthetadtheta_scalar = d2logdensity_dthetadtheta();
    he_log_smoothing_f_scalar = he_log_smoothing_f();
    he_log_smoothing_w_scalar = he_log_smoothing_w();
    he_logsmoothing_mat = he_logsmoothing();
    he_log_theta_scalar = he_log_theta();



    I_mat.row(kE+mE-1+kbetaR+kbetaF) << he_s_par_u_mat.row(0), he_log_theta_scalar;
    I_mat.col(kE+mE-1+kbetaR+kbetaF) = I_mat.row(kE+mE-1+kbetaR+kbetaF).transpose();

    IS_mat.block(0,0,kE+mE-1+kbetaR+kbetaF, kE+mE-1+kbetaR+kbetaF) = he_s_u_mat;
    IS_mat.row(kE+mE-1+kbetaR+kbetaF) << he_s_par_u_mat.row(0), he_log_theta_scalar;
    IS_mat.col(kE+mE-1+kbetaR+kbetaF) = IS_mat.row(kE+mE-1+kbetaR+kbetaF).transpose();


    I_log_smoothing_mat.setZero();
    I_log_smoothing_mat(0,0) = he_log_smoothing_f_scalar;
    I_log_smoothing_mat(1,1) = he_log_smoothing_w_scalar;
    // TODO: he for smoothing_f and smoothing_w?
    I_log_smoothing_mat.block(2,2,p,p) = he_logsmoothing_mat;



    // for d (u, log_theta) / d (par)
    he_s_par_u_log_theta_mat.setZero();
    he_s_par_u_log_theta_mat.block(0, 0, 1, kE) = he_alpha_f_log_smoothing_f_vec.transpose();
    he_s_par_u_log_theta_mat.block(1, 0, 1, kE) = he_alpha_f_log_smoothing_w_vec.transpose();
    he_s_par_u_log_theta_mat.block(2, kE+mE-1, p, kbetaR) = he_betaR_logsmoothing_mat.transpose();
    // he_s_par_u_log_theta_mat.block(1, kE+mE-1, p, kbetaR) = he_betaR_logsmoothing_mat.transpose();



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

  Eigen::MatrixXd prepare_NCV (Eigen::VectorXd nei_vec) {
    Eigen::MatrixXd Hunpen_nei(kE+mE-1+kbetaR+kbetaF + 1, kE+mE-1+kbetaR+kbetaF + 1);
    Hunpen_nei.setZero();
    Eigen::MatrixXd tmp(kE+mE-1+kbetaR+kbetaF + 1, kE+mE-1+kbetaR+kbetaF + 1);
    tmp.setZero();

    // compute Hunpen_nei
    for (size_t i = 0; i < nei_vec.size(); i++) {
      int index_int = static_cast<int>(nei_vec(i)) - 1; // start from 0 in C++
      I_alpha_f_i_mat = I_alpha_f_i(index_int);
      I_betaR_i_mat = I_betaR_i(index_int);
      he_betaF_i_mat = he_betaF_i(index_int);
      he_index_par_i_mat = he_index_par_i(index_int);
      he_con_index_par_i_mat = he_con_index_par_i(index_int);
      he_alpha_f_index_par_i_mat = he_alpha_f_index_par_i(index_int);
      he_alpha_f_con_index_par_i_mat = he_alpha_f_con_index_par_i(index_int);
      he_alpha_f_betaF_i_mat = he_alpha_f_betaF_i(index_int);
      he_alpha_f_betaR_i_mat = he_alpha_f_betaR_i(index_int);
      he_con_index_par_betaF_i_mat = he_con_index_par_betaF_i(index_int);
      he_betaR_betaF_i_mat = he_betaR_betaF_i(index_int);
      he_con_index_par_betaR_i_mat = he_con_index_par_betaR_i(index_int);

      he_alpha_f_log_theta_i_vec = he_alpha_f_log_theta_i(index_int);
      he_betaR_log_theta_i_vec = he_betaR_log_theta_i(index_int);
      he_betaF_log_theta_i_vec = he_betaF_log_theta_i(index_int);
      he_con_index_par_log_theta_i_vec = he_con_index_par_log_theta_i(index_int);

      tmp.setZero();
      tmp.block(0, 0, kE, kE)  = I_alpha_f_i_mat;

      tmp.block(0, kE, kE, mE-1) = he_alpha_f_con_index_par_i_mat;
      tmp.block(kE, 0, mE-1, kE) = he_alpha_f_con_index_par_i_mat.transpose();
      tmp.block(kE, kE, mE-1, mE-1) = he_con_index_par_i_mat;


      tmp.block(kE+mE-1, kE+mE-1, kbetaR, kbetaR) = I_betaR_i_mat;
      tmp.block(kE+kbetaR+mE-1, kE+kbetaR+mE-1, kbetaF, kbetaF) = he_betaF_i_mat;

      tmp.block(0,kE+mE-1,kE,kbetaR) = he_alpha_f_betaR_i_mat;
      tmp.block(kE,kE+mE-1,mE-1,kbetaR) = he_con_index_par_betaR_i_mat;
      tmp.block(0,kE+mE-1+kbetaR,kE,kbetaF) = he_alpha_f_betaF_i_mat;
      tmp.block(kE,kE+mE-1+kbetaR,mE-1,kbetaF) = he_con_index_par_betaF_i_mat;
      tmp.block(kE+mE-1, kE+mE-1+kbetaR, kbetaR, kbetaF) = he_betaR_betaF_i_mat;

      tmp.block(kE+mE-1,0,kbetaR,kE) = he_alpha_f_betaR_i_mat.transpose();
      tmp.block(kE+mE-1,kE,kbetaR,mE-1) = he_con_index_par_betaR_i_mat.transpose();
      tmp.block(kE+mE-1+kbetaR,0,kbetaF,kE) = he_alpha_f_betaF_i_mat.transpose();
      tmp.block(kE+mE-1+kbetaR,kE,kbetaF,mE-1) = he_con_index_par_betaF_i_mat.transpose();
      tmp.block(kE+mE-1+kbetaR, kE+mE-1, kbetaF, kbetaR) = he_betaR_betaF_i_mat.transpose();
      
      tmp.row(kE+mE-1+kbetaR+kbetaF) << he_alpha_f_log_theta_i_vec.transpose(), he_con_index_par_log_theta_i_vec.transpose(), he_betaR_log_theta_i_vec.transpose(), he_betaF_log_theta_i_vec.transpose(), he_log_theta_i(index_int);
      tmp.col(kE+mE-1+kbetaR+kbetaF) = tmp.row(kE+mE-1+kbetaR+kbetaF).transpose();

      Hunpen_nei += tmp;
    }

    return Hunpen_nei;

  }

  // ********* Derivatives *************

  // FUNCTIONS
  // 1. density function
  // d log(exponential family density) / d mu
  Eigen::VectorXd dlogdensity_dmu () {
    Eigen::VectorXd out(n);
    for (int i = 0; i < n; i++) {
      out(i) = y(i) / mu(i) - (theta + y(i)) / (theta + mu(i));
    }
    return out;
  }
  // d^2 log(exponential family density) / d mu^2
  Eigen::VectorXd d2logdensity_dmudmu () {
    Eigen::VectorXd out(n);
    for (int i = 0; i < n; i++) {
      out(i) = - y(i) / pow(mu(i), 2) + (theta + y(i)) / pow(theta + mu(i), 2);
    }
    return out;
  }
  // d log(exponential family density) / d theta
  double dlogdensity_dtheta () {
    double out = 0.0;
    // std::cout << "x" << 3.5 << std::endl;

    // TO DO: optimize it. Use property of gamma function...
    for (int i = 0; i < n; i++) {
      out += log_theta - log(theta + mu(i)) + (mu(i) - y(i))/(theta+mu(i)) + lgamma1st(theta+y(i)) - lgamma1st(theta);
    }
    return out;
  }
    // d^2 log(exponential family density) / d theta^2
  double d2logdensity_dthetadtheta () {
    double out = 0.0;

    for (int i = 0; i < n; i++) {
      out += 1/theta - 1/(theta + mu(i)) - (mu(i) - y(i)) / ((theta + mu(i))*(theta + mu(i))) + lgamma2nd(y(i) + theta) - lgamma2nd(theta);
    }
    return out;
  }

  Eigen::VectorXd d2logdensity_dmudtheta () {
    Eigen::VectorXd out(n);
    for (int i = 0; i < n; i++) {
      out(i) = (y(i) - mu(i)) / pow(theta+mu(i), 2);
    }
    return out;
  }



  // // 2. mean model
  // d log(mu) / d alpha_f
  Eigen::MatrixXd dlogmu_df () {
    return Bf_matrix;
  }
  // d mu / d alpha_f
  Eigen::MatrixXd dmu_df () {
    Eigen::MatrixXd out(n, kE);
    for (int i = 0; i < n; i++) {
      out.row(i) = dlogmu_df_mat.row(i) * mu(i);
    }
    return out;
  }

  // d log mu / d index_par
  Eigen::MatrixXd  dlogmu_dindex_par () {
    Eigen::MatrixXd out(n, mE);
    Eigen::MatrixXd out_tmp(mE, kE);
    Eigen::VectorXd Bf1st_tmp(kE);
    Eigen::VectorXd bx(kx);

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


  // d log(mu) / d betaR
  Eigen::MatrixXd dlogmu_dbetaR () {
    return Xrand;
  }
  // d mu / d betaR
  Eigen::MatrixXd dmu_dbetaR () {
    Eigen::MatrixXd out(n, kbetaR);
    for (int i = 0; i < n; i++) {
      out.row(i) = dlogmu_dbetaR_mat.row(i) * mu(i);
    }
    return out;
  }
  // d mu / d index_par
  Eigen::MatrixXd dmu_dindex_par () {
    Eigen::MatrixXd out(n, mE);
    for (int i = 0; i < n; i++) {
      out.row(i) = dlogmu_dindex_par_mat.row(i) * mu(i);
    }
    return out;
  }
  // d log(mu) / d betaF
  Eigen::MatrixXd dlogmu_dbetaF () {
    return Xfix;
  }
  // d mu / d betaR
  Eigen::MatrixXd dmu_dbetaF () {
    Eigen::MatrixXd out(n, kbetaF);
    for (int i = 0; i < n; i++) {
      out.row(i) = dlogmu_dbetaF_mat.row(i) * mu(i);
    }
    return out;
  }



  // 3. Re-parameterization

  Eigen::MatrixXd dindex_par_dcon_index_par () {
    // index_par_denominator = sqrt(con_index_par_long.dot(con_index_par_long));
    // index_par = con_index_par_long/index_par_denominator;

    Eigen::DiagonalMatrix<double, Eigen::Dynamic> D(mE);
    D.diagonal().setConstant(1.0/index_par_denominator);
    Eigen::MatrixXd Ddense = D.toDenseMatrix();


    Eigen::MatrixXd deriv_g2 = (con_index_par_long * (BtBindex * con_index_par_long).transpose()) * pow(index_par_denominator, -3);

    Eigen::MatrixXd deriv_g = Ddense - deriv_g2;
    // Remove the first column
    Eigen::MatrixXd out = Bindex * (deriv_g.block(0, 1, mE, mE - 1));
    return out;
  }



  std::vector<Eigen::MatrixXd> d2index_par_dcon_index_pardcon_index_par () {
    std::vector<Eigen::MatrixXd> out;
    std::vector<Eigen::MatrixXd> outreal;

    double tmp1 = pow(index_par_denominator, -3); // pow(tmp, -1.5)
    double tmp2 = pow(index_par_denominator, -5); // pow(tmp, -2.5)
    Eigen::VectorXd BtBindexcon_index_par_long = BtBindex * con_index_par_long;
    Eigen::MatrixXd outlarge(mE, mE);


    for (int s = 0; s < mE; s++) {
      if (s == 0) {
        outlarge = -1.0 * tmp1 * BtBindex + 3.0 * tmp2 * BtBindexcon_index_par_long * BtBindexcon_index_par_long.transpose();
      } else {
        Eigen::MatrixXd m1(mE, mE);
        m1.setZero();
        m1.row(s) = BtBindexcon_index_par_long.transpose()*tmp1;
        m1.col(s) = m1.col(s) + BtBindexcon_index_par_long*tmp1;
        Eigen::MatrixXd m2 = -1.0 * tmp1* BtBindex + 3.0 * tmp2 * BtBindexcon_index_par_long * BtBindexcon_index_par_long.transpose();
        outlarge = -1.0 * m1 + m2;
      }

      out.push_back(outlarge.block(1, 1, mE-1, mE-1));
    }

    Eigen::MatrixXd outrealtmp(mE-1, mE-1);
    for (int s = 0; s < mE; s++) {
      outrealtmp.setZero();
      for (int i = 0; i < mE; i++) {
        outrealtmp += Bindex(s, i) * out.at(i);
      }
      outreal.push_back(outrealtmp);
    }

    return outreal;
  }


  // *** GRADIENT ***
  Eigen::VectorXd gr_alpha_f () {
    Eigen::VectorXd out = - dmu_df_mat.transpose() * dlogdensity_dmu_vec + smoothing_f * Sf * alpha_f + smoothing_w * Sw * alpha_f;
    return out;
  }
  Eigen::MatrixXd K_alpha_f () {
    // dmu_df_mat: n * kE
    // dlogdensity_dmu_vec: n * 1
    // out: kE * n
    Eigen::MatrixXd out(kE, n);
    for (int i = 0; i < n; i++) {
      out.col(i) = dmu_df_mat.row(i).transpose() * dlogdensity_dmu_vec(i);
    }
    return out;
  }

  Eigen::VectorXd gr_betaR () {
    Eigen::VectorXd out = - dmu_dbetaR_mat.transpose() * dlogdensity_dmu_vec; // + smoothing * betaR;
    int begin = 0;
    for (int i = 0; i < p; i++) {
      // Smooth Penalty
      int ki = static_cast<int>(r(i));
      // for (int j = 0; j < ki; j++) betaRi(j) = betaR(begin + j);
      // out += smoothing(i) * betaRi;
      for (int j = 0; j < ki; j++) out(begin + j) += smoothing(i) * betaR(begin + j);
      begin += ki;
    }
    return out;
  }

  Eigen::MatrixXd K_betaR () {
    // dmu_dbetaR_mat: n * kbetaR
    // dlogdensity_dmu_vec: n * 1
    // out: kbetaR * n
    Eigen::MatrixXd out(kbetaR, n);
    for (int i = 0; i < n; i++) {
      out.col(i) = dmu_dbetaR_mat.row(i).transpose() * dlogdensity_dmu_vec(i);
    }
    return out;
  }

  Eigen::VectorXd gr_betaF () {
    Eigen::VectorXd out = - dmu_dbetaF_mat.transpose() * dlogdensity_dmu_vec;
    return out;
  }


  Eigen::MatrixXd K_betaF () {
    // dmu_dbetaF_mat: n * kbetaF
    // dlogdensity_dmu_vec: n * 1
    // out: kbetaF * n
    Eigen::MatrixXd out(kbetaF, n);
    for (int i = 0; i < n; i++) {
      out.col(i) = dmu_dbetaF_mat.row(i).transpose() * dlogdensity_dmu_vec(i);
    }
    return out;
  }

  Eigen::VectorXd gr_index_par () {
    Eigen::VectorXd out = - dmu_dindex_par_mat.transpose() * dlogdensity_dmu_vec;
    return out;
  }

  Eigen::MatrixXd K_index_par () {
    // dmu_dindex_par_mat: n * mE
    // dlogdensity_dmu_vec: n * 1
    // out: mE * n
    Eigen::MatrixXd out(mE, n);
    for (int i = 0; i < n; i++) {
      out.col(i) = dmu_dindex_par_mat.row(i).transpose() * dlogdensity_dmu_vec(i);
    }
    return out;
  }


  Eigen::VectorXd gr_con_index_par () {
    Eigen::VectorXd out = dindex_par_dcon_index_par_mat.transpose() * gr_index_par_vec;
    return out;
  }

  Eigen::MatrixXd K_con_index_par () {
    // K_index_par_mat: mE * n
    // dindex_par_dcon_index_par_mat: mE * (mE-1)
    // out: (mE-1) * n
    return dindex_par_dcon_index_par_mat.transpose() * K_index_par_mat;
  }


  double gr_log_smoothing_f () {
    // return 0.5 * smoothing_f * alpha_f.dot(Sf * alpha_f) - 0.5 * ((kx-1)*kl);
    return 0.5 * smoothing_f * alpha_f.dot(Sf * alpha_f) - 0.5 * smoothing_f * (Sinv * SfR).trace();
  }

  double gr_log_smoothing_w () {
    // return 0.5 * smoothing_w * alpha_f.dot(Sw * alpha_f) - 0.5 * (kx*(kl-2));
    return 0.5 * smoothing_w * alpha_f.dot(Sw * alpha_f) - 0.5 * smoothing_w * (Sinv * SwR).trace();
  }
  double gr_log_theta () {
    return -1.0 * theta * dlogdensity_dtheta_scalar;
  }


  Eigen::VectorXd K_log_theta () {
    // out: n * 1. DO NOT FORGET TO TRANSPOSE WHEN USING
    Eigen::VectorXd out(n);
    for (int i = 0; i < n; i++) {
      out(i) = -1.0 * theta * (log_theta - log(theta + mu(i)) + (mu(i) - y(i))/(theta+mu(i)) + lgamma1st(theta+y(i)) - lgamma1st(theta));
    }
    return out;
  }


  Eigen::VectorXd gr_logsmoothing () {
    Eigen::VectorXd out(p);
    int begin = 0;
    for (int i = 0; i < p; i++) {
      // Smooth Penalty
      int ki = static_cast<int>(r(i));
      Eigen::VectorXd betaRi(ki);
      for (int j = 0; j < ki; j++) betaRi(j) = betaR(begin + j);
      out(i) = 0.5 * smoothing(i) * betaRi.dot(betaRi) - 0.5*ki;
      begin += ki;
    }
    return out;
  }


  double he_log_theta () {
    return -1.0 * theta * dlogdensity_dtheta_scalar - theta * theta * d2logdensity_dthetadtheta_scalar;
  }
  double he_log_theta_i (int i) {
    double tmp1 = log_theta - log(theta + mu(i)) + (mu(i) - y(i))/(theta+mu(i)) + lgamma1st(theta+y(i)) - lgamma1st(theta);
    double tmp2 = 1/theta - 1/(theta + mu(i)) - (mu(i) - y(i)) / ((theta + mu(i))*(theta + mu(i))) + lgamma2nd(y(i) + theta) - lgamma2nd(theta);
    return -1.0 * theta * tmp1 - theta * theta * tmp2;
  }

  double he_log_smoothing_f () {
    return 0.5 * smoothing_f * alpha_f.dot(Sf * alpha_f) - 0.5 * smoothing_f * (Sinv * SfR).trace();
  }
  double he_log_smoothing_w () {
    return 0.5 * smoothing_w * alpha_f.dot(Sw * alpha_f) - 0.5 * smoothing_w * (Sinv * SwR).trace();
  }
  Eigen::MatrixXd he_logsmoothing () {
    Eigen::MatrixXd out(p, p);
    out.setZero();
    int begin = 0;
    for (int i = 0; i < p; i++) {
      // Smooth Penalty
      int ki = static_cast<int>(r(i));
      Eigen::VectorXd betaRi(ki);
      for (int j = 0; j < ki; j++) betaRi(j) = betaR(begin + j);
      out(i, i) = 0.5 * smoothing(i) * betaRi.dot(betaRi);
      begin += ki;
    }
    return out;
  }


  // *** Hessian ***
  Eigen::MatrixXd he_alpha_f () {
    Eigen::MatrixXd out1(kE, kE);
    Eigen::MatrixXd out2(kE, kE);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_df_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_df_mat.row(i).transpose() * dlogmu_df_mat.row(i));
    }
    return - out1 - out2 + smoothing_f*Sf + smoothing_w*Sw;
  }

  Eigen::MatrixXd I_alpha_f () { // hessian of negative likelihood without penalty
    Eigen::MatrixXd out1(kE, kE);
    Eigen::MatrixXd out2(kE, kE);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_df_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_df_mat.row(i).transpose() * dlogmu_df_mat.row(i));
    }
    return - out1 - out2;
  }

  Eigen::MatrixXd I_alpha_f_i (int i) {
    Eigen::MatrixXd out1(kE, kE);
    Eigen::MatrixXd out2(kE, kE);

    out1 = d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_df_mat.row(i);
    out2 = dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_df_mat.row(i).transpose() * dlogmu_df_mat.row(i));

    return - out1 - out2;
  }

  Eigen::MatrixXd he_betaR () {
    Eigen::MatrixXd out1(kbetaR, kbetaR);
    Eigen::MatrixXd out2(kbetaR, kbetaR);
    Eigen::MatrixXd Ones(kbetaR, kbetaR); // identity matrix with diagonal smoothing
    out1.setZero();
    out2.setZero();
    Ones.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat.row(i).transpose() * dmu_dbetaR_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dbetaR_mat.row(i).transpose() * dlogmu_dbetaR_mat.row(i));
    }
    int begin = 0;
    for (int i = 0; i < p; i++) {
      // Smooth Penalty
      int ki = static_cast<int>(r(i));
      for (int j = 0; j < ki; j++) Ones(begin + j, begin + j) = smoothing(i);
      begin += ki;
    }
    return - out1 - out2 + Ones;
  }

  Eigen::MatrixXd I_betaR () {  // hessian of negative likelihood without penalty
    Eigen::MatrixXd out1(kbetaR, kbetaR);
    Eigen::MatrixXd out2(kbetaR, kbetaR);
    Eigen::MatrixXd Ones(kbetaR, kbetaR); // identity matrix with diagonal smoothing
    out1.setZero();
    out2.setZero();
    Ones.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat.row(i).transpose() * dmu_dbetaR_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dbetaR_mat.row(i).transpose() * dlogmu_dbetaR_mat.row(i));
    }
    return - out1 - out2;
  }

  Eigen::MatrixXd I_betaR_i (int i) {
    Eigen::MatrixXd out1(kbetaR, kbetaR);
    Eigen::MatrixXd out2(kbetaR, kbetaR);

    out1 = d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat.row(i).transpose() * dmu_dbetaR_mat.row(i);
    out2 = dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dbetaR_mat.row(i).transpose() * dlogmu_dbetaR_mat.row(i));

    return - out1 - out2;
  }

  Eigen::MatrixXd he_betaF () {
    Eigen::MatrixXd out1(kbetaF, kbetaF);
    Eigen::MatrixXd out2(kbetaF, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_dbetaF_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dbetaF_mat.row(i).transpose() * dlogmu_dbetaF_mat.row(i));
    }
    return - out1 - out2;
  }

  Eigen::MatrixXd he_betaF_i (int i) {
    Eigen::MatrixXd out1(kbetaF, kbetaF);
    Eigen::MatrixXd out2(kbetaF, kbetaF);

    out1 = d2logdensity_dmudmu_vec(i) * dmu_dbetaF_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
    out2 = dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dbetaF_mat.row(i).transpose() * dlogmu_dbetaF_mat.row(i));

    return - out1 - out2;
  }



  Eigen::MatrixXd he_index_par () {
    Eigen::MatrixXd out1(mE, mE);
    Eigen::MatrixXd out2(mE, mE);
    Eigen::MatrixXd out_tmp(mE, mE);
    Eigen::MatrixXd xltxl_tmp(mE, mE);
    Eigen::VectorXd bx;
    out1.setZero();
    out2.setZero();

    Eigen::VectorXd Bf2nd_tmp(kE);
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

  Eigen::MatrixXd he_index_par_i (int i) {
    Eigen::MatrixXd out1(mE, mE);
    Eigen::MatrixXd out2(mE, mE);
    Eigen::MatrixXd out_tmp(mE, mE);
    Eigen::MatrixXd xltxl_tmp(mE, mE);
    Eigen::VectorXd bx;

    Eigen::VectorXd Bf2nd_tmp(kE);
    double tmp;

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


    out1 = d2logdensity_dmudmu_vec(i) * dmu_dindex_par_mat.row(i).transpose() * dmu_dindex_par_mat.row(i);
    out2 = dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dindex_par_mat.row(i).transpose() * dlogmu_dindex_par_mat.row(i) + mu(i) * out_tmp);


    return - out1 - out2;
  }

  Eigen::MatrixXd he_con_index_par () {
    Eigen::MatrixXd out1 = dindex_par_dcon_index_par_mat.transpose() * he_index_par_mat * dindex_par_dcon_index_par_mat;
    Eigen::MatrixXd out2(mE-1, mE-1);
    out2.setZero();
    for (int s = 0; s < mE; s++) {
      out2 = out2 + gr_index_par_vec(s) * d2index_par_dcon_index_pardcon_index_par_list.at(s);
    }
    return out1 + out2;
  }


  Eigen::MatrixXd he_con_index_par_i (int i) {
    Eigen::MatrixXd out1 = dindex_par_dcon_index_par_mat.transpose() * he_index_par_i_mat * dindex_par_dcon_index_par_mat;
    Eigen::MatrixXd out2(mE-1, mE-1);
    out2.setZero();
    Eigen::VectorXd tmp = - dmu_dindex_par_mat.row(i).transpose() * dlogdensity_dmu_vec(i);
    for (int s = 0; s < mE; s++) {
      out2 = out2 + tmp(s) * d2index_par_dcon_index_pardcon_index_par_list.at(s);
    }
    return out1 + out2;
  }


  Eigen::MatrixXd he_alpha_f_index_par () {
    Eigen::MatrixXd out1(kE, mE);
    Eigen::MatrixXd out2(kE, mE);
    Eigen::MatrixXd out_mat(kE, mE);
    Eigen::MatrixXd out_tmp(mE, kE); // out_tmp = out_mat.transpose()
    Eigen::VectorXd Bf1st_tmp(kE);
    Eigen::VectorXd bx;
    out1.setZero();
    out2.setZero();
    Eigen::VectorXd Bf1st;
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

      out1 += d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_dindex_par_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_df_mat.row(i).transpose() * dlogmu_dindex_par_mat.row(i) + mu(i)*out_mat);
    }
    return - out1 - out2;
  }

  Eigen::MatrixXd he_alpha_f_index_par_i (int i) {
    Eigen::MatrixXd out1(kE, mE);
    Eigen::MatrixXd out2(kE, mE);
    Eigen::MatrixXd out_mat(kE, mE);
    Eigen::MatrixXd out_tmp(mE, kE); // out_tmp = out_mat.transpose()
    Eigen::VectorXd Bf1st_tmp(kE);
    Eigen::VectorXd bx;
    Eigen::VectorXd Bf1st;

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

    out1 = d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_dindex_par_mat.row(i);
    out2 = dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_df_mat.row(i).transpose() * dlogmu_dindex_par_mat.row(i) + mu(i)*out_mat);

    return - out1 - out2;
  }


  Eigen::MatrixXd he_alpha_f_con_index_par () {
    Eigen::MatrixXd out = he_alpha_f_index_par_mat * dindex_par_dcon_index_par_mat;
    return out;
  }


  Eigen::MatrixXd he_alpha_f_con_index_par_i (int i) {
    Eigen::MatrixXd out = he_alpha_f_index_par_i_mat * dindex_par_dcon_index_par_mat;
    return out;
  }

  Eigen::MatrixXd he_alpha_f_betaF () {
    Eigen::MatrixXd out1(kE, kbetaF);
    Eigen::MatrixXd out2(kE, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * dlogmu_df_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
    }
    return - out1 - out2;
  }

  Eigen::MatrixXd he_alpha_f_betaF_i (int i) {
    Eigen::MatrixXd out1(kE, kbetaF);
    Eigen::MatrixXd out2(kE, kbetaF);
    out1 = d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
    out2 = dlogdensity_dmu_vec(i) * dlogmu_df_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);

    return - out1 - out2;
  }

  Eigen::MatrixXd he_alpha_f_betaR () {
    Eigen::MatrixXd out1(kE, kbetaR);
    Eigen::MatrixXd out2(kE, kbetaR);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_dbetaR_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * dlogmu_df_mat.row(i).transpose() * dmu_dbetaR_mat.row(i);
    }
    return - out1 - out2;
  }


  Eigen::MatrixXd he_alpha_f_betaR_i (int i) {
    Eigen::MatrixXd out1(kE, kbetaR);
    Eigen::MatrixXd out2(kE, kbetaR);

    out1 = d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_dbetaR_mat.row(i);
    out2 = dlogdensity_dmu_vec(i) * dlogmu_df_mat.row(i).transpose() * dmu_dbetaR_mat.row(i);

    return - out1 - out2;
  }


  Eigen::MatrixXd he_con_index_par_betaF () {
    Eigen::MatrixXd out1(mE-1, kbetaF);
    Eigen::MatrixXd out2(mE-1, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += dindex_par_dcon_index_par_mat.transpose() * (d2logdensity_dmudmu_vec(i) * dmu_dindex_par_mat.row(i).transpose() * dmu_dbetaF_mat.row(i));
      out2 += dindex_par_dcon_index_par_mat.transpose() * (dlogdensity_dmu_vec(i) * dlogmu_dindex_par_mat.row(i).transpose() * dmu_dbetaF_mat.row(i));
    }
    return - out1 - out2;
  }


  Eigen::MatrixXd he_con_index_par_betaF_i (int i) {
    Eigen::MatrixXd out1(mE-1, kbetaF);
    Eigen::MatrixXd out2(mE-1, kbetaF);
    out1 = dindex_par_dcon_index_par_mat.transpose() * (d2logdensity_dmudmu_vec(i) * dmu_dindex_par_mat.row(i).transpose() * dmu_dbetaF_mat.row(i));
    out2 = dindex_par_dcon_index_par_mat.transpose() * (dlogdensity_dmu_vec(i) * dlogmu_dindex_par_mat.row(i).transpose() * dmu_dbetaF_mat.row(i));

    return - out1 - out2;
  }



  Eigen::MatrixXd he_betaR_betaF () {
    Eigen::MatrixXd out1(kbetaR, kbetaF);
    Eigen::MatrixXd out2(kbetaR, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * dlogmu_dbetaR_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
    }
    return - out1 - out2;
  }

  Eigen::MatrixXd he_betaR_betaF_i (int i) {
    Eigen::MatrixXd out1(kbetaR, kbetaF);
    Eigen::MatrixXd out2(kbetaR, kbetaF);
    out1 = d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
    out2 = dlogdensity_dmu_vec(i) * dlogmu_dbetaR_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);

    return - out1 - out2;
  }



  Eigen::MatrixXd he_con_index_par_betaR () {
    Eigen::MatrixXd out1(mE-1, kbetaR);
    Eigen::MatrixXd out2(mE-1, kbetaR);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += dindex_par_dcon_index_par_mat.transpose() * (d2logdensity_dmudmu_vec(i) * dmu_dindex_par_mat.row(i).transpose() * dmu_dbetaR_mat.row(i));
      out2 += dindex_par_dcon_index_par_mat.transpose() * (dlogdensity_dmu_vec(i) * dlogmu_dindex_par_mat.row(i).transpose() * dmu_dbetaR_mat.row(i));
    }
    return - out1 - out2;
  }


  Eigen::MatrixXd he_con_index_par_betaR_i (int i) {
    Eigen::MatrixXd out1(mE-1, kbetaR);
    Eigen::MatrixXd out2(mE-1, kbetaR);
    out1 = dindex_par_dcon_index_par_mat.transpose() * (d2logdensity_dmudmu_vec(i) * dmu_dindex_par_mat.row(i).transpose() * dmu_dbetaR_mat.row(i));
    out2 = dindex_par_dcon_index_par_mat.transpose() * (dlogdensity_dmu_vec(i) * dlogmu_dindex_par_mat.row(i).transpose() * dmu_dbetaR_mat.row(i));

    return - out1 - out2;
  }




  Eigen::VectorXd he_alpha_f_log_smoothing_f () {
    return smoothing_f * Sf * alpha_f;
  }

  Eigen::VectorXd he_alpha_f_log_smoothing_w () {
    return smoothing_w * Sw * alpha_f;
  }



  Eigen::MatrixXd he_betaR_logsmoothing () {
    Eigen::MatrixXd out(kbetaR, p);
    out.setZero();
    int begin = 0;
    for (int i = 0; i < p; i++) {
      int ki = static_cast<int>(r(i));
      for (int j = 0; j < ki; j++) out(begin + j, i) = smoothing(i) * betaR(begin + j);
      begin += ki;
    }
    return out;
  }
  Eigen::VectorXd he_alpha_f_log_theta () {
    // he_alpha_f_theta = dmu_df_mat.transpose() * d2logdensity_dmudtheta_vec;
    return -1.0*dmu_df_mat.transpose() * d2logdensity_dmudtheta_vec * theta;
  }
  Eigen::VectorXd he_alpha_f_log_theta_i (int i) {
    // he_alpha_f_theta = dmu_df_mat.transpose() * d2logdensity_dmudtheta_vec;
    return -1.0*dmu_df_mat.row(i).transpose() * d2logdensity_dmudtheta_vec(i) * theta;
  }
  Eigen::VectorXd he_betaR_log_theta () {
    return -1.0*dmu_dbetaR_mat.transpose() * d2logdensity_dmudtheta_vec * theta;
  }
  Eigen::VectorXd he_betaR_log_theta_i (int i) {
    return -1.0*dmu_dbetaR_mat.row(i).transpose() * d2logdensity_dmudtheta_vec(i) * theta;
  }
  Eigen::VectorXd he_betaF_log_theta () {
    return -1.0*dmu_dbetaF_mat.transpose() * d2logdensity_dmudtheta_vec * theta;
  }
  Eigen::VectorXd he_betaF_log_theta_i (int i) {
    return -1.0*dmu_dbetaF_mat.row(i).transpose() * d2logdensity_dmudtheta_vec(i) * theta;
  }

  Eigen::VectorXd he_con_index_par_log_theta () {
    return -1.0*theta * dindex_par_dcon_index_par_mat.transpose() * ( dmu_dindex_par_mat.transpose() * d2logdensity_dmudtheta_vec );
  }
  Eigen::VectorXd he_con_index_par_log_theta_i (int i) {
    return -1.0*theta * dindex_par_dcon_index_par_mat.transpose() * ( dmu_dindex_par_mat.row(i).transpose() * d2logdensity_dmudtheta_vec(i) );
  }





  // *********** LAML ***********
  double logdetH05() {
    // double out = 0.5 * log(he_s_u_mat.determinant());
    // return out;

    double logdetH05 = 0.0;
    Eigen::PartialPivLU<Eigen::MatrixXd> lu(he_s_u_mat);
    Eigen::MatrixXd LU = lu.matrixLU();
    // double c = lu.permutationP().determinant(); // -1 or 1
    double lii;
    for (int i = 0; i < (kE+mE-1+kbetaR+kbetaF); i++) {
      lii = LU(i,i);

      logdetH05 += log(abs(lii));
      // logdetH05 += log(abs(lii) + 1e-8); // add small value to avoid log(0)
      // logdetH05 += log(std::max(lii, 1e-8));
    }
    // logdetH05 += log(c);
    return logdetH05/2.0;

  }

  // copy constructor
  Model(const Model&) = default;
  
};





// optimize alpha_f and betaF for a given index, log_smoothing_f, log_smoothing_w, and log_theta
// Use newton method with eigenvalue modification
void PL(Model& modelobj, bool verbose){
    int maxitr = 50, itr = 0;
    const double eps = 1e-05;
    double mineig = 1e-03; // minimum eigenvalue of Hessian, to ensure it is PD
    int maxstephalve = 50, stephalve = 0;
    int resetitr = 0, maxreset = 1; // if step is nan, reset coefficients as 0. only once.
    int additr = 50; // allow further iterations after resetting coefficients as 0

    Eigen::VectorXd alpha_f = modelobj.alpha_f;
    Eigen::VectorXd betaR = modelobj.betaR;
    Eigen::VectorXd betaF = modelobj.betaF;

    int kE = modelobj.kE;
    int mE = modelobj.mE;
    int kbetaR = modelobj.kbetaR;
    int kbetaF = modelobj.kbetaF;
    int converge = 0;

    int paraSize = kE+kbetaR+kbetaF;
    // Optimize ALPHA_F
    double u;
    double u_tmp;



    // update steps
    Eigen::VectorXd step(paraSize);
    step.setZero();

    Eigen::MatrixXd H(paraSize, paraSize);
    Eigen::VectorXd g(paraSize);

    g.setZero();
    H.setZero();

    // START DEFINE lanczos algorithm for smallest eigenvalue
    // Code from https://github.com/mrcdr/lambda-lanczos/blob/master/src/samples/sample4_use_Eigen_library.cpp
    // the matrix-vector multiplication routine
    auto mv_mul = [&](const vector<double>& in, vector<double>& out) {
      auto eigen_in = Eigen::Map<const Eigen::VectorXd>(&in[0], in.size());
      auto eigen_out = Eigen::Map<Eigen::VectorXd>(&out[0], out.size());

      // eigen_out = H * eigen_in; // Easy version
      eigen_out.noalias() += H * eigen_in; // Efficient version
    };

    LambdaLanczos<double> engine(mv_mul, paraSize, false, 1); // Find 1 minimum eigenvalue
    std::vector<double> smallest_eigenvalues;
    std::vector<std::vector<double>> smallest_eigenvectors;
    double smallest_eigval; // smallest eigenvalue
    // END DEFINE lanczos for smallest eigenvalue

    // eigen decomposition
    Eigen::VectorXd eigvals(paraSize);
    eigvals.setZero();
    Eigen::VectorXd invabseigvals(paraSize);
    invabseigvals.setZero();
    // double eigval;
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(H,false); // Only values, not vectors
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigvec(H,true); // Both values and vectors
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigvec; // Both values and vectors

    double delta = 1.0; // for New Q-Newton method
    double g_norm;

    for (int i = 0; i < paraSize; i++) g(i) = 1. + eps;

    while (itr < maxitr) {
      itr++;
      modelobj.derivative_f();

      g = modelobj.gr_inner_vec;
      g_norm = g.norm();
      if (g_norm < eps) break;
      modelobj.NegativeLogLikelihood();
      u = modelobj.NegLogL;

      H = modelobj.he_inner_mat;
      // std::cout << "H has NaN: " << containsNaN(H) << std::endl;

      engine.run(smallest_eigenvalues, smallest_eigenvectors);
      smallest_eigval = smallest_eigenvalues[0]; // the smallest eigenvalue
      if ((smallest_eigval < 1e-2) || std::isnan(smallest_eigval)) {
        // Do Q-Newton's Step
        eigvec.compute(H, Eigen::ComputeEigenvectors); // Compute eigenvalues and vectors
        eigvals = eigvec.eigenvalues().array();

        if (abs(eigvals.prod()) < 1e-3) {
          for (int iii = 0; iii < paraSize; iii++) eigvals(iii) += delta*g_norm;
        }
        // std::cout << "eigvals" << eigvals.transpose() << std::endl;
        // for (int i = 0; i < paraSize; i++) invabseigvals(i) = 1. / max(abs(eigvals(i)), mineig); // flip signs
        for (int i = 0; i < paraSize; i++) invabseigvals(i) = 1. / abs(eigvals(i)); // flip signs
        step = eigvec.eigenvectors() * (invabseigvals.asDiagonal()) * (eigvec.eigenvectors().transpose()) * g;
      } else {
        // smallest eigenvalue > 1e-3
        // regular Newton's step
        // step = H.llt().solve(g);
        step = H.ldlt().solve(g);
      }

      // check nan in step
      // Really needed
      if(hasNaN(step)){
        if (resetitr < maxreset){
          resetitr++;
          alpha_f.setZero(); // reset alpha_f
          betaR.setZero();
          betaF.setZero();
          // alpha_f.setOnes(); // reset alpha_f
          // betaR.setOnes();
          // betaF.setOnes();
          modelobj.setAlphaF(alpha_f);
          modelobj.setBetaR(betaR);
          modelobj.setBetaF(betaF);
          if(verbose) std::cout << "reset alpha_f and betaF as 0" << std::endl;
          itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
          continue; // do next iteration
        } else {
          converge = 99;
          break;
        }
      }

      alpha_f -= step.segment(0, kE);
      betaR -= step.segment(kE, kbetaR);
      betaF -= step.segment(kE+kbetaR, kbetaF);

      modelobj.setAlphaF(alpha_f);
      modelobj.setBetaR(betaR);
      modelobj.setBetaF(betaF);

      modelobj.NegativeLogLikelihood();

      u_tmp = modelobj.NegLogL;

      // halving if objective function increase.
      stephalve = 0;
      while ((u_tmp > u + 1e-8) & (stephalve < maxstephalve)){
        stephalve++;
        step /= 2.;

        alpha_f += step.segment(0, kE);
        betaR += step.segment(kE, kbetaR);
        betaF += step.segment(kE+kbetaR, kbetaF);
        modelobj.setAlphaF(alpha_f);
        modelobj.setBetaR(betaR);
        modelobj.setBetaF(betaF);

        modelobj.NegativeLogLikelihood();
        u_tmp = modelobj.NegLogL;
      }


      stephalve = 0;
      // Check feasibility of step. If u is nan then we went too far;
      // halve the step and try again
      while (std::isnan(u_tmp) & (stephalve < maxstephalve)) {
        stephalve++;
        step /= 2.; // This is still the step from the previous iteration

        alpha_f += step.segment(0, kE);
        betaR += step.segment(kE, kbetaR);
        betaF += step.segment(kE+kbetaR, kbetaF);
        modelobj.setAlphaF(alpha_f);
        modelobj.setBetaR(betaR);
        modelobj.setBetaF(betaF);
        modelobj.NegativeLogLikelihood();
        u_tmp = modelobj.NegLogL;
      }

      stephalve = 0;
      // if (stephalve > 0) std::cout << "Performed " << stephalve << " iterations of step-halving." << std::endl;
      if (std::isnan(u_tmp)) {
        // Step-halving didn't work
        // std::cout << "AlphaF: Step-halving failed with nan function value. Returning failure." << std::endl;
        converge = 99;
        break;
      }
    }
    if(itr == maxitr){
      // std::cout << "AlphaF: Newton method for updating alpha fails" << std::endl;
      converge = 99;
    }

    // if(verbose) std::cout << "-- AlphaF Gradient Max: " << g.maxCoeff() << std::endl;

    modelobj.derivative_coef();
    modelobj.NegativeLogLikelihood();

    Eigen::VectorXd gr_PL(mE-1);
    gr_PL = modelobj.gr_con_index_par_vec;

    Eigen::MatrixXd he_PL(mE-1, mE-1);




    modelobj.derivative_f();



    Eigen::MatrixXd mat_tmp_PL1(kE+kbetaR+kbetaF,mE-1);

    mat_tmp_PL1.block(0, 0, kE, mE-1) = modelobj.he_alpha_f_con_index_par_mat;
    mat_tmp_PL1.block(kE, 0, kbetaR, mE-1) = modelobj.he_con_index_par_betaR_mat.transpose();
    mat_tmp_PL1.block(kE+kbetaR, 0, kbetaF, mE-1) = modelobj.he_con_index_par_betaF_mat.transpose();


    Eigen::MatrixXd mat1 = modelobj.he_inner_mat.ldlt().solve(mat_tmp_PL1);
    he_PL = modelobj.he_con_index_par_mat - mat1.transpose() * mat_tmp_PL1;



    modelobj.PL_gradient = gr_PL;
    modelobj.PL_hessian = he_PL;
    modelobj.converge = converge; // 0: converge. 99: not converge
}



void Inner(Model& modelobj, bool verbose) {
    // newton method
    // int maxitr = 50, itr = 0;
    int maxitr = 50, itr = 0;
    const double eps = 1e-03;
    const double largereps = 1e-02;
    double mineig = 1e-03; // minimum eigenvalue of Hessian, to ensure it is PD
    // double mineig = 1e-02; // minimum eigenvalue of Hessian, to ensure it is PD
    int maxstephalve = 10, stephalve = 0;
    int increase_maxstephalve = 5; // increase maxstephalve if increasing. 0 means turn-off this feature
    int maxErangehalve = 20;
    int stephalve_inner = 0;
    int resetitr = 0, maxreset = 1; // if step is nan or always diverge, reset coefficients as 0. Only reset once
    int additr = 50; // allow further 50 iterations after resetting coefficients as 0.

    // check non-moving step following https://github.com/awstringer1/varcomptest/blob/main/src/reml-ad.cpp
    int maxnonmovingsteps = 5; // Maximum number of iterations for which we will tolerate no movement. 5 in awstringer1/varcomptest
    double stepeps = 1e-12; // If max(abs(step)) < stepeps then we say the iteration resulted in no movement.
    int stepcounter = 0; // Count the number of non-moving steps

    Eigen::VectorXd con_index_par = modelobj.con_index_par;
    Eigen::MatrixXd Bindex = modelobj.Bindex;
    Eigen::MatrixXd BtBindex = modelobj.BtBindex;

    int mE = modelobj.mE;

    int L = modelobj.L;

    int npar = mE-1;
    int n = modelobj.n;

    // catch double PL.fn
    double s;
    double s_tmp;

    // update steps
    Eigen::VectorXd step(npar);
    step.setZero();


    Eigen::MatrixXd H(npar, npar);
    Eigen::VectorXd g(npar);

    g.setZero();
    H.setZero();


    // START DEFINE lanczos algorithm for smallest eigenvalue
    // Code from https://github.com/mrcdr/lambda-lanczos/blob/master/src/samples/sample4_use_Eigen_library.cpp
    // the matrix-vector multiplication routine
    auto mv_mul = [&](const vector<double>& in, vector<double>& out) {
      auto eigen_in = Eigen::Map<const Eigen::VectorXd>(&in[0], in.size());
      auto eigen_out = Eigen::Map<Eigen::VectorXd>(&out[0], out.size());

      // eigen_out = H * eigen_in; // Easy version
      eigen_out.noalias() += H * eigen_in; // Efficient version
    };

    LambdaLanczos<double> engine(mv_mul, npar, false, 1); // Find 1 minimum eigenvalue
    std::vector<double> smallest_eigenvalues;
    std::vector<std::vector<double>> smallest_eigenvectors;
    double smallest_eigval; // smallest eigenvalue
    // END DEFINE lanczos for smalles eigenvalue


    // eigen decomposition
    Eigen::VectorXd eigvals(npar);
    eigvals.setZero();
    Eigen::VectorXd invabseigvals(npar);
    invabseigvals.setZero();
    // double eigval;
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(H,false); // Only values, not vectors
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigvec(H,true); // Both values and vectors
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigvec; // Both values and vectors
    // check range of E
    Eigen::VectorXd con_index_par_long(mE);
    con_index_par_long(0) = 1.0;
    double index_par_denominator;
    Eigen::VectorXd index_par;


    int kmin;
    int kmax;
    int krange;
    // int krangemin = 4;
    int krangemin = 3;
    bool krangewarning = false;
    std::vector<Eigen::MatrixXd> B_inner_list = modelobj.getB_inner_list();
    Eigen::MatrixXd B_inner = modelobj.getB_inner();
    Eigen::MatrixXd B_inner_tmp(n, L+1);
    B_inner_tmp.setZero();
    Eigen::VectorXd knots_f = modelobj.getknots_f();


    double delta = 1.0; // for New Q-Newton method
    double g_norm;

    // initialize gradient
    for (int i = 0; i < npar; i++) g(i) = 1. + eps;

    if(verbose) std::cout << "* Start optimize profile likelihood" << std::endl;

    // start newton's method
    while (itr < maxitr) {
      itr++;
      // modelobj.NegativeLogLikelihood();

      resetcon_label: // reset index par as all zero. If the current index par always leads to the divergence in updating alpha_f

      PL(modelobj, verbose);
      g = modelobj.PL_gradient;
      g_norm = g.norm();

      if (g_norm < eps) break;

      // std::cout << "** g_norm: " << g_norm << std::endl;

      s = modelobj.NegLogL;
      H = modelobj.PL_hessian;

      // std::cout << "g: " << g.transpose() << std::endl;
      // std::cout << "H: " << H << std::endl;
      // std::cout << "s: " << s << std::endl;
      // std::cout << "step: " << step.transpose() << std::endl;

      engine.run(smallest_eigenvalues, smallest_eigenvectors);
      smallest_eigval = smallest_eigenvalues[0]; // the smallest eigenvalue

      // std::cout << "smallest_eigval: " << smallest_eigval << std::endl;

      if ((smallest_eigval < 1e-2) || std::isnan(smallest_eigval)) {
        // Do Q-Newton's Step
        eigvec.compute(H, Eigen::ComputeEigenvectors); // Compute eigenvalues and vectors
        eigvals = eigvec.eigenvalues().array();

        // std::cout << "eigvals: " << eigvals.transpose() << std::endl;

        if (abs(eigvals.prod()) < 1e-3) {
          for (int iii = 0; iii < npar; iii++) eigvals(iii) += delta*g_norm;
        }

        for (int i = 0; i < npar; i++) invabseigvals(i) = 1. / abs(eigvals(i)); // flip signs
        // std::cout << "invabseigvals max" << invabseigvals.maxCoeff() << std::endl;
        step = eigvec.eigenvectors() * (invabseigvals.asDiagonal()) * (eigvec.eigenvectors().transpose()) * g;
      } else {
        // smallest eigenvalue > 1e-3
        // regular Newton's step
        // step = H.llt().solve(g);
        step = H.ldlt().solve(g);
      }
      // check nan in step
      // NOT really needed. checking here to align with alpha_f.
      if(hasNaN(step)){
        if (resetitr < maxreset){
          resetitr++;
          con_index_par.setOnes();
          modelobj.setCon_Index_Par(con_index_par);
          if(verbose) std::cout << "reset index because of nan step" << std::endl;
          itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
          continue; // do next iteration
        } else {
          break;
        }
      }

      con_index_par -= step;

      // ****** start checking range of E

      for (int j = 0; j < (mE - 1); j++) {
        con_index_par_long(j + 1) = con_index_par(j);
      }
      index_par_denominator = sqrt(con_index_par_long.dot(BtBindex * con_index_par_long));
      index_par = Bindex * con_index_par_long/index_par_denominator;

      B_inner_tmp.setZero();
      for (int j = 0; j < mE; j++) {
        B_inner_tmp += B_inner_list.at(j) * index_par(j);
      }


      kmin = knotindexEigen(B_inner_tmp.minCoeff(), knots_f);
      kmax = knotindexEigen(B_inner_tmp.maxCoeff(), knots_f);
      krange = kmax - kmin;
      stephalve = 0;

      while ((krange < krangemin) & (stephalve < maxErangehalve)){
        stephalve++;
        step /= 2.;

        con_index_par += step;


        for (int j = 0; j < (mE - 1); j++) {
          con_index_par_long(j + 1) = con_index_par(j);
        }
        index_par_denominator = sqrt(con_index_par_long.dot(BtBindex * con_index_par_long));
        index_par = Bindex * con_index_par_long/index_par_denominator;

        B_inner_tmp.setZero();
        for (int j = 0; j < mE; j++) {
          B_inner_tmp += B_inner_list.at(j) * index_par(j);
        }
        kmin = knotindexEigen(B_inner_tmp.minCoeff(), knots_f);
        kmax = knotindexEigen(B_inner_tmp.maxCoeff(), knots_f);
        krange = kmax - kmin;
        // if(verbose) std::cout << "E range krange" << krange << std::endl;
        // std::cout << "halving because of krange" << std::endl;
        // std::cout << "alpha_f" << convertToDouble(modelobj.alpha_f) << std::endl;
      }
      if (stephalve >= maxErangehalve) {
        // std::cout << "index_par" << index_par.transpose() << std::endl;
        // std::cout << "B_inner_tmp.minCoeff()" << B_inner_tmp.minCoeff() << std::endl;
        // std::cout << "B_inner_tmp.maxCoeff()" << B_inner_tmp.maxCoeff() << std::endl;
        // std::cout << "kmin: " << kmin << ", kmax: " << kmax << ", krange: " << krange << std::endl;
        krangewarning = true;
        if(verbose) std::cout << "E range krange: " << krange << " < " << krangemin << std::endl;
        // std::cout << "Range of weighted exposure is small. Consider increasing kE and resetting starting values." << std::endl;
      }
      // finish checking ******
      modelobj.setCon_Index_Par(con_index_par);

      PL(modelobj, verbose);

      // halving if the optimization for alpha_f fails
      stephalve = 0;
      while ((modelobj.converge != 0) & (stephalve < maxstephalve)){
        stephalve++;
        step /= 2.;

        con_index_par += step;
        // if(verbose) std::cout << "halving step because of divergence" << std::endl;
        modelobj.setCon_Index_Par(con_index_par);

        PL(modelobj, verbose);
      }
      if (modelobj.converge != 0) {
        if (resetitr < maxreset){
          resetitr++;

          con_index_par.setOnes();
          modelobj.setCon_Index_Par(con_index_par);

          if(verbose) std::cout << "reset index because of divergence" << std::endl;
          itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
          goto resetcon_label; // do next iteration
        } else {
          std::cout << "Optimization for alpha_f fails" << std::endl;
          break;
        }
      }

      s_tmp = modelobj.NegLogL;
      // std::cout << "s_tmp: " << s_tmp << std::endl;
      // halving if objective function increase.
      stephalve = 0;
      while ((s_tmp > s + 1e-3) & (stephalve < increase_maxstephalve)){
        stephalve++;
        step /= 2.;

        con_index_par += step;
        // if(verbose) std::cout << "halving step because of increase" << std::endl;

        modelobj.setCon_Index_Par(con_index_par);

        PL(modelobj, verbose);
        // when dealing with increase: halving if the optimization for alpha_f fails
        stephalve_inner = 0;
        while ((modelobj.converge != 0) & (stephalve_inner < increase_maxstephalve)){
          stephalve_inner++;
          step /= 2.;

          con_index_par += step;

          modelobj.setCon_Index_Par(con_index_par);

          PL(modelobj, verbose);
        }
        if (modelobj.converge != 0) {
          if (resetitr < maxreset){
            resetitr++;
            con_index_par.setOnes();
            modelobj.setCon_Index_Par(con_index_par);


            if(verbose) std::cout << "reset index because of divergence" << std::endl;
            itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
            goto resetcon_label; // do next iteration
          } else {
            std::cout << "Optimization for alpha_f fails" << std::endl;
            break;
          }
        }
        s_tmp = modelobj.NegLogL;
      }

      stephalve = 0;
      // Check feasibility of step. If u is nan then we went too far;
      // halve the step and try again
      while (std::isnan(s_tmp) & (stephalve < maxstephalve)) {
        stephalve++;
        step /= 2.; // This is still the step from the previous iteration
        // std::cout << "halving step because of nan" << std::endl;
        con_index_par += step;

        modelobj.setCon_Index_Par(con_index_par);

        PL(modelobj, verbose);
        // when dealing with NaN: halving if the optimization for alpha_f fails
        stephalve_inner = 0;
        while ((modelobj.converge != 0) & (stephalve_inner < maxstephalve)){
          stephalve_inner++;
          step /= 2.;


          con_index_par += step;

          modelobj.setCon_Index_Par(con_index_par);
          // std::cout << "halving step because of NaN" << std::endl;
          PL(modelobj, verbose);
        }
        if (modelobj.converge != 0) {
          if (resetitr < maxreset){
            resetitr++;

            con_index_par.setOnes();
            modelobj.setCon_Index_Par(con_index_par);


            if(verbose) std::cout << "reset index because of divergence" << std::endl;
            itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
            goto resetcon_label; // do next iteration
          } else {
            break; // break the current innner while. The next code to run: if (std::isnan(s_tmp)) {...}
          }
        }
        s_tmp = modelobj.NegLogL;
      }

      if (std::isnan(s_tmp)) {
        if (resetitr < maxreset){
          resetitr++;

          con_index_par.setOnes();
          modelobj.setCon_Index_Par(con_index_par);


          if(verbose) std::cout << "reset index because of nan function" << std::endl;
          itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
          goto resetcon_label; // do next iteration
        } else {
          // Step-halving didn't work
          std::cout << "index par: Step-halving failed with nan function value. Returning failure." << std::endl;
          break;
        }
      }
      // Count the number of iterations where we didn't move; if too many, we got stuck.
      // a part of the code follows https://github.com/awstringer1/varcomptest/blob/main/src/reml-ad.cpp
      if (step.lpNorm<Eigen::Infinity>() < stepeps) {
        stepcounter++;
        if (stepcounter > maxnonmovingsteps) {
          if (resetitr < maxreset){
            resetitr++;

            con_index_par.setOnes();
            modelobj.setCon_Index_Par(con_index_par);


            if(verbose) std::cout << "reset index because of non-moving" << std::endl;
            itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
            goto resetcon_label; // do next iteration
          } else {
            std::cout << "The algorithm hasn't moved for " << stepcounter << " steps; terminating. Please check the answer." << std::endl;
            break;
          }
        }
      } else {
        stepcounter = 0; // if move
      }

      // The last reset
      if((itr == maxitr) & (g.norm() >= largereps)) {
        if (resetitr < maxreset){
          resetitr++;

          con_index_par.setOnes();
          modelobj.setCon_Index_Par(con_index_par);


          if(verbose) std::cout << "reset index. The last one" << std::endl;
          itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
          goto resetcon_label; // do next iteration
        } else {
          if (verbose) {
            std::cout << "Newton method for updating weight function might fail. Profile Likelihood Gradient Max: " << g.maxCoeff() << std::endl;
          } else {
            // report less information if verbose is false
            if (g.maxCoeff() >= 1e-2) {
             std::cout << "Newton method for updating weight function might fail. Profile Likelihood Gradient Max: " << g.maxCoeff() << std::endl;
            }
          }
          break;
        }
      }
      // initilized count
      stephalve = 0;
    }
    if(krangewarning) std::cout << "Range of weighted exposure is small. Consider increasing kE or resetting starting values. Profile Likelihood Gradient Max: " << g.maxCoeff() << std::endl;
    if(verbose) std::cout << "* Finish middle opt. Profile Likelihood Gradient Max: " << g.maxCoeff() << std::endl;

    modelobj.PLg = g.maxCoeff();
}



#endif
