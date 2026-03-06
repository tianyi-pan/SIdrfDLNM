/** Include **/
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

#include <random> // for generating samples from standard normal distribution

#ifdef _OPENMP
  #include <omp.h>
#endif


#include "sidDLNMcppadClass.hpp"
#include "sidDLNMeigenClass.hpp"

// #include <LBFGSB.h>
// using namespace LBFGSpp;
// TODO: use https://github.com/yixuan/LBFGSpp


std::vector<Mat> convertListToVectorMat(List rList) {

  int n = rList.size();

  std::vector<Mat> matrixVector;


  for (int i = 0; i < n; ++i) {
    NumericMatrix rMatrix = as<NumericMatrix>(rList[i]);
    Eigen::Map<Eigen::MatrixXd> eigenMatrix(as<Eigen::Map<Eigen::MatrixXd>>(rMatrix));
    matrixVector.push_back(eigenMatrix.cast<Scalar>());
  }

  return matrixVector;
}



std::vector<Eigen::MatrixXd> convertListToVectorEigenMatrixXd(List rList) {

  int n = rList.size();

  std::vector<Eigen::MatrixXd> matrixVector;


  for (int i = 0; i < n; ++i) {
    NumericMatrix rMatrix = as<NumericMatrix>(rList[i]);
    Eigen::Map<Eigen::MatrixXd> eigenMatrix(as<Eigen::Map<Eigen::MatrixXd>>(rMatrix));
    matrixVector.push_back(eigenMatrix);
  }

  return matrixVector;
}



struct LAMLResult {
    double fn;
    Eigen::VectorXd gradient;
};


CppAD::ADFun<double> LAMLTape(Model& modelobj, Modelcppad& modelcppadobj){
    // if(!modelcppadobj.ifhastape) {
      int kE = modelobj.kE;
      int kbetaR = modelobj.kbetaR;
      int kbetaF = modelobj.kbetaF;
      int p = modelobj.p;
      int mE = modelobj.mE;

      Eigen::VectorXd R_alpha_f = modelobj.alpha_f;
      Eigen::VectorXd R_betaR = modelobj.betaR;
      Eigen::VectorXd R_betaF = modelobj.betaF;
      Eigen::VectorXd R_con_index_par = modelobj.con_index_par;

      double R_log_theta = modelobj.log_theta;
      double R_log_smoothing_f = modelobj.log_smoothing_f;
      double R_log_smoothing_w = modelobj.log_smoothing_w;
      Eigen::VectorXd R_logsmoothing = modelobj.logsmoothing;

      Vec alpha_f = R_alpha_f.cast<Scalar>();
      Vec betaR = R_betaR.cast<Scalar>();
      Vec betaF = R_betaF.cast<Scalar>();
      Vec con_index_par = R_con_index_par.cast<Scalar>();

      Scalar log_theta = R_log_theta;
      Scalar log_smoothing_f = R_log_smoothing_f;
      Scalar log_smoothing_w = R_log_smoothing_w;
      Vec logsmoothing = R_logsmoothing.cast<Scalar>();

      // alpha_f, con_index_par, betaR, betaF, log_theta, log_smoothing_f, log_smoothing_w, logsmoothing


      Vec at(kE+kbetaR+kbetaF + mE - 1 + 3 + p);
      at << alpha_f, con_index_par, betaR, betaF, log_theta, log_smoothing_f, log_smoothing_w, logsmoothing;


      // start reverse mode
      Vec result(1);
      CppAD::ADFun<double> gr;

      CppAD::Independent(at);
      modelcppadobj.setAlphaF(at.segment(0, kE));
      modelcppadobj.setCon_Index_Par(at.segment(kE, mE-1));
      modelcppadobj.setBetaR(at.segment(kE+mE-1, kbetaR));
      modelcppadobj.setBetaF(at.segment(kE+mE-1+kbetaR, kbetaF));
      modelcppadobj.setLogTheta(at(kE+kbetaR+kbetaF+mE-1));
      modelcppadobj.setLogSmoothingF(at(kE+kbetaR+kbetaF+mE-1+1));
      modelcppadobj.setLogSmoothingW(at(kE+kbetaR+kbetaF+mE-1+2));
      modelcppadobj.setLogsmoothing(at.segment(kE+kbetaR+kbetaF+mE-1+3, p));


      modelcppadobj.derivative_coef();
      modelcppadobj.derivative_he();
      result(0) = modelcppadobj.logdetH05();

      gr.Dependent(at, result);
      // gr.optimize();
      // modelcppadobj.gr = gr;
      return gr;
      // END reverse



      // Forward mode 1
      // std::vector<Scalar> at_std(at.data(), at.data() + (kE+kbetaR+kbetaF + 3 + p));

      // CppAD::Independent(at_std);

      // std::vector<Scalar> result(1);
      // Vec at_eigen = Eigen::Map<Vec>(at_std.data(), (kE+kbetaR+kbetaF + 3 + p));
      // modelcppadobj.setAlphaF(at_eigen.segment(0, kE));
      // modelcppadobj.setCon_Index_Par(at_eigen.segment(kE, mE-1));
      // modelcppadobj.setBetaR(at_eigen.segment(kE+mE-1, kbetaR));
      // modelcppadobj.setBetaF(at_eigen.segment(kE+mE-1+kbetaR, kbetaF));
      // modelcppadobj.setLogTheta(at_eigen(kE+kbetaR+kbetaF+mE-1));
      // modelcppadobj.setLogSmoothingF(at_eigen(kE+kbetaR+kbetaF+mE-1+1));
      // modelcppadobj.setLogSmoothingW(at_eigen(kE+kbetaR+kbetaF+mE-1+2));
      // modelcppadobj.setLogsmoothing(at_eigen.segment(kE+kbetaR+kbetaF+mE-1+3, p));

      // modelcppadobj.derivative_coef();
      // modelcppadobj.derivative_he();
      // result[0] = modelcppadobj.logdetH05();

      // CppAD::ADFun<double> gr(at_std, result);
      // return gr;
      // END forward 1


      // Forward mode 2
      // CppAD::Independent(at);

      // Vec result(1);
      // modelcppadobj.setAlphaF(at.segment(0, kE));
      // modelcppadobj.setCon_Index_Par(at.segment(kE, mE-1));
      // modelcppadobj.setBetaR(at.segment(kE+mE-1, kbetaR));
      // modelcppadobj.setBetaF(at.segment(kE+mE-1+kbetaR, kbetaF));
      // modelcppadobj.setLogTheta(at(kE+kbetaR+kbetaF+mE-1));
      // modelcppadobj.setLogSmoothingF(at(kE+kbetaR+kbetaF+mE-1+1));
      // modelcppadobj.setLogSmoothingW(at(kE+kbetaR+kbetaF+mE-1+2));
      // modelcppadobj.setLogsmoothing(at.segment(kE+kbetaR+kbetaF+mE-1+3, p));


      // modelcppadobj.derivative_coef();
      // modelcppadobj.derivative_he();
      // result(0) = modelcppadobj.logdetH05();

      // CppAD::ADFun<double> gr(at, result);
      // return gr;
      // END forward 2


      // modelcppadobj.ifhastape = true;





    // }
}

LAMLResult LAML(Model& modelobj, Modelcppad& modelcppadobj, bool ifgradient) {
    LAMLResult result;
    double u_LAML;
    int kE = modelobj.kE;
    int kbetaR = modelobj.kbetaR;
    int kbetaF = modelobj.kbetaF;
    int mE = modelobj.mE;
    int p = modelobj.p;
    modelobj.derivative_coef();
    modelobj.derivative_he();
    modelobj.derivative_full();
    modelobj.NegativeLogLikelihood();

    Eigen::VectorXd alpha_f = modelobj.alpha_f;
    Eigen::VectorXd betaR = modelobj.betaR;
    Eigen::VectorXd betaF = modelobj.betaF;
    Eigen::VectorXd con_index_par = modelobj.con_index_par;

    Eigen::VectorXd gr_s_u_vec = modelobj.gr_s_u_vec;
    Eigen::VectorXd gr_s_par_vec = modelobj.gr_s_par_vec;
    Eigen::MatrixXd he_s_u_mat = modelobj.he_s_u_mat;
    Eigen::MatrixXd he_s_par_u_mat = modelobj.he_s_par_u_mat;


    double log_theta = modelobj.log_theta;
    double log_smoothing_f = modelobj.log_smoothing_f;
    double log_smoothing_w = modelobj.log_smoothing_w;
    Eigen::VectorXd logsmoothing = modelobj.logsmoothing;
    Eigen::VectorXd gr(3 + p);

    if(ifgradient) {
      // First derivative of LAML

      CppAD::ADFun<double> cppadgr = LAMLTape(modelobj, modelcppadobj);



      Eigen::VectorXd at0(kE+kbetaR+kbetaF + mE-1 + 3 + p);
      at0 << alpha_f, con_index_par, betaR, betaF, log_theta, log_smoothing_f, log_smoothing_w, logsmoothing;

      // // reverse mode Jacobian
      // Eigen::VectorXd g_LAML(kE + mE-1 + 3 + kbetaR+kbetaF + p);
      // g_LAML.setZero();
      // g_LAML = cppadgr.Jacobian(at0);
      // // END reverse mode Jacobian

      // using reverse mode Reverse(1).
      // Should be faster than Jacobian when dim(output) is 1.
      std::vector<double> x(at0.data(), at0.data() + at0.size());
      // 1) Evaluate f at x (zero-order forward)
      cppadgr.Forward(0, x);
      // 2) Seed the single output with 1.0 (since output is scalar)
      std::vector<double> w(1, 1.0);
      // 3) One reverse sweep gives the whole gradient
      std::vector<double> grad = cppadgr.Reverse(1, w);
      // free cppadgr
      cppadgr = CppAD::ADFun<double>();
      // 4) Map to Eigen
      Eigen::Map<const Eigen::VectorXd> g_LAML_map(grad.data(), static_cast<Eigen::Index>(grad.size()));
      Eigen::VectorXd g_LAML = g_LAML_map;
      // END using reverse mode Reverse(1)


      // forward mode 1, 2
      // std::vector<Scalar> at0_std(at0.data(), at0.data() + (kE + mE-1 + kbetaR + kbetaF + 3 + p));
      // std::vector<double> dx((kE + mE-1 + kbetaR + kbetaF + 3 + p), 0.0);
      // std::vector<double> g_LAML_std((kE + mE-1 + kbetaR + kbetaF + 3 + p));
      // std::vector<double> dy(1);
      // for (size_t i = 0; i < (kE + mE-1 + kbetaR + kbetaF + 3 + p); ++i) {
      //     dx[i] = 1.0;
      //     dy = cppadgr.Forward(1, dx);
      //     g_LAML_std[i] = dy[0];
      //     dx[i] = 0.0;
      // }
      // Eigen::VectorXd g_LAML = Eigen::Map<Eigen::VectorXd>(g_LAML_std.data(), g_LAML_std.size());
      // END forward mode 1, 2



      // std::cout << "modelobj.logdetH05()" << modelobj.logdetH05() << std::endl;
      // std::cout << "modelobj.NegLogL" << modelobj.NegLogL << std::endl;
      // std::cout << "u_LAML: " << u_LAML << std::endl;

      // In R: grad[-(1:(kE+kw-1))] - H.full[-(1:(kE+kw-1)),(1:(kE+kw-1))] %*% as.vector(solve(H.alpha, grad[(1:(kE+kw-1))]))
      Eigen::VectorXd g1 = g_LAML.segment(0, kE + mE-1+kbetaR+kbetaF) + gr_s_u_vec;
      Eigen::VectorXd g2 = g_LAML.segment(kE+mE-1+kbetaR+kbetaF, 3+p) + gr_s_par_vec;
      gr = g2 - he_s_par_u_mat * he_s_u_mat.ldlt().solve(g1);

    } else {
      gr = Eigen::VectorXd::Zero(3 + p);
    }

    u_LAML = modelobj.logdetH05() + modelobj.NegLogL - modelobj.n/2.0 * log(2*3.141592653589793238462643383279);

    result.fn = u_LAML;
    result.gradient = gr;
    return result;
}




// [[Rcpp::export]]
List sidDLNMbuild(const Eigen::VectorXd R_y,
                   const List R_B_inner_list,
                   const Eigen::VectorXd R_knots_f,
                   const Eigen::VectorXd R_knots_w,
                   const Eigen::MatrixXd R_Sw,
                   const Eigen::MatrixXd R_SwR,
                   const Eigen::MatrixXd R_Sf,
                   const Eigen::MatrixXd R_SfR,
                   const Eigen::MatrixXd R_Xrand,
                   const Eigen::MatrixXd R_Xfix,
                   const Eigen::MatrixXd R_Zf,
                   const Eigen::VectorXd R_Xoffset,
                   const Eigen::VectorXd R_r,
                   const Eigen::MatrixXd R_Bindex,
                   Eigen::VectorXd R_alpha_f,
                   Eigen::VectorXd R_con_index_par,
                   double R_log_theta,
                   double R_log_smoothing_f,
                   double R_log_smoothing_w,
                   Eigen::VectorXd R_betaR,
                   Eigen::VectorXd R_betaF,
                   Eigen::VectorXd R_logsmoothing) {
    // convert
    Vec y = R_y.cast<Scalar>();
    std::vector<Mat> B_inner_list = convertListToVectorMat(R_B_inner_list);
    std::vector<Eigen::MatrixXd> R_B_inner_Eigen_list = convertListToVectorEigenMatrixXd(R_B_inner_list);
    Vec knots_f = R_knots_f.cast<Scalar>();
    Vec knots_w = R_knots_w.cast<Scalar>();
    Mat Sw = R_Sw.cast<Scalar>();
    Mat SwR = R_SwR.cast<Scalar>();
    Mat Sf = R_Sf.cast<Scalar>();
    Mat SfR = R_SfR.cast<Scalar>();
    Mat Xrand = R_Xrand.cast<Scalar>();
    Mat Xfix = R_Xfix.cast<Scalar>();
    Mat Zf = R_Zf.cast<Scalar>();
    Vec Xoffset = R_Xoffset.cast<Scalar>();
    Vec r = R_r.cast<Scalar>();
    Mat Bindex = R_Bindex.cast<Scalar>();
    Vec alpha_f = R_alpha_f.cast<Scalar>();
    Vec con_index_par = R_con_index_par.cast<Scalar>();
    Scalar log_theta = R_log_theta;
    Scalar log_smoothing_f = R_log_smoothing_f;
    Scalar log_smoothing_w = R_log_smoothing_w;
    Vec betaR = R_betaR.cast<Scalar>();
    Vec betaF = R_betaF.cast<Scalar>();
    Vec logsmoothing = R_logsmoothing.cast<Scalar>();


    Modelcppad* modelcppadobj_ptr = new Modelcppad(y, B_inner_list, knots_f, knots_w, Sw, SwR, Sf, SfR,
                                    Xrand, Xfix, Zf, Xoffset, r, Bindex,
                                    alpha_f, con_index_par, log_theta, log_smoothing_f, log_smoothing_w,
                                    betaR, betaF, logsmoothing);
    Rcpp::XPtr<Modelcppad> ptrcppad(modelcppadobj_ptr);

    Model* modelobj_ptr = new Model(R_y, R_B_inner_Eigen_list, R_knots_f, R_knots_w, R_Sw, R_SwR, R_Sf, R_SfR,
                                R_Xrand, R_Xfix, R_Zf, R_Xoffset, R_r, R_Bindex,
                                R_alpha_f, R_con_index_par, R_log_theta, R_log_smoothing_f, R_log_smoothing_w,
                                R_betaR, R_betaF, R_logsmoothing);
    Rcpp::XPtr<Model> ptr(modelobj_ptr);


    return List::create(Named("address.eigen") = ptr,
                        Named("address.cppad") = ptrcppad);
}






// [[Rcpp::export]]
List sidDLNMopt(SEXP ptr,
                SEXP ptrcppad,
                Eigen::VectorXd R_alpha_f,
                Eigen::VectorXd R_con_index_par,
                double R_log_theta,
                double R_log_smoothing_f,
                double R_log_smoothing_w,
                Eigen::VectorXd R_betaR,
                Eigen::VectorXd R_betaF,
                Eigen::VectorXd R_logsmoothing,
                bool ifgradient,
                bool verbose) {


    Rcpp::XPtr<Model> modelobj_ptr(ptr);
    Rcpp::XPtr<Modelcppad> modelcppadobj_ptr(ptrcppad);

    Model& modelobj = *modelobj_ptr;
    Modelcppad& modelcppadobj = *modelcppadobj_ptr;

    modelobj.setAlphaF(R_alpha_f);
    modelobj.setCon_Index_Par(R_con_index_par);
    modelobj.setBetaR(R_betaR);
    modelobj.setBetaF(R_betaF);
    modelobj.setLogTheta(R_log_theta);
    modelobj.setLogSmoothingF(R_log_smoothing_f);
    modelobj.setLogSmoothingW(R_log_smoothing_w);
    modelobj.setLogsmoothing(R_logsmoothing);

    // test PL
    // PL(modelobj, verbose);

    LAMLResult LAMLresult;
    // // Inner opt.
    Inner(modelobj, verbose);
    // // get gr of LAML
    LAMLresult = LAML(modelobj, modelcppadobj, ifgradient);

    return List::create(Named("LAML.fn") = LAMLresult.fn,
                        Named("LAML.gradient") = LAMLresult.gradient,
                        Named("alpha_f.mod") = modelobj.alpha_f,
                        Named("con_index_par.mod") = modelobj.con_index_par,
                        Named("betaR.mod") = modelobj.betaR,
                        Named("betaF.mod") = modelobj.betaF,
                        Named("PLg") = modelobj.PLg,
                        Named("address") = modelobj_ptr
                        );
}




// [[Rcpp::export]]
List sidDLNMCI(SEXP ptr,
                const int Rci,
                const int rseed,
                bool ifeta,
                bool verbose,
                Rcpp::Nullable<Rcpp::NumericMatrix> R_he_input = R_NilValue) {

  Rcpp::XPtr<Model> modelobj_ptr(ptr);
  Model& modelobj = *modelobj_ptr;


  Eigen::VectorXd R_alpha_f = modelobj.alpha_f;
  Eigen::VectorXd R_con_index_par = modelobj.con_index_par;
  Eigen::VectorXd R_betaR = modelobj.betaR;
  Eigen::VectorXd R_betaF = modelobj.betaF;

  Eigen::MatrixXd R_Bindex = modelobj.Bindex;
  Eigen::MatrixXd R_BtBindex = modelobj.BtBindex;

  int n = modelobj.n;
  int kE = modelobj.kE;
  int mE = modelobj.mE;
  int kbetaR = modelobj.kbetaR;
  int kbetaF = modelobj.kbetaF;
  Eigen::MatrixXd R_B_inner(n, kE);

  int paraSize = kE+mE-1+kbetaR+kbetaF;
  int paraSizefull;

  // hessian
  Eigen::MatrixXd R_he;
  Eigen::VectorXd R_index_par(mE);

  // Vectors for sampling
  Eigen::VectorXd R_alpha_f_sample(kE);
  Eigen::VectorXd R_betaR_sample(kbetaR);
  Eigen::VectorXd R_betaF_sample(kbetaF);

  Eigen::VectorXd R_con_index_par_sample(mE-1);
  Eigen::VectorXd R_index_par_sample(mE);

  // Matrices to save results
  Eigen::MatrixXd con_index_par_sample_mat(Rci, mE-1);
  Eigen::MatrixXd index_par_sample_mat(Rci, mE);
  Eigen::MatrixXd alpha_f_sample_mat(Rci, kE);
  Eigen::MatrixXd betaR_sample_mat(Rci, kbetaR);
  Eigen::MatrixXd betaF_sample_mat(Rci, kbetaF);

  // components for eta
  Eigen::VectorXd R_E;
  Eigen::VectorXd eta_sample;
  Eigen::MatrixXd eta_sample_mat;
  Eigen::VectorXd eta_E_sample;
  Eigen::VectorXd eta_other_sample;
  Eigen::MatrixXd eta_E_sample_mat;
  Eigen::MatrixXd eta_other_sample_mat;

  std::vector<Eigen::MatrixXd> R_B_inner_list = modelobj.getB_inner_list();
  Eigen::VectorXd R_knots_f = modelobj.getknots_f();
  Eigen::VectorXd R_knots_w = modelobj.getknots_w();
  Eigen::MatrixXd R_Zf = modelobj.getZf();
  Eigen::MatrixXd R_Xfix = modelobj.getXfix();
  Eigen::MatrixXd R_Xrand = modelobj.getXrand();
  Eigen::VectorXd R_Xoffset = modelobj.getXoffset();

  if(ifeta) {
    R_E.resize(n);
    eta_sample.resize(n);
    eta_sample_mat.resize(Rci, n);
    eta_E_sample.resize(n);
    eta_other_sample.resize(n);
    eta_E_sample_mat.resize(Rci, n);
    eta_other_sample_mat.resize(Rci, n);
  }



  Eigen::VectorXd R_con_index_par_long(mE); // con_index_par_long = c(1,con_index_par)
  R_con_index_par_long(0) = 1.0;



  paraSizefull = paraSize;

  int L = modelobj.L;
  int kl = modelobj.kl;
  int kx = modelobj.kx;
  Eigen::VectorXd Lseq(L+1); // Lseq = 0:L
  Eigen::MatrixXd Blag(kl, L+1); // basis for lag
  for (int l = 0; l < (L+1); l++) {
    Lseq(l) = (double) l;
  }

  for (int j = 0; j < (L+1); j++) {
    Blag.col(j) = Bsplinevec(Lseq(j), R_knots_w, 4);
  }
  Eigen::VectorXd Bf(kE);
  Eigen::VectorXd Bf_tmp(kE);

  // Joint
  if (R_he_input.isNotNull()) {
    R_he = Rcpp::as<Eigen::MatrixXd>(R_he_input.get());
  } else {
    R_he = modelobj.he_s_u_mat;
  }

  Eigen::VectorXd R_u_mod(paraSizefull);
  // Hessian
  // cholesky of inverse Hessian
  Eigen::MatrixXd R_he_u_L(paraSize, paraSize);
  Eigen::MatrixXd R_he_u_L_inv(paraSizefull, paraSize);
  Eigen::VectorXd zjoint(paraSize);
  Eigen::VectorXd samplejoint(paraSizefull);


  R_u_mod << R_alpha_f, R_con_index_par, R_betaR, R_betaF;

  // cholesky of inverse Hessian
  R_he_u_L = R_he.llt().matrixL();

  R_he_u_L_inv = (invertL(R_he_u_L)).transpose();


  // std::random_device rd;
  // std::mt19937 gen(rd());
  std::mt19937 gen(rseed);
  std::normal_distribution<> dist(0, 1);


  for(int i = 0; i < Rci; i++)
  {
    // Jointly sample
    for (int j = 0; j < paraSize; j++) {
      zjoint(j) = dist(gen);
    }
    samplejoint = R_u_mod + R_he_u_L_inv * zjoint;
    // get alpha_f
    R_alpha_f_sample = samplejoint.segment(0, kE);


    // get index par
    R_con_index_par_sample = samplejoint.segment(kE, mE-1);
    for (int j = 0; j < (mE - 1); j++) {
      R_con_index_par_long(j + 1) = R_con_index_par_sample(j);
    }
    R_index_par_sample = R_Bindex * R_con_index_par_long/sqrt(R_con_index_par_long.dot(R_BtBindex * R_con_index_par_long));

    // get betaR
    R_betaR_sample = samplejoint.segment(kE+mE-1, kbetaR);
    // get betaF
    R_betaF_sample = samplejoint.segment(kE+mE-1+kbetaR, kbetaF);

    // save
    con_index_par_sample_mat.row(i) = R_con_index_par_sample.transpose();

    if(ifeta) {
      R_B_inner.setZero();
      for (int ii = 0; ii < mE; ii++) {
        R_B_inner += R_B_inner_list.at(ii) * R_index_par_sample(ii);
      }

      for (int ii = 0; ii < n; ii++) {
        Bf.setZero();
        Bf_tmp.setZero();
        for (int j = 0; j < (L+1); j++) {
          Eigen::VectorXd bx = BsplinevecCon(R_B_inner(ii, j), R_knots_f, 4, R_Zf);
          for (int iii = 0; iii < kx; iii++) {
            // TODO: optimize the computation of Bf_tmp
            Bf_tmp.segment(iii*kl, kl) = bx(iii) * Blag.col(j);
          }
          Bf += Bf_tmp;
        }

        eta_E_sample(ii) = Bf.dot(R_alpha_f_sample);
        eta_other_sample(ii) = R_Xfix.row(ii).dot(R_betaF_sample) + R_Xrand.row(ii).dot(R_betaR_sample); // + R_Xoffset(ii);
        eta_sample(ii) = eta_E_sample(ii) + eta_other_sample(ii) + R_Xoffset(ii);
      }

      eta_E_sample = eta_E_sample + eta_other_sample.mean()*Eigen::VectorXd::Ones(n);
      eta_other_sample = eta_other_sample - eta_other_sample.mean()*Eigen::VectorXd::Ones(n);

      eta_E_sample_mat.row(i) = eta_E_sample.transpose();
      eta_other_sample_mat.row(i) = eta_other_sample.transpose();
      eta_sample_mat.row(i) = eta_sample.transpose();
    }

    // save
    alpha_f_sample_mat.row(i) = R_alpha_f_sample.transpose();
    betaR_sample_mat.row(i) = R_betaR_sample.transpose();
    betaF_sample_mat.row(i) = R_betaF_sample.transpose();
    index_par_sample_mat.row(i) = R_index_par_sample.transpose();
  }

  // point estimate for eta
  Eigen::VectorXd eta_point(n);
  Eigen::VectorXd eta_E_point(n);
  Eigen::VectorXd eta_other_point(n);

  if(ifeta) {
    for (int j = 0; j < (mE - 1); j++) {
      R_con_index_par_long(j + 1) = R_con_index_par(j);
    }
    R_index_par = R_Bindex * R_con_index_par_long/sqrt(R_con_index_par_long.dot(R_BtBindex * R_con_index_par_long));
    R_B_inner.setZero();
    for (int ii = 0; ii < mE; ii++) {
      R_B_inner += R_B_inner_list.at(ii) * R_index_par(ii);
    }

    for (int ii = 0; ii < n; ii++) {
      Bf.setZero();
      Bf_tmp.setZero();
      for (int j = 0; j < (L+1); j++) {
        Eigen::VectorXd bx = BsplinevecCon(R_B_inner(ii, j), R_knots_f, 4, R_Zf);
        for (int iii = 0; iii < kx; iii++) {
          // TODO: optimize the computation of Bf_tmp
          Bf_tmp.segment(iii*kl, kl) = bx(iii) * Blag.col(j);
        }
        Bf += Bf_tmp;
      }

      eta_E_point(ii) = Bf.dot(R_alpha_f);
      eta_other_point(ii) = R_Xfix.row(ii).dot(R_betaF) + R_Xrand.row(ii).dot(R_betaR);
      eta_point(ii) = eta_E_point(ii) + eta_other_point(ii) + R_Xoffset(ii);
    }

    eta_E_point = eta_E_point + eta_other_point.mean()*Eigen::VectorXd::Ones(n);
    eta_other_point = eta_other_point - eta_other_point.mean()*Eigen::VectorXd::Ones(n);
  }


  return List::create(Named("index_sample") = index_par_sample_mat,
                    Named("alpha_f_sample") = alpha_f_sample_mat,
                    Named("betaR_sample") = betaR_sample_mat,
                    Named("betaF_sample") = betaF_sample_mat,
                    Named("eta_sample_mat") = eta_sample_mat,
                    Named("eta_E_sample_mat") = eta_E_sample_mat,
                    Named("eta_other_sample_mat") = eta_other_sample_mat,
                    Named("eta_point") = eta_point,
                    Named("eta_E_point") = eta_E_point,
                    Named("eta_other_point") = eta_other_point,
                    Named("Hessian_inner") = R_he);


}

CppAD::ADFun<double> LAMLTapeAIC(Model& modelobj, Modelcppad& modelcppadobj){
      int kE = modelobj.kE;
      int kbetaR = modelobj.kbetaR;
      int kbetaF = modelobj.kbetaF;
      int p = modelobj.p;
      int mE = modelobj.mE;

      Eigen::VectorXd R_alpha_f = modelobj.alpha_f;
      Eigen::VectorXd R_betaR = modelobj.betaR;
      Eigen::VectorXd R_betaF = modelobj.betaF;
      Eigen::VectorXd R_con_index_par = modelobj.con_index_par;

      double R_log_theta = modelobj.log_theta;
      double R_log_smoothing_f = modelobj.log_smoothing_f;
      double R_log_smoothing_w = modelobj.log_smoothing_w;
      Eigen::VectorXd R_logsmoothing = modelobj.logsmoothing;

      Vec alpha_f = R_alpha_f.cast<Scalar>();
      Vec betaR = R_betaR.cast<Scalar>();
      Vec betaF = R_betaF.cast<Scalar>();
      Vec con_index_par = R_con_index_par.cast<Scalar>();

      Scalar log_theta = R_log_theta;
      Scalar log_smoothing_f = R_log_smoothing_f;
      Scalar log_smoothing_w = R_log_smoothing_w;
      Vec logsmoothing = R_logsmoothing.cast<Scalar>();

      // // dynamic parameters:
      // int n_theta = kE+mE-1+kbetaR+kbetaF + 1;
      // Eigen::VectorXd v =  Eigen::VectorXd::Zero(2*n_theta); // all zero with length kE+mE-1+kbetaR+kbetaF + 1
      
      
      
      // Vec v_ad = v.cast<Scalar>();
      // Vec vi_ad(n_theta);
      // Vec vj_ad(n_theta);
      
      // Mat IS_mat(kE+mE-1+kbetaR+kbetaF + 1, kE+mE-1+kbetaR+kbetaF + 1);
      // Mat R(kE+mE-1+kbetaR+kbetaF + 1, kE+mE-1+kbetaR+kbetaF + 1);
      // Mat I_mat(kE+mE-1+kbetaR+kbetaF + 1, kE+mE-1+kbetaR+kbetaF + 1);
      // IS_mat.setZero();
      // R.setZero();
      // I_mat.setZero();

      Vec at(kE+kbetaR+kbetaF + mE);
      at << alpha_f, con_index_par, betaR, betaF, log_theta; //, log_smoothing_f, log_smoothing_w, logsmoothing;


      // Forward mode 2
      // CppAD::Independent(at, 0, true, v_ad);
      // Vec result(1);
      CppAD::Independent(at, 0, true);
      Vec result((kE+mE-1+kbetaR+kbetaF + 1) * (kE+mE-1+kbetaR+kbetaF + 2) / 2);

      modelcppadobj.setAlphaF(at.segment(0, kE));
      modelcppadobj.setCon_Index_Par(at.segment(kE, mE-1));
      modelcppadobj.setBetaR(at.segment(kE+mE-1, kbetaR));
      modelcppadobj.setBetaF(at.segment(kE+mE-1+kbetaR, kbetaF));
      modelcppadobj.setLogTheta(at(kE+kbetaR+kbetaF+mE-1));
      // modelcppadobj.setLogSmoothingF(at(kE+kbetaR+kbetaF+mE-1+1));
      // modelcppadobj.setLogSmoothingW(at(kE+kbetaR+kbetaF+mE-1+2));
      // modelcppadobj.setLogsmoothing(at.segment(kE+kbetaR+kbetaF+mE-1+3, p));

      modelcppadobj.derivative_coef();
      // modelcppadobj.derivative_he();
      modelcppadobj.prepare_AIC();


      // IS_mat.block(0,0,kE+mE-1+kbetaR+kbetaF, kE+mE-1+kbetaR+kbetaF) = modelcppadobj.he_s_u_mat;
      // IS_mat.row(kE+mE-1+kbetaR+kbetaF) << modelcppadobj.he_alpha_f_log_theta_vec.transpose(), modelcppadobj.he_con_index_par_log_theta_vec.transpose(), modelcppadobj.he_betaR_log_theta_vec.transpose(), modelcppadobj.he_betaF_log_theta_vec.transpose(), modelcppadobj.he_log_theta_scalar;
      // IS_mat.col(kE+mE-1+kbetaR+kbetaF) = IS_mat.row(kE+mE-1+kbetaR+kbetaF).transpose();

      // Eigen::LLT<Mat> llt(IS_mat.ldlt().solve(Mat::Identity(kE+mE-1+kbetaR+kbetaF + 1, kE+mE-1+kbetaR+kbetaF + 1)).selfadjointView<Eigen::Upper>());
      // R = llt.matrixU();  // Upper-triangular R so that inv(IS_mat) = R^T * R; nrow(R) = kE+mE-1+kbetaR+kbetaF + 1
      // I_mat = modelcppadobj.I_mat;
      int k = 0;
      for (int j = 0; j < (kE+mE-1+kbetaR+kbetaF + 1); ++j) {
        for (int i = 0; i <= j; ++i) {
          // result[k++] = R(i, j);
          result[k++] = modelcppadobj.I_mat(i, j);
        }
      }
      // for (int i = 0; i < n_theta; ++i) {
      //   vi_ad[i] = v_ad[i];               // first half = vi
      //   vj_ad[i] = v_ad[i + n_theta];     // second half = vj
      // }
      // result[0] = vi_ad.dot(modelcppadobj.I_mat * vj_ad);
      CppAD::ADFun<double> gr(at, result);

      return gr;
      // END forward 2

}





CppAD::ADFun<double> LAMLTapeAICproposed(Model& modelobj, Modelcppad& modelcppadobj){
      int kE = modelobj.kE;
      int kbetaR = modelobj.kbetaR;
      int kbetaF = modelobj.kbetaF;
      int p = modelobj.p;
      int mE = modelobj.mE;

      Eigen::VectorXd R_alpha_f = modelobj.alpha_f;
      Eigen::VectorXd R_betaR = modelobj.betaR;
      Eigen::VectorXd R_betaF = modelobj.betaF;
      Eigen::VectorXd R_con_index_par = modelobj.con_index_par;

      double R_log_theta = modelobj.log_theta;
      double R_log_smoothing_f = modelobj.log_smoothing_f;
      double R_log_smoothing_w = modelobj.log_smoothing_w;
      Eigen::VectorXd R_logsmoothing = modelobj.logsmoothing;

      Vec alpha_f = R_alpha_f.cast<Scalar>();
      Vec betaR = R_betaR.cast<Scalar>();
      Vec betaF = R_betaF.cast<Scalar>();
      Vec con_index_par = R_con_index_par.cast<Scalar>();

      Scalar log_theta = R_log_theta;
      Scalar log_smoothing_f = R_log_smoothing_f;
      Scalar log_smoothing_w = R_log_smoothing_w;
      Vec logsmoothing = R_logsmoothing.cast<Scalar>();




      // Mat IS_mat(kE+mE-1+kbetaR+kbetaF + 1, kE+mE-1+kbetaR+kbetaF + 1);
      // Mat I_mat(kE+mE-1+kbetaR+kbetaF + 1, kE+mE-1+kbetaR+kbetaF + 1);
      Mat Khat(kE+mE-1+kbetaR+kbetaF + 1, kE+mE-1+kbetaR+kbetaF + 1);
      // Mat R(kE+mE-1+kbetaR+kbetaF + 1, kE+mE-1+kbetaR+kbetaF + 1);
      // IS_mat.setZero();
      // I_mat.setZero();
      Khat.setZero();
      // R.setZero();

      Mat IdenMat = Mat::Identity(kE+mE-1+kbetaR+kbetaF + 1, kE+mE-1+kbetaR+kbetaF + 1);

      Vec at(kE+kbetaR+kbetaF + mE);
      at << alpha_f, con_index_par, betaR, betaF, log_theta; //, log_smoothing_f, log_smoothing_w, logsmoothing;


      // Forward mode 2
      CppAD::Independent(at);
      // Vec result((kE+mE-1+kbetaR+kbetaF + 1) * (kE+mE-1+kbetaR+kbetaF + 2) / 2);
      Vec result((kE+mE-1+kbetaR+kbetaF + 1) * (kE+mE-1+kbetaR+kbetaF + 1));

      modelcppadobj.setAlphaF(at.segment(0, kE));
      modelcppadobj.setCon_Index_Par(at.segment(kE, mE-1));
      modelcppadobj.setBetaR(at.segment(kE+mE-1, kbetaR));
      modelcppadobj.setBetaF(at.segment(kE+mE-1+kbetaR, kbetaF));
      modelcppadobj.setLogTheta(at(kE+kbetaR+kbetaF+mE-1));
      // modelcppadobj.setLogSmoothingF(at(kE+kbetaR+kbetaF+mE-1+1));
      // modelcppadobj.setLogSmoothingW(at(kE+kbetaR+kbetaF+mE-1+2));
      // modelcppadobj.setLogsmoothing(at.segment(kE+kbetaR+kbetaF+mE-1+3, p));

      modelcppadobj.derivative_coef();
      modelcppadobj.prepare_AIC_proposed();


      // IS_mat.block(0,0,kE+mE-1+kbetaR+kbetaF, kE+mE-1+kbetaR+kbetaF) = modelcppadobj.he_s_u_mat;
      // IS_mat.row(kE+mE-1+kbetaR+kbetaF) << modelcppadobj.he_alpha_f_log_theta_vec.transpose(), modelcppadobj.he_con_index_par_log_theta_vec.transpose(), modelcppadobj.he_betaR_log_theta_vec.transpose(), modelcppadobj.he_betaF_log_theta_vec.transpose(), modelcppadobj.he_log_theta_scalar;
      // IS_mat.col(kE+mE-1+kbetaR+kbetaF) = IS_mat.row(kE+mE-1+kbetaR+kbetaF).transpose();
      // I_mat = modelcppadobj.I_mat;
      Khat = modelcppadobj.Khat;

      // Mat Vbeta = IS_mat.ldlt().solve(Mat::Identity(kE+mE-1+kbetaR+kbetaF + 1, kE+mE-1+kbetaR+kbetaF + 1));
      // Mat Vbetahat = Vbeta * Khat * Vbeta;

      Eigen::LDLT<Mat> ldlt(Khat.selfadjointView<Eigen::Upper>()); // use Vbetahat = Vbeta * I_mat * Vbeta, instead of Vbeta itself. IMPORTANT for the proposed correction term!!!
      Mat U = ldlt.matrixU(); // Upper-triangular U so that Vbetahat = P^T * U^T * D * U * P
      // // R = sqrt(D) * U * P
      Vec d = ldlt.vectorD();
      Vec d_mask_sqrt(kE+mE-1+kbetaR+kbetaF + 1);
      for (int i = 0; i < kE+mE-1+kbetaR+kbetaF + 1; ++i) {
        d_mask_sqrt(i) = CppAD::CondExpGt(d(i), Scalar(4e-32), CppAD::sqrt(d(i)), Scalar(0));
      }
      Mat R = d_mask_sqrt.asDiagonal() * U * (ldlt.transpositionsP() * IdenMat); // Upper-triangular R so that Vbetahat = R^T * R; nrow(R) = kwopt+kE+mE-1+kbetaR+kbetaF + 1
      // ldlt.transpositionsP() is like an operator, and we need to multiply a matrix after it to get the permuted matrix.

      // Eigen::LLT<Mat> llt(Vbetahat); // use Vbetahat = Vbeta * I_mat * Vbeta, instead of Vbeta itself. IMPORTANT for the proposed correction term!!!
      // Mat R = llt.matrixU();  // Upper-triangular R so that Vbetahat = R^T * R; nrow(R) = kwopt+kE+mE-1+kbetaR+kbetaF + 1

      int k = 0;
      for (int j = 0; j < (kE+mE-1+kbetaR+kbetaF + 1); ++j) {
        for (int i = 0; i < (kE+mE-1+kbetaR+kbetaF + 1); ++i) {
          result[k++] = R(i, j);
        }
      }
      CppAD::ADFun<double> gr(at, result);

      return gr;
      // END forward 2

}




// [[Rcpp::export]]
List ConditionalAICsidDLNM(SEXP ptr, SEXP ptrcppad, Eigen::MatrixXd Vrho, bool verbose = false, bool marginal = false) {

  Rcpp::XPtr<Model> modelobj_ptr(ptr);
  Model& modelobj = *modelobj_ptr;

  Rcpp::XPtr<Modelcppad> modelcppadobj_ptr(ptrcppad);
  Modelcppad& modelcppadobj = *modelcppadobj_ptr;


  int kE = modelobj.kE;
  int kbetaR = modelobj.kbetaR;
  int kbetaF = modelobj.kbetaF;
  int mE = modelobj.mE;
  int p = modelobj.p;
  int n = modelobj.n;
  Eigen::VectorXd r = modelobj.r;


  Eigen::VectorXd alpha_f = modelobj.alpha_f;
  Eigen::VectorXd betaR = modelobj.betaR;
  Eigen::VectorXd betaF = modelobj.betaF;
  Eigen::VectorXd con_index_par = modelobj.con_index_par;

  double log_theta = modelobj.log_theta;
  double log_smoothing_f = modelobj.log_smoothing_f;
  double log_smoothing_w = modelobj.log_smoothing_w;
  Eigen::VectorXd logsmoothing = modelobj.logsmoothing;

  modelobj.prepare_AIC();
  // negative log likelihood without penalty
  double l = modelobj.NegLogL_l;

  // I
  Eigen::MatrixXd R_I;
  R_I = modelobj.I_mat;

  Eigen::MatrixXd IS_mat = modelobj.IS_mat;

  Eigen::MatrixXd Vbeta = IS_mat.ldlt().solve(Eigen::MatrixXd::Identity(IS_mat.rows(), IS_mat.cols())); // solve(IS_mat)
  Eigen::MatrixXd Khat = modelobj.Khat;


  Eigen::MatrixXd J = -1.0 * Vbeta * modelobj.he_s_par_u_log_theta_mat.transpose();
  Eigen::MatrixXd Smat = IS_mat - R_I; // penalty matrix

  // the widely used version of conditional AIC proposed by Hastie and Tibshirani (1990).
  // See Wood et al. 2016 JASA
  Eigen::MatrixXd mat_AIC = Vbeta * R_I;
  // double edf1 = (2 * mat_AIC - mat_AIC * mat_AIC).trace();
  double edf_conventional = mat_AIC.trace();
  double AIC_conventional = 2.0*l + 2.0*edf_conventional;

  // proposed conditional AIC
  Eigen::MatrixXd mat_cAIC = Khat * Vbeta;
  double edf_cAIC = mat_cAIC.trace();
  double AIC_cAIC = 2.0*l + 2.0*edf_cAIC;

  

  if(marginal) {
    // std::cout << "Start building tape" << std::endl;
    CppAD::ADFun<double> cppadgr = LAMLTapeAIC(modelobj, modelcppadobj);
    // cppadgr.optimize();
    // std::cout << "Finished building tape" << std::endl;

    // std::cout << "Start optimizing tape" << std::endl;
    // cppadgr.optimize();
    // std::cout << "Finished optimizing tape" << std::endl;

    int a = (kE+mE-1+kbetaR+kbetaF + 1) * (kE+mE-1+kbetaR+kbetaF + 2) / 2;
    int afull = (kE+mE-1+kbetaR+kbetaF + 1) * (kE+mE-1+kbetaR+kbetaF + 1);
    // int b = (kE+kbetaR+kbetaF + mE-1 + 3 + p);
    int b = (kE+kbetaR+kbetaF + mE);

    Eigen::MatrixXd dR(a, b);
    dR.setZero();

    // forward mode
    // 0-order forward: set at0 (evaluation point).
    Eigen::VectorXd x0(b);
    x0 << alpha_f, con_index_par, betaR, betaF, log_theta; //, log_smoothing_f, log_smoothing_w, logsmoothing;
    std::vector<double> at0(b);
    for (size_t i = 0; i < x0.size(); ++i) at0[i] = x0[i];
    // std::cout << "start 0-order forward" << std::endl;
    // cppadgr.Forward(0, at0); 
    // std::cout << "Finished 0-order forward" << std::endl;
    // std::cout << "start 1-order forward" << std::endl;

    // CppAD::ADFun<double> cppadgr = LAMLTapeAIC(modelobj, modelcppadobj);
    // cppadgr.optimize();
    // cppadgr.Forward(0, at0);
    // int n_theta = kE+mE-1+kbetaR+kbetaF + 1;
    // std::vector<double> dyn(2 * n_theta);
    
    // int k = 0;
    // for (int j = 0; j < n_theta; ++j) {
    //   for (int i = 0; i <= j; ++i) {
    //     // Eigen::VectorXd vi = Eigen::VectorXd::Unit(kE+mE-1+kbetaR+kbetaF + 1, i);
    //     // Eigen::VectorXd vj = Eigen::VectorXd::Unit(kE+mE-1+kbetaR+kbetaF + 1, j);
    //     std::fill(dyn.begin(), dyn.end(), 0.0);
    //     dyn[i] = 1.0;                    // vi
    //     dyn[j + n_theta] = 1.0;          // vj
    //     cppadgr.new_dynamic(dyn);
    //     std::vector<double> w(1, 1.0);
    //     auto g = cppadgr.Reverse(1, w);
    //     dR.row(k) = Eigen::Map<Eigen::VectorXd>(g.data(), b);
    //     k++;
    //     if(verbose) std::cout << "Computed dR row " << k << " / " << a << std::endl;
    //   }
    // }
    std::vector<double> input(b, 0.0);
    for (int j = 0; j < b; ++j) {
      // std::cout << j << std::endl;
      input[j] = 1.0;
      std::vector<double> dy = cppadgr.Forward(1, input);
      for (int i = 0; i < a; ++i) dR(i, j) = dy[i];
      input[j] = 0.0;
      if(verbose) std::cout << "Computed dR column " << j + 1 << " / " << b << std::endl;
    }
    // std::cout << "Finished 1-order forward" << std::endl;
    // end forward mode

    // free cppadgr.
    cppadgr = CppAD::ADFun<double>();

    // std::cout << "Start building tape proposed" << std::endl;
    CppAD::ADFun<double> cppadgr_proposed = LAMLTapeAICproposed(modelobj, modelcppadobj);
    // std::cout << "Finished building tape proposed" << std::endl;

    Eigen::MatrixXd dR_proposed(afull,b);
    // std::cout << "start 0-order forward proposed" << std::endl;
    // cppadgr_proposed.Forward(0, at0);
    // std::cout << "Finished 0-order forward proposed" << std::endl;
    // std::cout << "start 1-order forward proposed" << std::endl;
    std::vector<double> input_proposed(b, 0.0);
    for (int j = 0; j < b; ++j) {
      // std::cout << j << std::endl;
      input_proposed[j] = 1.0;
      std::vector<double> dy_proposed = cppadgr_proposed.Forward(1, input_proposed);
      for (int i = 0; i < afull; ++i) dR_proposed(i, j) = dy_proposed[i];
      input_proposed[j] = 0.0;
      if(verbose) std::cout << "Computed dR_proposed column " << j + 1 << " / " << b << std::endl;
    }
    // std::cout << "Finished 1-order forward proposed" << std::endl;
    // free cppadgr_proposed.
    cppadgr_proposed = CppAD::ADFun<double>();


    // Vpp in wood
    // // find dR / d rho, using implicit derivative
    // // Select the last (p+2) columns of dR
    // Eigen::MatrixXd dR_rho = dR.rightCols(p + 2); // dim: a*(p+2)
    // // Find the first (b - p - 2) columns of dR
    // Eigen::MatrixXd dR_beta = dR.leftCols(b - p - 2); // dim: a*(b-p-2)
    // // dim(J): (b-p-2) * (p+2)
    // Eigen::MatrixXd dRdrho = dR_rho + dR_beta * J; // dim: a*(p+2)

    Eigen::MatrixXd dRdrho = dR * J;
    Eigen::LLT<Eigen::MatrixXd> llt(Vbeta.selfadjointView<Eigen::Upper>());
    Eigen::MatrixXd R = llt.matrixU();  // Upper-triangular R so that Vbeta = R^T * R; nrow(R) = kwopt+kE+mE-1+kbetaR+kbetaF + 1


    // Eigen::LLT<Eigen::MatrixXd> llthat(Vbetahat.selfadjointView<Eigen::Upper>());
    // Eigen::MatrixXd Rhat = llthat.matrixU();  // Upper-triangular R so that Vbetahat = Rhat^T * Rhat; nrow(Rhat) = kwopt+kE+mE-1+kbetaR+kbetaF + 1
    // Eigen::MatrixXd Rhatinv = Rhat.ldlt().solve(Eigen::MatrixXd::Identity(Rhat.rows(), Rhat.cols())); // solve(Rhat)

    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigvec;
    // eigvec.compute(R_I, Eigen::ComputeEigenvectors);
    // Eigen::VectorXd eigvals = eigvec.eigenvalues().array();
    // Eigen::VectorXd invabseigvals(eigvals.size());
    // for (int i = 0; i < eigvals.size(); i++) invabseigvals(i) = 1. / abs(eigvals(i));
    // Eigen::MatrixXd Iinv = eigvec.eigenvectors() * (invabseigvals.asDiagonal()) * (eigvec.eigenvectors().transpose());
    // Eigen::MatrixXd Iinv = R_I.completeOrthogonalDecomposition().pseudoInverse();
    // A list with element dR / d rho

    std::vector<Eigen::MatrixXd> dRdrho_list;
    std::vector<Eigen::MatrixXd> dRhatdrho_list;
    std::vector<Eigen::MatrixXd> dVbetadrho_list;

    int begin = 0;
    for (int i = 0; i < p + 2; i++) {
      Eigen::VectorXd dRdrho_coli = dRdrho.col(i);
      Eigen::MatrixXd dIdrho(kE+mE-1+kbetaR+kbetaF + 1, kE+mE-1+kbetaR+kbetaF + 1);
      dIdrho.setZero();
      int k = 0;
      for (int jj = 0; jj < kE+mE-1+kbetaR+kbetaF + 1; ++jj) {
          for (int ii = 0; ii <= jj; ++ii) {
              dIdrho(ii, jj) = dRdrho_coli(k);
              k++;
          }
      }


      dIdrho = dIdrho.selfadjointView<Eigen::Upper>();
      Eigen::MatrixXd dISdrho = dIdrho;
      if (i == 0) {
        dISdrho.block(0, 0, kE, kE) += modelobj.smoothing_f * modelobj.getSf();
      }
      if (i == 1) {
        dISdrho.block(0, 0, kE, kE) += modelobj.smoothing_w * modelobj.getSw();
      }

      if (i >= 2) {
        int ki = static_cast<int>(r(i-2));
        dISdrho.block(kE+mE-1+begin,kE+mE-1+begin,ki,ki) += Smat.block(kE+mE-1+begin,kE+mE-1+begin,ki,ki);
        begin += ki;
      }



      Eigen::MatrixXd RinvtdVbetadrhoRinv = -1.0 * R * dISdrho * R.transpose();
      Eigen::MatrixXd RinvtdVbetadrhoRinvU = RinvtdVbetadrhoRinv.triangularView<Eigen::Upper>();
      RinvtdVbetadrhoRinvU.diagonal() *= 0.5;
      dRdrho_list.push_back(RinvtdVbetadrhoRinvU * R);

      Eigen::MatrixXd dVbetadrho = -1.0 * Vbeta * dISdrho * Vbeta;
      dVbetadrho_list.push_back(dVbetadrho);

      // // some components for dRhat / d rho
      // Eigen::MatrixXd dVbetadrho = -1.0 * Vbeta * dISdrho * Vbeta;
      // Eigen::MatrixXd dVbetahatdrho = dVbetadrho * R_I * Vbeta + Vbeta * dIdrho * Vbeta + Vbeta * R_I * dVbetadrho;

      // Eigen::MatrixXd RhatinvtdVbetahatdrhoRhatinv = -1.0 * Rhatinv.transpose() * dVbetahatdrho * Rhatinv;
      // // Step 1: solve Y = dVbetahatdrho * Rhat^{-1}   (right solve with Rhat)
      // // Eigen::MatrixXd Y = Rhat.triangularView<Eigen::Upper>().solve(dVbetahatdrho.transpose());
      // // Y.transposeInPlace();   // ensure Y = dV * Rhat^{-1}
      // // Step 2: solve M = Rhat^{-T} * Y   (left solve with Rhat^T)
      // // Eigen::MatrixXd RhatinvtdVbetahatdrhoRhatinv = Rhat.transpose().triangularView<Eigen::Lower>().solve(Y);

      // Eigen::MatrixXd RhatinvtdVbetahatdrhoRhatinvU = RhatinvtdVbetahatdrhoRhatinv.triangularView<Eigen::Upper>();
      // RhatinvtdVbetahatdrhoRhatinvU.diagonal() *= 0.5;


      // // some components for dRhat / d rho
      // Eigen::MatrixXd dIinvdrho = -1.0 * Iinv * dIdrho * Iinv;
      // // W = IS_mat * Iinv * IS_mat
      // Eigen::MatrixXd dWdrho = dISdrho * Iinv * IS_mat + IS_mat * dIinvdrho * IS_mat + IS_mat * Iinv * dISdrho;

      // Eigen::MatrixXd RhatinvtdVbetahatdrhoRhatinv = -1.0 * Rhat * dWdrho * Rhat.transpose();

      // Eigen::MatrixXd RhatinvtdVbetahatdrhoRhatinvU = RhatinvtdVbetahatdrhoRhatinv.triangularView<Eigen::Upper>();
      // RhatinvtdVbetahatdrhoRhatinvU.diagonal() *= 0.5;



      // dRhatdrho_list.push_back(RhatinvtdVbetahatdrhoRhatinvU * Rhat);

    }

    // calculate Vpp
    Eigen::MatrixXd Vpp(kE+mE-1+kbetaR+kbetaF + 1, kE+mE-1+kbetaR+kbetaF + 1);
    Vpp.setZero();
    for (size_t j = 0; j < kE+mE-1+kbetaR+kbetaF + 1; ++j)
    {
      for (size_t m = 0; m < kE+mE-1+kbetaR+kbetaF + 1; ++m)
      {
        double Vppjm = 0.0;
        for (size_t i = 0; i < kE+mE-1+kbetaR+kbetaF + 1; ++i)
        {
          for (size_t l = 0; l < p+2; ++l)
          {
            for (size_t k = 0; k < p+2; ++k)
            {
              Vppjm += dRdrho_list.at(k)(i,j) * Vrho(k,l) * dRdrho_list.at(l)(i,m);
            }
          }
        }
        Vpp(j, m) = Vppjm;
      }
    }
    // END vpp in wood



    // Vpp proposed

    Eigen::MatrixXd IdenMat = Eigen::MatrixXd::Identity(kE+mE-1+kbetaR+kbetaF + 1, kE+mE-1+kbetaR+kbetaF + 1);
    Eigen::LDLT<Eigen::MatrixXd> ldlt(Khat.selfadjointView<Eigen::Upper>());
    Eigen::MatrixXd U = ldlt.matrixU(); // Upper-triangular U so that Vbetahat = P^T * U^T * D * U * P
    // // R = sqrt(D) * U * P
    Eigen::VectorXd d = ldlt.vectorD();
    Eigen::VectorXd d_mask_sqrt(kE+mE-1+kbetaR+kbetaF + 1);
    for (int i = 0; i < kE+mE-1+kbetaR+kbetaF + 1; ++i) {
      if(d(i) > 4e-32) {
        d_mask_sqrt(i) = std::sqrt(d(i));
      } else {
        d_mask_sqrt(i) = 0.0;
      }
    }
    Eigen::MatrixXd RKhat = d_mask_sqrt.asDiagonal() * U * (ldlt.transpositionsP() * IdenMat);



    // calculate Vpp proposed

    // // Select the last (p+2) columns of dR
    // Eigen::MatrixXd dR_rho_proposed = dR_proposed.rightCols(p + 2); // dim: a*(p+2)
    // // Find the first (b - p - 2) columns of dR
    // Eigen::MatrixXd dR_beta_proposed = dR_proposed.leftCols(b - p - 2); // dim: a*(b-p-2)
    // // dim(J): (b-p-2) * (p+2)
    // Eigen::MatrixXd dRdrho_proposed = dR_rho_proposed + dR_beta_proposed * J; // dim: a*(p+2)

    Eigen::MatrixXd dRdrho_proposed = dR_proposed * J;

    // A list with element dR / d rho
    std::vector<Eigen::MatrixXd> dRdrho_proposed_list;
    for (int i = 0; i < p + 2; i++) {
      Eigen::VectorXd dRdrho_proposed_coli = dRdrho_proposed.col(i);
      Eigen::MatrixXd dRdrho_proposed_coli_mat(kE+mE-1+kbetaR+kbetaF + 1, kE+mE-1+kbetaR+kbetaF + 1);
      dRdrho_proposed_coli_mat.setZero();
      int k = 0;
      for (int jj = 0; jj < kE+mE-1+kbetaR+kbetaF + 1; ++jj) {
          for (int ii = 0; ii < kE+mE-1+kbetaR+kbetaF + 1; ++ii) {
              dRdrho_proposed_coli_mat(ii, jj) = dRdrho_proposed_coli(k++);
          }
      }
      dRdrho_proposed_list.push_back(dRdrho_proposed_coli_mat * Vbeta + RKhat * dVbetadrho_list.at(i));
    }

    // calculate Vpp
    Eigen::MatrixXd Vpp_proposed(kE+mE-1+kbetaR+kbetaF + 1, kE+mE-1+kbetaR+kbetaF + 1);
    Vpp_proposed.setZero();
    for (size_t j = 0; j < kE+mE-1+kbetaR+kbetaF + 1; ++j)
    {
      for (size_t m = 0; m < kE+mE-1+kbetaR+kbetaF + 1; ++m)
      {
        double Vppjm = 0.0;
        for (size_t i = 0; i < kE+mE-1+kbetaR+kbetaF + 1; ++i)
        {
          for (size_t l = 0; l < p+2; ++l)
          {
            for (size_t k = 0; k < p+2; ++k)
            {
              Vppjm += dRdrho_proposed_list.at(k)(i,j) * Vrho(k,l) * dRdrho_proposed_list.at(l)(i,m);
            }
          }
        }
        Vpp_proposed(j, m) = Vppjm;
      }
    }
    // Vpp proposed END


    // // Vpp proposed new
    // Eigen::MatrixXd Vpp_proposed(kE+mE-1+kbetaR+kbetaF + 1, kE+mE-1+kbetaR+kbetaF + 1);
    // Vpp_proposed.setZero();
    // for (size_t j = 0; j < kE+mE-1+kbetaR+kbetaF + 1; ++j)
    // {
    //   for (size_t m = 0; m < kE+mE-1+kbetaR+kbetaF + 1; ++m)
    //   {
    //     double Vppjm = 0.0;
    //     for (size_t i = 0; i < kE+mE-1+kbetaR+kbetaF + 1; ++i)
    //     {
    //       for (size_t l = 0; l < p+2; ++l)
    //       {
    //         for (size_t k = 0; k < p+2; ++k)
    //         {
    //           Vppjm += dRhatdrho_list.at(k)(i,j) * Vrho(k,l) * dRhatdrho_list.at(l)(i,m);
    //         }
    //       }
    //     }
    //     Vpp_proposed(j, m) = Vppjm;
    //   }
    // }
    // // END Vpp proposed new
  

    // Wood et al. 2016 JASA
    Eigen::MatrixXd Vbetap = Vbeta + J * Vrho * J.transpose() + Vpp;
    Eigen::MatrixXd mat_AIC_woodetal = Vbetap * R_I;
    double edf_woodetal = mat_AIC_woodetal.trace();
    double AIC_woodetal = 2.0*l + 2.0*edf_woodetal;

    // Wood et al with corrected Vpp
    Eigen::MatrixXd mat_AIC_woodetal_correct = Vbeta * Khat + (J * Vrho * J.transpose() + Vpp_proposed) * R_I;
    double edf_woodetal_correct = mat_AIC_woodetal_correct.trace();
    double AIC_woodetal_correct = 2.0*l + 2.0*edf_woodetal_correct;


    // proposed AIC
    // Eigen::MatrixXd mat_corrected = Vbetap * Smat + Vbeta * R_I - Eigen::MatrixXd::Identity(IS_mat.rows(), IS_mat.cols());
    Eigen::MatrixXd mat_corrected = (J * Vrho * J.transpose() + Vpp_proposed) * Smat;
    Eigen::MatrixXd mat_proposed = mat_AIC_woodetal_correct + mat_corrected;
    double edf_proposed = mat_proposed.trace();
    double AIC_proposed = 2.0*l + 2.0*edf_proposed;

    return List::create(Named("l") = -1.0*l,
                        Named("edf_conventional") = edf_conventional,
                        Named("AIC_conventional") = AIC_conventional,
                        Named("edf_cAIC") = edf_cAIC,
                        Named("AIC_cAIC") = AIC_cAIC,
                        // wood et al 2016
                        Named("mat_AIC_woodetal") = mat_AIC_woodetal,
                        Named("edf_woodetal") = edf_woodetal,
                        Named("AIC_woodetal") = AIC_woodetal,
                        // wood et al 2016 with corrected Vpp
                        Named("mat_AIC_woodetal_correct") = mat_AIC_woodetal_correct,
                        Named("edf_woodetal_correct") = edf_woodetal_correct,
                        Named("AIC_woodetal_correct") = AIC_woodetal_correct,
                        // proposed
                        Named("mat_proposed") = mat_proposed,
                        Named("edf_proposed") = edf_proposed,
                        Named("AIC_proposed") = AIC_proposed
                        );
  } else {
    return List::create(Named("l") = -1.0*l,
                        Named("edf_conventional") = edf_conventional,
                        Named("AIC_conventional") = AIC_conventional,
                        Named("edf_cAIC") = edf_cAIC,
                        Named("AIC_cAIC") = AIC_cAIC
                        );
  }
}



// [[Rcpp::export]]
List NCVsidDLNM(SEXP ptr, const List nei_list, bool verbose = false, int nthreads = 1) {
  Rcpp::XPtr<Model> modelobj_ptr(ptr);
  Model& modelobj = *modelobj_ptr;


  int kE = modelobj.kE;
  int kbetaR = modelobj.kbetaR;
  int kbetaF = modelobj.kbetaF;
  int mE = modelobj.mE;
  int p = modelobj.p;
  int n = modelobj.n;
  Eigen::VectorXd r = modelobj.r;

  Eigen::VectorXd alpha_f = modelobj.alpha_f;
  Eigen::VectorXd betaR = modelobj.betaR;
  Eigen::VectorXd betaF = modelobj.betaF;
  Eigen::VectorXd con_index_par = modelobj.con_index_par;
  double log_theta = modelobj.log_theta;




  Eigen::VectorXd beta_mod(kE+mE-1+kbetaR+kbetaF+1);
  beta_mod << alpha_f, con_index_par, betaR, betaF, log_theta;
  Eigen::MatrixXd beta_nei(beta_mod.size(), nei_list.size());

  modelobj.derivative_coef();
  modelobj.derivative_he();
  modelobj.derivative_full();
  modelobj.prepare_AIC(); // for gunpen_nei
  Eigen::MatrixXd Hpen = modelobj.IS_mat;

  Eigen::MatrixXd Kleft = modelobj.Kleft;
  Eigen::VectorXd gnei(Kleft.rows());


  Eigen::VectorXd nei_vec;
  Eigen::MatrixXd Hunpen_nei;
  Eigen::MatrixXd Hlambdanei(kE+mE-1+kbetaR+kbetaF+1, kE+mE-1+kbetaR+kbetaF+1); // Hpen - Hunpen_nei

  Eigen::VectorXd Dnei(nei_list.size());
  Eigen::VectorXd Pnei(nei_list.size());
  Eigen::VectorXd BMAP(nei_list.size());


  // Eigen::VectorXd Dfull(nei_list.size());
  for (size_t i = 0; i < nei_list.size(); i++) {
    // Dfull(i) = modelobj.NegativeLogLikelihood_l_i(i);

    // std::cout << "NCV for neighbor " << i+1 << " / " << nei_list.size() << std::endl;
    nei_vec = as<Eigen::VectorXd>(nei_list[i]);
    // std::cout << "nei_vec" << nei_vec.transpose() << std::endl;

    // compute Hunpen_nei
    Hunpen_nei = modelobj.prepare_NCV(nei_vec);
    Hlambdanei = Hpen - Hunpen_nei;

    // compute gunpen_nei from K
    gnei.setZero();
    for (size_t j = 0; j < nei_vec.size(); j++) {
      int j_int = static_cast<int>(nei_vec(j)) - 1; // R to C++ index
      gnei += Kleft.col(j_int); // col(j_int);
    }

    beta_nei.col(i) = beta_mod - Hlambdanei.completeOrthogonalDecomposition().solve(gnei);


    Eigen::VectorXd beta_i = beta_nei.col(i);

    // for NCV loss
    Dnei(i) = modelobj.get_D_i(beta_i.segment(0, kE), beta_i.segment(kE, mE-1), beta_i.segment(kE+mE-1, kbetaR),
                     beta_i.segment(kE+mE-1+kbetaR, kbetaF), 
                     beta_i(kE+mE-1+kbetaR+kbetaF),
                     static_cast<int>(i));
  }

  int MCR = 1000;
  Eigen::VectorXd p_i_vec(MCR);

  // // compute BMA \int p(y_t | beta) p(beta | Mk) d beta where p(beta | Mk) is determined the penalty term
  // Eigen::MatrixXd I_mat = modelobj.I_mat.block(0, 0, kE+mE-1+kbetaR+kbetaF, kE+mE-1+kbetaR+kbetaF);
  // Eigen::MatrixXd S_mat = Hpen - I_mat; // S_mat inv is the covariance matrix prior

  // p_i_vec.setZero();
  // Eigen::MatrixXd S_mat_L_inv;
  // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigvecS(S_mat,true); // Both values and vectors
  // Eigen::VectorXd eigvalsS = eigvecS.eigenvalues().array();
  // Eigen::VectorXd invabseigvalsS(eigvalsS.size());
  // for (int ii = 0; ii < eigvalsS.size(); ii++) invabseigvalsS(ii) = 1. / max(abs(eigvalsS(ii)), 1e-3);
  // S_mat_L_inv = eigvecS.eigenvectors() * (invabseigvalsS.cwiseSqrt().asDiagonal());

  // for (size_t i = 0; i < nei_list.size(); i++) {
  //   for(int r = 0; r < MCR; r++) {
  //     // Jointly sample
  //     Eigen::VectorXd zjoint(kE+mE-1+kbetaR+kbetaF);
  //     p_i_vec.setZero();

  //     for (int j = 0; j < kE+mE-1+kbetaR+kbetaF; j++) {
  //       zjoint(j) = dist(gen);
  //     }

  //     Eigen::VectorXd samplejoint = S_mat_L_inv * zjoint; // mean zero and covariance S_mat_inv
  //     // get alpha_f
  //     Eigen::VectorXd R_alpha_f_sample = samplejoint.segment(0, kE);


  //     // get index par
  //     Eigen::VectorXd R_con_index_par_sample = samplejoint.segment(kE, mE-1);
  //     Eigen::VectorXd R_con_index_par_long(mE); // con_index_par_long = c(1,con_index_par)
  //     R_con_index_par_long(0) = 1.0;
  //     for (int j = 0; j < (mE - 1); j++) {
  //       R_con_index_par_long(j + 1) = R_con_index_par_sample(j);
  //     }
  //     Eigen::VectorXd R_index_par_sample = modelobj.Bindex * R_con_index_par_long/sqrt(R_con_index_par_long.dot(modelobj.BtBindex * R_con_index_par_long));

  //     // get betaR
  //     Eigen::VectorXd R_betaR_sample = samplejoint.segment(kE+mE-1, kbetaR);
  //     // get betaF
  //     Eigen::VectorXd R_betaF_sample = samplejoint.segment(kE+mE-1+kbetaR, kbetaF);

  //     double p_i = modelobj.get_p_i(R_alpha_f_sample, R_con_index_par_sample, R_betaR_sample, R_betaF_sample, static_cast<int>(i));

  //     p_i_vec(r) = p_i;
  //   }
  //   BMAP(i) = p_i_vec.mean();
  // }
  // double logPy = BMAP.array().log().sum();




  // compute predictive density
  // MCR = 1000;
  // p_i_vec.resize(MCR);
  // Eigen::MatrixXd PneiMat(MCR, nei_list.size());

  // single core version
  if (nthreads == 1) {

    std::mt19937 gen(123);
    std::normal_distribution<> dist(0, 1);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigvec(Hlambdanei); // Both values and vectors
    Eigen::LLT<Eigen::MatrixXd> cholSolver(Hlambdanei);
    p_i_vec.setZero();
    Eigen::MatrixXd R_he_u_L_inv;
    for (size_t i = 0; i < nei_list.size(); i++) {
      Eigen::VectorXd beta_i = beta_nei.col(i);
      // modelobj.setAlphaF(beta_i.segment(0, kE));
      modelobj.setAlphaFBetaRBetaFCon_Index_Par(beta_i.segment(0, kE), beta_i.segment(kE+mE-1, kbetaR), beta_i.segment(kE+mE-1+kbetaR, kbetaF), beta_i.segment(kE, mE-1), beta_i(kE+mE-1+kbetaR+kbetaF));
      // modelobj.setCon_Index_Par(beta_i.segment(kE, mE-1));
      // modelobj.setBetaR(beta_i.segment(kE+mE-1, kbetaR));
      // modelobj.setBetaF(beta_i.segment(kE+mE-1+kbetaR, kbetaF));

      modelobj.derivative_coef();
      modelobj.derivative_he();
      modelobj.derivative_full();
      modelobj.prepare_AIC();

      nei_vec = as<Eigen::VectorXd>(nei_list[i]);
      Hlambdanei = modelobj.IS_mat - modelobj.prepare_NCV(nei_vec);

      cholSolver.compute(Hlambdanei);
      if(cholSolver.info()!=Eigen::Success) {
        eigvec.compute(Hlambdanei); // Compute eigenvalues and vectors
        Eigen::VectorXd eigvals = eigvec.eigenvalues().array();
        Eigen::VectorXd invabseigvals(eigvals.size());
        for (int ii = 0; ii < eigvals.size(); ii++) invabseigvals(ii) = 1. / max(abs(eigvals(ii)), 1e-3);
        R_he_u_L_inv = eigvec.eigenvectors() * (invabseigvals.cwiseSqrt().asDiagonal());
        // if(verbose) {
        //   std::cout << "Warning: HLambdanei is not positive definite for neighbor " << i + 1 << " . Using eigen decomposition to compute the inverse of Cholesky factor." << std::endl;
        // }
      } else {
        Eigen::MatrixXd chol_L = cholSolver.matrixL();
        R_he_u_L_inv = invertL(chol_L).transpose();
      }

      Eigen::VectorXd zjoint(kE+kbetaR+kbetaF + mE);



      p_i_vec.setZero();
      for(int r = 0; r < MCR; r++) {
        // Jointly sample
        for (int j = 0; j < kE+kbetaR+kbetaF + mE; j++) {
          zjoint(j) = dist(gen);
        }
        Eigen::VectorXd samplejoint = beta_i + R_he_u_L_inv * zjoint;
        // get alpha_f
        Eigen::VectorXd R_alpha_f_sample = samplejoint.segment(0, kE);


        // get index par
        Eigen::VectorXd R_con_index_par_sample = samplejoint.segment(kE, mE-1);
        Eigen::VectorXd R_con_index_par_long(mE); // con_index_par_long = c(1,con_index_par)
        R_con_index_par_long(0) = 1.0;
        for (int j = 0; j < (mE - 1); j++) {
          R_con_index_par_long(j + 1) = R_con_index_par_sample(j);
        }
        Eigen::VectorXd R_index_par_sample = modelobj.Bindex * R_con_index_par_long/sqrt(R_con_index_par_long.dot(modelobj.BtBindex * R_con_index_par_long));

        // get betaR
        Eigen::VectorXd R_betaR_sample = samplejoint.segment(kE+mE-1, kbetaR);
        // get betaF
        Eigen::VectorXd R_betaF_sample = samplejoint.segment(kE+mE-1+kbetaR, kbetaF);

        double R_log_theta_sample = samplejoint(kE+mE-1+kbetaR+kbetaF);

        double p_i = modelobj.get_p_i(R_alpha_f_sample, R_con_index_par_sample, R_betaR_sample, R_betaF_sample, R_log_theta_sample, static_cast<int>(i));

        p_i_vec(r) = p_i;

      }
      // PneiMat.col(i) = p_i_vec;
      Pnei(i) = p_i_vec.mean();

      // print if is multiple of 100
      if(verbose) {
        if ((i + 1) % 100 == 0) {
          std::cout << "calculate NCV loss and predictive density: for neighbor " << i + 1 << " / " << nei_list.size() << std::endl;
        }
      }
    }
  }
  // end single core version

  // multi-thread version

  if(nthreads > 1) {

    // openMP version
    const size_t N = nei_list.size();
    std::vector<Eigen::VectorXd> nei_vecs(N);
    for (size_t i = 0; i < N; ++i) nei_vecs[i] = as<Eigen::VectorXd>(nei_list[i]);

    #ifdef _OPENMP
      omp_set_num_threads(nthreads);
    #endif

    #ifdef _OPENMP
      if(verbose) {
        std::cout << "OpenMP is ON" << std::endl;
        std::cout << "Using " << nthreads << " threads for NCV loss and predictive density calculation." << std::endl;
      }
    #else
      std::cout << "Warnings: OpenMP is OFF. Install OpenMP to enable parallel processing." << std::endl;
    #endif

    Eigen::setNbThreads(1);
    
    #pragma omp parallel
    {
    
    Model local_model(modelobj);  // one per thread
    Eigen::VectorXd zjoint(kE+kbetaR+kbetaF + mE);
    Eigen::VectorXd local_p_i_vec(MCR);
    
    #pragma omp for schedule(static)
    for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(N); ++i) {
      // Model local_model(modelobj);


      std::mt19937 gen(123);
      std::normal_distribution<> dist(0, 1);



      const Eigen::Ref<const Eigen::VectorXd> beta_i = beta_nei.col(i);

      local_model.setAlphaFBetaRBetaFCon_Index_Par(beta_i.segment(0, kE), beta_i.segment(kE+mE-1, kbetaR), beta_i.segment(kE+mE-1+kbetaR, kbetaF), beta_i.segment(kE, mE-1), beta_i(kE+mE-1+kbetaR+kbetaF));
      local_model.derivative_coef();
      local_model.derivative_he();
      local_model.derivative_full();
      local_model.prepare_AIC();

      Eigen::MatrixXd local_Hlambdanei = local_model.IS_mat - local_model.prepare_NCV(nei_vecs[i]);
      Eigen::MatrixXd R_he_u_L_inv;

      Eigen::LLT<Eigen::MatrixXd> cholSolver(local_Hlambdanei);
      if(cholSolver.info()!=Eigen::Success) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigvec(local_Hlambdanei); // Both values and vectors
        Eigen::VectorXd eigvals = eigvec.eigenvalues().array();
        Eigen::VectorXd invabseigvals(eigvals.size());
        for (int ii = 0; ii < eigvals.size(); ii++) invabseigvals(ii) = 1. / max(abs(eigvals(ii)), 1e-3);
        R_he_u_L_inv = eigvec.eigenvectors() * (invabseigvals.cwiseSqrt().asDiagonal());
        // if(verbose) {
        //   std::cout << "Warning: local_Hlambdanei is not positive definite for neighbor " << i + 1 << " . Using eigen decomposition to compute the inverse of Cholesky factor." << std::endl;
        // }
      } else {
        Eigen::MatrixXd chol_L = cholSolver.matrixL();
        R_he_u_L_inv = invertL(chol_L).transpose();
      }

      // Eigen::VectorXd zjoint(kE+kbetaR+kbetaF + mE);
      // Eigen::VectorXd local_p_i_vec(MCR);
      local_p_i_vec.setZero();

      for(int r = 0; r < MCR; r++) {
        // Jointly sample
        for (int j = 0; j < kE+kbetaR+kbetaF + mE; j++) {
          zjoint(j) = dist(gen);
        }
        Eigen::VectorXd samplejoint = beta_i + R_he_u_L_inv * zjoint;
        // get alpha_f
        Eigen::VectorXd R_alpha_f_sample = samplejoint.segment(0, kE);


        // get index par
        Eigen::VectorXd R_con_index_par_sample = samplejoint.segment(kE, mE-1);
        Eigen::VectorXd R_con_index_par_long(mE); // con_index_par_long = c(1,con_index_par)
        R_con_index_par_long(0) = 1.0;
        for (int j = 0; j < (mE - 1); j++) {
          R_con_index_par_long(j + 1) = R_con_index_par_sample(j);
        }
        Eigen::VectorXd R_index_par_sample = local_model.Bindex * R_con_index_par_long/sqrt(R_con_index_par_long.dot(local_model.BtBindex * R_con_index_par_long));

        // get betaR
        Eigen::VectorXd R_betaR_sample = samplejoint.segment(kE+mE-1, kbetaR);
        // get betaF
        Eigen::VectorXd R_betaF_sample = samplejoint.segment(kE+mE-1+kbetaR, kbetaF);

        double R_log_theta_sample = samplejoint(kE+mE-1+kbetaR+kbetaF);

        double p_i = local_model.get_p_i(R_alpha_f_sample, R_con_index_par_sample, R_betaR_sample, R_betaF_sample, R_log_theta_sample, static_cast<int>(i));

        local_p_i_vec(r) = p_i;
      }
      Pnei(i) = local_p_i_vec.mean();

      // print if is multiple of 100
      if(verbose) {
        if ((i + 1) % 100 == 0) {
          std::cout << "calculate NCV loss and predictive density: for neighbor " << i + 1 << " / " << nei_list.size() << std::endl;
        }
      }
    }

    }

  }

  if(verbose) {
      if ((nei_list.size() + 1) % 100 != 0) {
        std::cout << "calculate NCV loss and predictive density: for neighbor " << nei_list.size() << " / " << nei_list.size() << std::endl;
      }
    }


  return List::create(Named("beta_nei") = beta_nei,
                      Named("beta_mod") = beta_mod,
                      Named("Dnei") = Dnei,
                      Named("Pnei") = Pnei
                      );

}
