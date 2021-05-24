#include <RcppArmadillo.h>
#include <vector>
#include <limits>
#include <algorithm>
#include <numeric>      // std::iota
#include <cmath>        // std::abs

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

struct cv_aux_out{
  arma::mat MSFEs;
  arma::mat sparsity;
  arma::cube Phi;
  arma::cube varxB;
};

struct prox_out{
  arma::mat PHI;
  arma::mat B;
};

struct hvarx{
  arma::mat PHI;
  arma::mat B;
  arma::mat phi0;
};

// Useful 2-norm function
// [[Rcpp::export]]
double norm2(NumericVector x){
  arma::vec xx = x;
  double g=arma::norm(xx,2);
  return (as<double>(wrap(g)));
}

uvec ind(int n2,int m){

  std::vector<int> subs;

  for(int i =0 ; i<n2;++i)

  {

    subs.push_back(i);

  }

  subs.erase(subs.begin()+m);

  return(conv_to<uvec>::from(subs));
}

uvec bbsubs(int j,int k,int p)
{
  uvec bb(p);
  bb(0)=j;
  for(int i=1;i<p;++i)
  {
    bb(i)=j+k*(i);

  }
  return(bb);

}


uvec vsubscppelem(int p,int pmax)
{
  uvec vs(pmax-p+1);
  for(int i=pmax;i>=p;--i)
  {
    vs(i-p)=i-1;
  }
  return(vs);
}


rowvec proxcppelem(arma::colvec v2,int L,double lambda,uvec res1, arma::colvec w)
{
  arma::colvec r =v2;
  for(int i=(L-1); i>=0;--i)

  {


    uvec res=vsubscppelem(i+1,L);



    if(norm(r(res)/(lambda*w(i)),"fro")<1+1e-8)
    {
      r(res)=zeros(res.n_elem);
    }
    else{
      r(res)=r(res)-lambda*w(i)*r(res)/(norm(r(res),"fro"));
    }

  }


  return(trans(r));


}

arma::colvec prox2HVARnew(arma::colvec v, double lambda, int k,int p)
{
  uvec res1=ind(p,0);
  arma::colvec w(p);
  w.ones();
  arma::colvec v2(v.n_elem);
  arma::rowvec v3(p);
  for(int i=0;i<k;++i)
  {
    uvec bb=bbsubs(i,k,p);
    arma::colvec v1=v(bb);
    v3=proxcppelem(v1,p,lambda,res1,w);
    v2(bb)=v3;
  }
  return(v2);
}

// Prox function for HVAR
// [[Rcpp::export]]
arma::rowvec prox2HVAR(arma::colvec v, double lambda, int k,int p)
{
  uvec res1=ind(p,0);
  arma::colvec w(p);
  w.ones();
  arma::rowvec v2(v.n_elem);
  arma::rowvec v3(p);
  for(int i=0;i<k;++i)
  {
    uvec bb=bbsubs(i,k,p);
    arma::colvec v1=v(bb);
    v3=proxcppelem(v1,p,lambda,res1,w);
    v2(bb)=v3;
  }
  return(v2);
}

arma::rowvec prox2(arma::colvec v,double lambda, int k,int p,uvec res1, arma::colvec w)
{
  arma::rowvec v2(v.n_elem);
  arma::rowvec v3(p);
  for(int i=0;i<k;++i)
  {
    uvec bb=bbsubs(i,k,p);
    arma::colvec v1=v(bb);
    v3=proxcppelem(v1,p,lambda,res1,w);
    v2(bb)=v3;
  }
  return(v2);
}

// [[Rcpp::export]]
arma::mat FistaElem(const arma::mat& Y,const arma::mat& Z, arma::mat phi, const int p,const int k,double lambda, const double eps,const double tk)
{
  double j=1;
  arma::mat phiFin=phi;
  arma::rowvec phiR=phi.row(0);
  arma::rowvec phiOLD=phiR;
  arma::rowvec phiOLDOLD=phiOLD;
  arma::rowvec v=phiOLD;
  uvec res1=ind(p,0);
  arma::colvec w(p);
  w.ones();


  for(int i=0;i<k;++i)
  {
    j=1;
    double thresh=10*eps;
    phiR=phi.row(i);
    phiOLD=phiR;
    phiOLDOLD=phiOLD;
    v=phiR;
    while(thresh>eps)
    { 
      v=phiOLD+((j-2)/(j+1))*(phiOLD-phiOLDOLD);

      phiR=prox2(vectorise(v)+tk*vectorise((trans(Y.col(i))-v*Z)*trans(Z)),tk*lambda,k,p,res1,w);

      thresh=max(abs(phiR-v));
      phiOLDOLD=phiOLD;
      phiOLD=phiR;
      j+=1;

    }
    phiFin.row(i)=phiR;
  }
  return(phiFin);




}

// Lamba loop
// [[Rcpp::export]]
arma::cube gamloopElem(NumericVector beta_, const arma::mat& Y,const arma::mat& Z, arma::colvec gammgrid, const double eps,const arma::colvec YMean2, const arma::colvec ZMean2, arma::mat B1, const int k, const int p){

  arma::mat B1F2=B1;

  vec eigval;
  arma::mat eigvec;
  const mat Zt=Z*trans(Z);
  eig_sym(eigval, eigvec, Zt);

  double tk=1/max(eigval);


  IntegerVector dims=beta_.attr("dim");


  const int ngridpts=dims[2];
  arma::cube bcube(beta_.begin(),dims[0],dims[1],ngridpts,false);
  arma::cube bcube2(dims[0],dims[1]+1,ngridpts);
  bcube2.fill(0);
  arma::colvec nu=zeros<colvec>(k);

  int i;

  for (i=0; i<ngridpts;++i) {


    B1F2=bcube.slice(i);
    B1 = FistaElem(Y,Z,B1F2,p,k,gammgrid[i],eps,tk);

    nu = YMean2 - B1 *ZMean2;



    bcube2.slice(i) = mat(join_horiz(nu, B1));
  }

  return(bcube2);
}


// Soft thresholding
// [[Rcpp::export]]
double ST1a(double z,double gam){

  if(z>0 && gam<fabs(z)) return(z-gam);
  if(z<0 && gam<fabs(z)) return(z+gam);
  if(gam>=fabs(z)) return(0);
  else return(0);

}


// [[Rcpp::export]]
arma::colvec ST3a(arma::colvec z ,double gam)
{

  int n=z.size();

  arma::colvec z1(n);
  for( int i=0; i<n;++i)
  {
    double z11=z(i);
    z1(i)=ST1a(z11,gam);

  }

  return(z1);
}

// Lasso Fista Function
arma::mat FistaLV(const arma::mat& Y, const arma::mat& Z, arma::mat& B, const double gam, const double eps, double tk, int k,int p)
{
  B=trans(B);
  arma::colvec B1=B.col(0);

  double j = 1;

  for( int i =0; i<k; ++i)
  {
    B1=B.col(i);
    arma::colvec BOLD=B.col(i);
    arma::colvec BOLDOLD=BOLD;
    double thresh=10*eps;
    j=1;

    while(thresh>eps)
    {

      arma::colvec v=BOLD+((j-2)/(j+1))*(BOLD-BOLDOLD);

      B1=ST3a(vectorise(v)+tk*vectorise((trans(Y.col(i))-trans(v)*Z)*trans(Z)),gam*tk);
      thresh=max(abs(B1-v));
      BOLDOLD=BOLD;
      BOLD=B1;
      j+=1;
      if(j>10000){
        break;
      }


    }

    B.col(i)=B1;

  }


  B=trans(B);

  return(B);

}


// [[Rcpp::export]]
arma::cube gamloopFista(NumericVector beta_, const arma::mat& Y,const arma::mat& Z,const  arma::colvec gammgrid, const double eps,const arma::colvec& YMean2, const arma::colvec& ZMean2, arma::mat& B1, int k, int p,double tk, int k1,int s){

  arma::mat b2=B1;
  arma::mat B1F2=B1;
  // const int ngridpts=gammgrid.size();
  IntegerVector dims=beta_.attr("dim");

  arma::cube bcube(beta_.begin(),dims[0],dims[1],dims[2],false);
  arma::cube bcube2(dims[0],dims[1]+1,dims[2]);
  bcube2.fill(0);

  arma::colvec nu=zeros<colvec>(dims[0]);
  double gam =0;

  int i;
  //loop through candidate lambda values
  for (i=0; i<dims[2];++i) {
    gam=gammgrid[i];

    arma::mat B1F2=bcube.slice(i);
    B1 = FistaLV(Y,Z,B1F2,gam,eps,tk,k1,p);

    nu = YMean2 - B1 *ZMean2;
    bcube2.slice(i) = mat(join_horiz(nu, B1));
  }

  return(bcube2);
}


// Lamba loop
// [[Rcpp::export]]
arma::cube gamloopElem2(arma::cube &bcube, const arma::mat& Y,const arma::mat& Z, arma::colvec gammgrid, const double eps,
                        const arma::colvec YMean2, const arma::colvec ZMean2, arma::mat &B1, const int k, const int p, const double tk,
                        const int flag_restart_opt = 1){

  arma::mat B1F2 = B1;
  // IntegerVector dims=beta_.attr("dim");
  int ngridpts = bcube.n_slices;
  arma::cube bcube2(k, k*p + 1, ngridpts);
  bcube2.fill(0);
  colvec nu=zeros<colvec>(k);

  int i;
  B1 = bcube.slice(0);
  for (i=0; i<ngridpts;++i) {
    B1F2 = (flag_restart_opt == 1) ? B1 : bcube.slice(i); // using previous iter value if fresh opt else use given starting value
    B1 = FistaElem(Y, Z, B1F2, p, k, gammgrid[i], eps, tk);
    nu = YMean2 - B1*ZMean2;
    bcube2.slice(i) = mat(join_horiz(nu, B1));
  }

  return(bcube2);
}


// [[Rcpp::export]]
arma::cube gamloopFista2(arma::cube bcube, const arma::mat& Y,const arma::mat& Z,const  arma::colvec gammgrid, const double eps,
                         const arma::colvec& YMean2, const arma::colvec& ZMean2, arma::mat& B1, int k, int p, double tk){

  arma::mat b2=B1;
  arma::mat B1F2=B1;
  int ngridpts = bcube.n_slices;
  arma::cube bcube2(k, k*p + 1, ngridpts);
  bcube2.fill(0);
  arma::colvec nu=zeros<colvec>(k);

  // const int ngridpts=gammgrid.size();
  // IntegerVector dims=beta_.attr("dim");
  // arma::cube bcube(beta_.begin(),dims[0],dims[1],dims[2],false);
  // arma::cube bcube2(dims[0],dims[1]+1,dims[2]);


  double gam =0;

  int i;
  //loop through candidate lambda values
  for (i=0; i<ngridpts;++i) {
    gam=gammgrid[i];
    arma::mat B1F2 = bcube.slice(i);
    B1 = FistaLV(Y, Z, B1F2, gam, eps, tk, k, p);
    nu = YMean2 - B1 *ZMean2;
    bcube2.slice(i) = mat(join_horiz(nu, B1));
  }

  return(bcube2);
}


// [[Rcpp::export]]
arma::cube lassoVARFistcpp(const arma::cube& beta, const arma::mat& trainY, const arma::mat& trainZ, const arma::colvec& lambda,
                           const double& tol, const int& p){

  int n = trainY.n_rows;
  int k = trainY.n_cols;
  int g = beta.n_slices;
  arma::mat YMean = mean(trainY);
  arma::mat ZMean = mean(trainZ.t());

  arma::mat Y = zeros(n, k);
  arma::mat trainZt = trainZ.t();
  arma::mat Z = zeros(n, k*p);
  for (int i = 0; i < n; i++) {
    Y.row(i) = trainY.row(i)  - YMean ;
    Z.row(i) = trainZt.row(i)  - ZMean ;
  }
  Z = Z.t();

  // Step size
  vec eigval;
  arma::mat eigvec;
  const arma::mat Zt=Z*trans(Z);
  eig_sym(eigval, eigvec, Zt);
  double tk=1/max(eigval);

  arma::mat betaini1 = beta.subcube(0, 1, 0, k-1, k*p + 1 - 1, 0);
  arma::cube betaini = beta.subcube(0, 1, 0, k-1, k*p + 1 - 1, g-1);

  arma::cube betaout = gamloopFista2(betaini, Y, Z, lambda, tol, YMean.t(), ZMean.t(), betaini1, k, p, tk);

  return(betaout);
}

// [[Rcpp::export]]
arma::cube HVARElemAlgcpp(const arma::cube &beta, const arma::mat& trainY, const arma::mat& trainZ, const arma::colvec& lambda,
                          const double& tol, const int& p, const int flag_restart_opt = 0){
  // Prelimaries
  int n = trainY.n_rows;
  int k = trainY.n_cols;
  int g = beta.n_slices;
  arma::mat YMean = mean(trainY);
  arma::mat ZMean = mean(trainZ.t());



  arma::mat Y = zeros(n, k);
  arma::mat trainZt = trainZ.t();
  arma::mat Z = zeros(n, k*p);
  for (int i = 0; i < n; i++) {
    Y.row(i) = trainY.row(i)  - YMean ;
    Z.row(i) = trainZt.row(i)  - ZMean ;
  }
  Z = Z.t();

  // Step size
  vec eigval;
  arma::mat eigvec;
  const arma::mat Zt=Z*trans(Z);
  eig_sym(eigval, eigvec, Zt);
  double tk=1/max(eigval);

  arma::mat betaini1 = beta.subcube(0, 1, 0, k-1, k*p + 1 - 1, 0);
  arma::cube betaini = beta.subcube(0, 1, 0, k-1, k*p + 1 - 1, g-1);
  arma::cube betaout = gamloopElem2(betaini, Y, Z, lambda, tol, YMean.t(), ZMean.t(), betaini1, k, p, tk, flag_restart_opt);

  return(betaout);
}

cv_aux_out HVAR_cvaux_cpp(const arma::mat & Y, const arma::mat& Z, int t, const arma::colvec& gamm, const double eps, int p, const double estim,
                          arma::cube Phi, const int flag_restart_opt = 0){

  int k = Y.n_cols;
  int g = gamm.size();
  arma::mat trainY = Y.submat(0, 0, t-1, k-1);
  arma::mat trainZ = Z.submat(0, 0, k*p-1, t-1);
  arma::mat testY = Y.submat(t, 0, t, k-1); // 1 row
  arma::mat testZ = Z.submat(0, t, k*p-1, t); // 1 col
  arma::mat testZc = ones(k*p+1, 1);
  testZc.submat(1, 0, k*p + 1 - 1, 0) = testZ;
  arma::mat MSFE = zeros(1, g);
  arma::mat sparsity = zeros(1, g);
  arma::mat pred;

  if(estim==1){
    Phi  = lassoVARFistcpp(Phi, trainY, trainZ, gamm, eps, p);
  }

  if(estim==2){
    Phi  = HVARElemAlgcpp(Phi, trainY, trainZ, gamm, eps, p, flag_restart_opt);
  }

  // Forecasts
  for (int igran = 0; igran < g; igran++) {
    arma::mat betagran = Phi.subcube(0, 0, igran, k-1,  k*p + 1 - 1, igran);
    arma::vec nz = nonzeros(betagran);
    pred = betagran*testZc;
    arma::mat diff = testY.t() - pred;
    MSFE.col(igran) = mean(pow(diff, 2));
    sparsity(igran) = nz.size() - k;
  }


  cv_aux_out my_cv_aux;
  my_cv_aux.MSFEs = MSFE;
  my_cv_aux.sparsity = sparsity;
  my_cv_aux.Phi = Phi;


  return(my_cv_aux);
}

cv_aux_out HVAR_cvaux_cpp(const arma::mat & Y, const arma::mat& Z, int t, const arma::colvec& gamm, const double eps, int p, const double estim){
  int k = Y.n_cols;
  int g = gamm.size();
  arma::cube Phi = zeros(k, k*p+1, g);
  return(HVAR_cvaux_cpp(Y, Z, t, gamm, eps, p, estim, Phi, 1));
}

// [[Rcpp::export]]
Rcpp::List HVAR_cvaux_loop_cpp(const arma::mat & Y, const arma::mat& Z, const arma::colvec& tseq, const arma::colvec& gamm, const double eps, int p, const double estim){

  // int k = Y.n_cols;
  int tlength = tseq.size();
  int g = gamm.size();
  Rcpp::List myList;
  arma::mat MSFEcv = zeros(tlength, g);
  arma::mat sparsitycv = zeros(tlength, g);
  cv_aux_out get_cv_out;

  for(int it = 0; it < tlength; it++){
    if (it == 0){
      get_cv_out = HVAR_cvaux_cpp(Y,  Z, tseq[it], gamm,  eps, p, estim);
    }
    else {
      get_cv_out = HVAR_cvaux_cpp(Y,  Z, tseq[it], gamm,  eps, p, estim, get_cv_out.Phi);
    }
    MSFEcv.submat(it, 0, it, g-1) = get_cv_out.MSFEs;
    sparsitycv.submat(it, 0, it, g-1) = get_cv_out.sparsity;
  }

  Rcpp::List results=Rcpp::List::create(
    Rcpp::Named("MSFEcv") = MSFEcv,
    Rcpp::Named("sparsitycv") = sparsitycv);
  return(results);
}


prox_out proxHVAR_cpp(const arma::mat& Ydata, const arma::mat& Zdata, const arma::mat& Xdata,
                      const double& lambdaPhi, const double& lambdaB,
                      const int& k, const int p, const int& kX, const int& s,
                      const double& eps, const double& max_iter, const double& alpha,
                      const arma::mat& Binit, const arma::mat& Phiinit){
  // Binit needs to be of dimension d TIMES kX*s
  // make sure that Ydata and Xdata are  matrices; be careful when kX =1 and s = 1

  // Number of time series
  int d = Ydata.n_cols;
  arma::mat Phi = zeros(k*p, d);
  arma::mat B = zeros(kX*s, d);

  // Compute step size
  arma::mat ZX = join_cols(Zdata, Xdata);
  vec svdval;
  mat svdU;
  mat svdV;
  arma::svd_econ(svdU, svdval, svdV, ZX);
  double tk = 1/pow(max(svdval), 2);

  arma::mat BinitT = arma::trans(Binit);
  arma::colvec B1 = BinitT.col(0);
  arma::mat PhiinitT = arma::trans(Phiinit);
  arma::colvec Phi1 = PhiinitT.col(0);

  for( int i =0; i < d; ++i){
    // Initializations
    B1 = BinitT.col(i);
    arma::colvec BOLD = B1;
    arma::colvec BOLDOLD = B1;

    Phi1 = PhiinitT.col(i);
    arma::colvec PhiOLD = Phi1;
    arma::colvec PhiOLDOLD = Phi1;

    // Convergence
    double threshPhi = eps*10;
    double threshB = eps*10;
    double thresh = eps*10;
    double it = 3;
    // arma::colvec B1;



    while((thresh > eps) & (it < max_iter) ){

      arma::colvec phi = PhiOLD + ((it-2)/(it+1))*(PhiOLD - PhiOLDOLD);
      arma::colvec beta = BOLD + ((it-2)/(it+1))*(BOLD - BOLDOLD);
      Phi1 = prox2HVARnew(arma::vectorise(phi) + tk*arma::vectorise((arma::trans(Ydata.col(i))-arma::trans(beta)*Xdata-arma::trans(phi)*Zdata)*arma::trans(Zdata)),lambdaPhi*tk, k, p);
      B1 = prox2HVARnew(arma::vectorise(beta) + tk*arma::vectorise((arma::trans(Ydata.col(i))-arma::trans(beta)*Xdata-arma::trans(phi)*Zdata)*arma::trans(Xdata)),lambdaB*tk, kX, s);

      Phi1 = (1/(1+alpha))*Phi1;
      B1 = (1/(1+alpha))*B1;

      threshPhi = max(abs(Phi1  - phi));
      threshB = max(abs(B1 - beta));
      thresh = std::max(threshB, threshPhi);
      PhiOLDOLD = PhiOLD;
      PhiOLD = Phi1;
      BOLDOLD = BOLD;
      BOLD = B1;

      it = it + 1;
    }

    B.col(i) = B1;
    Phi.col(i) = Phi1;
  }


  prox_out my_out;
  my_out.PHI = Phi;
  my_out.B = B;

  return(my_out);
}




prox_out proxBasic_cpp(const arma::mat& Ydata, const arma::mat& Zdata, const arma::mat& Xdata,
                       const double& lambdaPhi, const double& lambdaB,
                       const int& k, const int& p, const int& kX, const int& s,
                       const double& eps, const double& max_iter,
                       const double& alpha, const arma::mat& Binit, const arma::mat& Phiinit){

  // make sure that Ydata and Xdata are  matrices; be careful when kX =1 and s = 1

  int d = Ydata.n_cols;
  arma::mat Phi = zeros(k*p, d);
  arma::mat B = zeros(kX*s, d);

  // Compute step size
  arma::mat ZX = join_cols(Zdata, Xdata);
  vec svdval;
  mat svdU;
  mat svdV;
  arma::svd_econ(svdU, svdval, svdV, ZX);
  double tk = 1/pow(max(svdval), 2);

  arma::mat BinitT = arma::trans(Binit);
  arma::colvec B1 = BinitT.col(0);
  arma::mat PhiinitT = arma::trans(Phiinit);
  arma::colvec Phi1 = PhiinitT.col(0);

  for( int i =0; i < d; ++i){
    // Initializations
    B1 = BinitT.col(i);
    arma::colvec BOLD = B1;
    arma::colvec BOLDOLD = B1;

    Phi1 = PhiinitT.col(i);
    arma::colvec PhiOLD = Phi1;
    arma::colvec PhiOLDOLD = Phi1;

    // Convergence
    double threshPhi = eps*10;
    double threshB = eps*10;
    double thresh = eps*10;
    double it = 3;
    // arma::colvec B1;

    while((thresh > eps) & (it < max_iter) ){
      arma::colvec phi = PhiOLD + ((it-2)/(it+1))*(PhiOLD - PhiOLDOLD);
      arma::colvec beta = BOLD + ((it-2)/(it+1))*(BOLD - BOLDOLD);
      Phi1 = ST3a(arma::vectorise(phi) + tk*arma::vectorise((arma::trans(Ydata.col(i))-arma::trans(beta)*Xdata-arma::trans(phi)*Zdata)*arma::trans(Zdata)),lambdaPhi*tk);
      B1 = ST3a(arma::vectorise(beta) + tk*arma::vectorise((arma::trans(Ydata.col(i))-arma::trans(beta)*Xdata-arma::trans(phi)*Zdata)*arma::trans(Xdata)),lambdaB*tk);

      Phi1 = (1/(1+alpha))*Phi1;
      B1 = (1/(1+alpha))*B1;

      threshPhi = max(abs(Phi1  - phi));
      threshB = max(abs(B1 - beta));
      thresh = std::max(threshB, threshPhi);

      PhiOLDOLD = PhiOLD;
      PhiOLD = Phi1;
      BOLDOLD = BOLD;
      BOLD = B1;
      it = it + 1;
    }

    B.col(i) = B1;
    Phi.col(i) = Phi1;

  }

  prox_out my_out;
  my_out.PHI = Phi;
  my_out.B = B;

  return(my_out);
}


hvarx HVARX_NEW_cpp(const arma::mat& fullY, const arma::mat& fullZ, const arma::mat& fullX, const int& k, const int& kX, const int& p, const int& s,
                    const double& lambdaPhi, const double& lambdaB, const double& eps, const double& max_iter, const double& alpha, const int& type,
                    const arma::mat& Binit, const arma::mat& Phiinit){
  // type is an integer here 1: Lasso, 2: HLag

  // Inputs
  // fullY : matrix of responses of dimension nxd
  // fullX : matrix of predictors of dimension (kXs)xn
  // k : numeric, number of responses (columns in fullY)
  // kX : numeric, number of predictors (rows in fullX is kXs)
  // s : numeric, number of lagged values for the predictors
  // lambdaB : numeric, sparsity parameter
  // eps : numeric, convergence tolerance
  // max.iter : numeric, maximum number of iterations
  // alpha : numeric, small l2 penalty or not?
  // type : numeric, type of penalization c"HLag", "L1", "L2")

  // Preliminaries : de-mean the data
  int n = fullY.n_rows;
  arma::mat YMean;
  arma::mat Ydemean;
  arma::mat XMean;
  arma::mat Xdemean;
  arma::mat ZMean;
  arma::mat Zdemean;

  YMean = mean(fullY); // 1xk matrix as in R
  ZMean = mean(fullZ.t());
  XMean = mean(fullX.t());

  Ydemean = fullY - ones(n, 1)*YMean; // de-meaned Y
  Zdemean = fullZ.t() - ones(n, 1)*ZMean;
  Zdemean = Zdemean.t();
  Xdemean = fullX.t() - ones(n, 1)*XMean;
  Xdemean = Xdemean.t();

  prox_out out;
  // Proximal Gradient Algorithm
  if(type==1){ //L1
    out = proxBasic_cpp(Ydemean, Zdemean, Xdemean, lambdaPhi, lambdaB, k, p, kX, s, eps, max_iter, alpha, Binit, Phiinit);
  }

  if(type==2){ //HLag
    out = proxHVAR_cpp(Ydemean, Zdemean, Xdemean, lambdaPhi, lambdaB, k, p, kX, s, eps, max_iter, alpha, Binit, Phiinit);
  }


  arma::mat PHI = out.PHI;
  PHI = PHI.t();
  arma::mat B = out.B; // matrix of dimension (d)x(kX s)
  B = B.t();

  arma::mat phi0;
  arma::mat resids;



  phi0 = YMean - (B*XMean.t()).t() -  (PHI*ZMean.t()).t(); // matrix of dimension (d)x1


  resids = fullY - ones(n, 1)*phi0 -  fullX.t()*B.t() -  fullZ.t()*PHI.t();


  hvarx my_hvarx;
  my_hvarx.PHI = PHI;
  my_hvarx.B = B;
  my_hvarx.phi0 = phi0; // 1 x k

  return(my_hvarx);
}

// [[Rcpp::export]]
Rcpp::List HVARX_NEW_export_cpp(const arma::mat& fullY, const arma::mat& fullZ, const arma::mat& fullX, const int& k, const int& kX, const int& p, const int& s,
                    const double& lambdaPhi, const double& lambdaB, const double& eps, const double& max_iter, const double& alpha, const int& type,
                    const arma::mat& Binit, const arma::mat& Phiinit){
  // type is an integer here 1: Lasso, 2: HLag

  // Inputs
  // fullY : matrix of responses of dimension nxd
  // fullX : matrix of predictors of dimension (kXs)xn
  // k : numeric, number of responses (columns in fullY)
  // kX : numeric, number of predictors (rows in fullX is kXs)
  // s : numeric, number of lagged values for the predictors
  // lambdaB : numeric, sparsity parameter
  // eps : numeric, convergence tolerance
  // max.iter : numeric, maximum number of iterations
  // alpha : numeric, small l2 penalty or not?
  // type : numeric, type of penalization c"HLag", "L1", "L2")

  // Preliminaries : de-mean the data
  int n = fullY.n_rows;
  arma::mat YMean;
  arma::mat Ydemean;
  arma::mat XMean;
  arma::mat Xdemean;
  arma::mat ZMean;
  arma::mat Zdemean;

  YMean = mean(fullY); // 1xk matrix as in R
  ZMean = mean(fullZ.t());
  XMean = mean(fullX.t());

  Ydemean = fullY - ones(n, 1)*YMean; // de-meaned Y
  Zdemean = fullZ.t() - ones(n, 1)*ZMean;
  Zdemean = Zdemean.t();
  Xdemean = fullX.t() - ones(n, 1)*XMean;
  Xdemean = Xdemean.t();

  prox_out out;
  // Proximal Gradient Algorithm
  if(type==1){ //L1
    out = proxBasic_cpp(Ydemean, Zdemean, Xdemean, lambdaPhi, lambdaB, k, p, kX, s, eps, max_iter, alpha, Binit, Phiinit);
  }

  if(type==2){ //HLag
    out = proxHVAR_cpp(Ydemean, Zdemean, Xdemean, lambdaPhi, lambdaB, k, p, kX, s, eps, max_iter, alpha, Binit, Phiinit);
  }

  arma::mat PHI = out.PHI;
  PHI = PHI.t();
  arma::mat B = out.B; // matrix of dimension (d)x(kX s)
  B = B.t();

  arma::mat phi0;
  arma::mat resids;



  phi0 = YMean - (B*XMean.t()).t() -  (PHI*ZMean.t()).t(); // matrix of dimension (d)x1

  resids = fullY - ones(n, 1)*phi0 -  fullX.t()*B.t() -  fullZ.t()*PHI.t();


  List results = List::create(Named("fullY") = fullY, Named("fullZ") = fullZ, Named("fullX") = fullX,
                                    Named("Phi") = PHI, Named("B") = B, Named("phi0") = phi0, Named("resids") = resids);

  return(results);
}

cv_aux_out HVARX_cvaux_cpp(const arma::mat& Y, const arma::mat& Z, const arma::mat& X, int t,
                           const arma::colvec& lambdaPhiseq,  const arma::colvec& lambdaBseq,
                           const double eps, const double max_iter, int k, int kX, int p, int s, 
                           const double alpha, const double estim, 
                           arma::cube myBinit, arma::cube myPhiinit, const int flag_restart_opt = 0){


  int g = lambdaPhiseq.size(); // 10*10
  arma::mat trainY = Y.submat(0, 0, t-1, k-1);
  arma::mat trainX = X.submat(0, 0, kX*s-1, t-1);
  arma::mat trainZ = Z.submat(0, 0, k*p-1, t-1);
  arma::mat testY = Y.submat(t, 0, t, k-1); // 1 row
  arma::mat testX = X.submat(0, t, kX*s-1, t); // 1 col
  arma::mat testZ = Z.submat(0, t, k*p-1, t); // 1 col

  hvarx hvarx_out;
  arma::mat Phifit;
  arma::cube Phicv = zeros(k, k*p, g);
  arma::mat Bfit;
  arma::cube Bcv = zeros(k, kX*s, g);
  arma::mat phi0fit;

  arma::mat MSFE = zeros(1, g);
  arma::mat sparsity = zeros(1, g);
  arma::mat pred;


  Bfit = myBinit.slice(0);
  Phifit = myPhiinit.slice(0);
  for (int igran = 0; igran < g; igran++) {

    if (flag_restart_opt == 0){
      Bfit = myBinit.slice(igran);
      Phifit = myPhiinit.slice(igran);
    }

    hvarx_out = HVARX_NEW_cpp(trainY, trainZ, trainX, k, kX, p, s, lambdaPhiseq[igran], lambdaBseq[igran], eps, max_iter, alpha,  estim,  Bfit, Phifit);


    Phifit = hvarx_out.PHI;
    Phicv.slice(igran) = Phifit;
    Bfit = hvarx_out.B;
    Bcv.slice(igran) = Bfit;
    phi0fit = hvarx_out.phi0;

    pred = Bfit*testX + Phifit*testZ + phi0fit.t();

    arma::mat diff = testY.t() - pred;
    MSFE.col(igran) = mean(pow(diff, 2));
    arma::vec Phinz = nonzeros(Phifit);
    arma::vec Bnz = nonzeros(Bfit);
    sparsity(igran) = Phinz.size() + Bnz.size();
  }

  cv_aux_out hvarx_cv_out;
  hvarx_cv_out.MSFEs = MSFE;
  hvarx_cv_out.sparsity = sparsity;
  hvarx_cv_out.Phi = Phicv;
  hvarx_cv_out.varxB = Bcv;
  return(hvarx_cv_out);
}


// Overloading function for the cases in which no starting Binit and Phiinit are given.
// in these cases we just set everything to zero.
cv_aux_out HVARX_cvaux_cpp(const arma::mat& Y, const arma::mat& Z, const arma::mat& X, int t,
                           const arma::colvec& lambdaPhiseq,  const arma::colvec& lambdaBseq,
                           const double eps, const double max_iter, int k, int kX, int p, int s, 
                           const double alpha, const double estim){
  int g = lambdaPhiseq.size(); // 10*10
  arma::cube myBinit = zeros(k, kX*s, g);
  arma::cube myPhiinit = zeros(k, k*p, g);
  return(HVARX_cvaux_cpp(Y, Z, X, t, lambdaPhiseq, lambdaBseq, eps, max_iter, k, kX, p, s, alpha, estim, myBinit, myPhiinit, 1));
}

// [[Rcpp::export]]
Rcpp::List HVARX_cvaux_cpp_loop(const arma::mat& Y, const arma::mat& Z, const arma::mat& X,
                                const arma::colvec& lambdaPhiseq,  const arma::colvec& lambdaBseq,
                                const double eps, const double max_iter, int k, int kX, int p, int s, const double alpha, const double estim,
                                const arma::colvec& tseq){

  // int k = Y.n_cols;
  int tlength = tseq.size();
  int g = lambdaPhiseq.size();
  arma::mat MSFEcv = zeros(tlength, g);
  arma::mat sparsitycv = zeros(tlength, g);
  cv_aux_out get_cv_out;

  for(int it = 0; it < tlength; it++){
    if (it == 0){
      get_cv_out = HVARX_cvaux_cpp(Y, Z, X, tseq[it], lambdaPhiseq, lambdaBseq, eps, max_iter, k,  kX,  p,  s,  alpha,  estim);
    }
    else {
      // For some reasons it would give me different results if I would also take the previous B as starting values ...
      // TODO: find out why!
      get_cv_out = HVARX_cvaux_cpp(Y, Z, X, tseq[it], lambdaPhiseq, lambdaBseq, eps, max_iter, k,  kX,  p,  s,  alpha,  estim, get_cv_out.varxB, get_cv_out.Phi);
    }
    MSFEcv.submat(it, 0, it, g-1) = get_cv_out.MSFEs;
    sparsitycv.submat(it, 0, it, g-1) = get_cv_out.sparsity;
  }

  Rcpp::List results=Rcpp::List::create(
    Rcpp::Named("MSFEcv") = MSFEcv,
    Rcpp::Named("sparsitycv") = sparsitycv);
  return(results);
}



bool moveup_LGSearch_cpp(arma::mat &param){
  int n = param.n_rows;
  int k = param.n_cols;

  int i = 0;
  int j = 0;
  for (int i=0; i<n; i++){
    for (int j=0; j<k; j++){
      if (param(i, j) != 0) return(false);
    }
  }
  return(true);
}

// estim 1 = lasso
// estim 2 = hvar
// [[Rcpp::export]]
double LGSearch_cpp(double gstart, arma::mat &Y, arma::mat &Z, arma::cube beta, int estim, int k, int p) {
  double lambdah = gstart;
  double lambdal = 0.0;
  arma::mat param;

  while (std::abs(lambdah - lambdal) > 0.00001){
    double lambda = (lambdah + lambdal)/2;
    arma::colvec lvec(1);
    lvec.fill(lambda);
    if (estim == 1){
      beta = lassoVARFistcpp(beta, Y, Z, lvec, 0.0001, p);
      param = beta.slice(0).cols(1, k*p);
    }
    else if (estim == 2){
      beta = HVARElemAlgcpp(beta, Y, Z, lvec, 0.0001, p);
      param = beta.slice(0).cols(1, k*p);
    }

    bool move_up = moveup_LGSearch_cpp(param);
    if (move_up) lambdah = lambda;
    else lambdal = lambda;

  }

  return(lambdah);
}