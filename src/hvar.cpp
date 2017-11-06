#include <RcppArmadillo.h>
#include <vector>
#include <limits>
#include <algorithm>
#include <numeric>      // std::iota

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

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
