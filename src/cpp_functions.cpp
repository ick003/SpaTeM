#include <Rcpp.h>
using namespace Rcpp;


//' Function to suppress rows from a matrix
//' @param x matrix which rows need to be erased
//' @param rowID row index
//' @examples
//' x<- matrix(1:9, ncol = 3)
//' row_erase(x,1:2)
//' @export
// [[Rcpp::export]]
NumericMatrix row_erase (NumericMatrix x, IntegerVector rowID) {
  rowID = rowID.sort();

  NumericMatrix x2(Dimension(x.nrow()- rowID.size(), x.ncol()));
  int iter = 0;
  int del = 1; // to count deleted elements
  for (int i = 0; i < x.nrow(); i++) {
    if (i != rowID[del - 1]) {
      x2.row(iter) = x.row(i);
      iter++;
    } else {
      del++;
    }
  }
  return x2;
}

// declare row_erase
NumericMatrix row_erase(NumericMatrix x, IntegerVector rowID);

//' Function to get the bonds based on neighbourhood
//' @param Location matrix of coordinates.
//' @param NN maximum number of neighbours.
//' @examples
//' Obs.loc = as.matrix(expand.grid(1:4,1:4))
//' getBonds(Obs.loc, NN = 4)
//' @export
// [[Rcpp::export]]
NumericMatrix getBonds(NumericMatrix Location, double NN, double th = 2) {

  int nvert = Location.nrow();

  NumericMatrix bond(NN*nvert,2);


  //  NumericMatrix res = Location;
  int k = 0;
  for (int i = 0; i < nvert; i++){

    for( int j = i+1; j < nvert; j++){
      double xj = (Location(j,0) - Location(i,0))*(Location(j,0) - Location(i,0));
      double yj = (Location(j,1) - Location(i,1))*(Location(j,1) - Location(i,1));
      double dist = xj + yj;

      if(dist < th){
        bond(k,0) = i+1;
        bond(k,1) = j+1;
        k++;
      }
    }
  }
  IntegerVector idx(NN*nvert - k);
  for (int i = 0; i < NN*nvert-k; i++){
    idx(i) = i + k;
  }

  bond = row_erase(bond, idx);

  return bond;
}


//' Function to select bonds based on colors and parameters
//' @param Bds matrix of neighbouring bonds.
//' @param Cols vector of colors per vertex.
//' @param Betas vector of parameters for the Potts model.
//' @examples
//' Obs.loc = expand.grid(1:10,1:10)
//' Bds = getBonds(Obs.loc)
//' Bds = Bds[which(Bds[,1]>0),]
//' Betas = c(0.8, 0.8, 0.2)
//' Cols = sample(1:3, 100, replace = TRUE)
//' Bd = selectBonds(Bds, Cols, Betas)
//' par(mfrow = c(1,1))
//' col = grey.colors(3)
//' image(matrix(Cols, ncol = 10), col = col)
//' segments((Obs.loc[Bd[,1],1]-1)/9, (Obs.loc[Bd[,1],2]-1)/9, (Obs.loc[Bd[,2],1]-1)/9, (Obs.loc[Bd[,2],2]-1)/9)
//' segments((Obs.loc[,1]-1)/9, (Obs.loc[,2]-1)/9, (Obs.loc[,1]-1)/9, (Obs.loc[,2]-1)/9, lwd = 4)
//' @export
// [[Rcpp::export]]
NumericMatrix selectBonds(NumericMatrix Bds, NumericVector Cols, NumericVector Betas) {

  int nbds = Bds.nrow();

  NumericMatrix Bd(nbds,2);

  NumericMatrix dC(nbds,2);

// Same color bonds
int k = 0;
  for (int i = 0; i < nbds; i++){
    double diffCol = Cols(Bds(i,0)-1) - Cols(Bds(i,1)-1);
    if(diffCol == 0){
      double pr = -expm1(-Betas(Cols(Bds(i,1)-1)-1));
      double u = unif_rand();
      dC(i,0) = u;
      dC(i,1) = pr;
      if(u < pr){
        Bd(k,0) = Bds(i,0);
        Bd(k,1) = Bds(i,1);
        k++;
      }
    }
  }

  IntegerVector idx(nbds - k);
  for (int i = 0; i < nbds-k; i++){
    idx(i) = i + k;
  }

  Bd = row_erase(Bd, idx);

  return Bd;
}


//' Function to identify components based on bonds.
//' @param Bds matrix of bonds.
//' @param nvert number of vertices.
//' @examples
//' Obs.loc = expand.grid(1:10,1:10)
//' Bds = getBonds(Obs.loc)
//' Bds = Bds[which(Bds[,1]>0),]
//' Betas = c(0.8, 0.8, 0.2)
//' Cols = sample(1:3, 100, replace = TRUE)
//' Bd = selectBonds(Bds, Cols, Betas)
//' Bd = Bd[which(Bd[,1]>0),]
//' par(mfrow = c(1,1))
//' col = grey.colors(3)
//' image(matrix(Cols, ncol = 10), col = col)
//' segments((Obs.loc[Bd[,1],1]-1)/9, (Obs.loc[Bd[,1],2]-1)/9, (Obs.loc[Bd[,2],1]-1)/9, (Obs.loc[Bd[,2],2]-1)/9)
//' segments((Obs.loc[,1]-1)/9, (Obs.loc[,2]-1)/9, (Obs.loc[,1]-1)/9, (Obs.loc[,2]-1)/9, lwd = 33, col = "black")
//' segments((Obs.loc[,1]-1)/9, (Obs.loc[,2]-1)/9, (Obs.loc[,1]-1)/9, (Obs.loc[,2]-1)/9, lwd = 32, col = "white")
//' Comp = setComp(Bd, 100)
//' col.comp = rainbow(100)
//' text((Obs.loc[,1]-1)/div,(Obs.loc[,2]-1)/div, labels = 1:100, col = col.comp[Comp])
//' @export
// [[Rcpp::export]]
NumericVector setComp(NumericMatrix Bds, double nvert) {

  int nbds = Bds.nrow();

  NumericVector Comp(nvert);
  NumericVector Comp_f(nvert);
  // Same color bonds

  for (int i = 0; i < nvert; i++){
    Comp(i) = i+1;
  }
double diffComp = 1;

  while(diffComp > 0){
    for (int i = 0; i < nvert; i++){
      for (int j = 0; j < nbds; j++){
        double bdx = Bds(j,0)-1;
        double bdy = Bds(j,1)-1;
          if(bdx == i){
            Comp(bdy) = min(NumericVector::create(Comp(i), Comp(bdy)));
            Comp(i) = min(NumericVector::create(Comp(i), Comp(bdy)));
            }
          if(bdy == i){
            Comp(bdx) = min(NumericVector::create(Comp(i), Comp(bdx)));
            Comp(i) = min(NumericVector::create(Comp(i), Comp(bdx)));
          }
        }
    }
    diffComp = 0;
    for (int i = 0; i < nvert; i++){
      diffComp = diffComp + (Comp(i)- Comp_f(i)) * (Comp(i)- Comp_f(i));
    }
    Comp_f = Comp;
  }
  return Comp;
}


//' Function to summarize the components' statistics (size and colors).
//' @param Bds matrix of bonds.
//' @param nvert number of vertices.
//' @examples
//' Obs.loc = expand.grid(1:10,1:10)
//' Bds = getBonds(Obs.loc)
//' Bds = Bds[which(Bds[,1]>0),]
//' Betas = c(0.8, 0.8, 0.2)
//' Cols = sample(1:3, 100, replace = TRUE)
//' Bd = selectBonds(Bds, Cols, Betas)
//' Bd = Bd[which(Bd[,1]>0),]
//' par(mfrow = c(1,1))
//' col = grey.colors(3)
//' image(matrix(Cols, ncol = 10), col = col)
//' segments((Obs.loc[Bd[,1],1]-1)/9, (Obs.loc[Bd[,1],2]-1)/9, (Obs.loc[Bd[,2],1]-1)/9, (Obs.loc[Bd[,2],2]-1)/9)
//' segments((Obs.loc[,1]-1)/9, (Obs.loc[,2]-1)/9, (Obs.loc[,1]-1)/9, (Obs.loc[,2]-1)/9, lwd = 33, col = "black")
//' segments((Obs.loc[,1]-1)/9, (Obs.loc[,2]-1)/9, (Obs.loc[,1]-1)/9, (Obs.loc[,2]-1)/9, lwd = 32, col = "white")
//' Comp = setComp(Bd, 100)
//' col.comp = rainbow(100)
//' text((Obs.loc[,1]-1)/div,(Obs.loc[,2]-1)/div, labels = 1:100, col = col.comp[Comp])
//' @export
// [[Rcpp::export]]
NumericMatrix getCompStat(NumericVector Comp, NumericVector CompID, NumericVector Cols) {

  int nComp = CompID.length();
  int nVert = Cols.length();

  NumericMatrix compDF(nComp, 3);

  for (int i = 0; i < nComp; i++){
    compDF(i,0) = CompID(i);
    for (int j = 0; j < nVert; j++){
      if(Comp(j) == CompID(i)){
        compDF(i,1)++ ;
        compDF(i,2) = Cols(j);
      }
    }
  }
  return compDF;
}

//' Function to flip components spins.
//' @param CompID vector of components IDs.
//' @param CompIDList Vector of possible components ID.
//' @param ncolor number of spins.
//' @examples
//' Obs.loc = expand.grid(1:10,1:10)
//' Bds = getBonds(Obs.loc)
//' Bds = Bds[which(Bds[,1]>0),]
//' Betas = c(0.8, 0.8, 0.2)
//' Cols = sample(1:3, 100, replace = TRUE)
//' Bd = selectBonds(Bds, Cols, Betas)
//' Bd = Bd[which(Bd[,1]>0),]
//' par(mfrow = c(1,1))
//' col = grey.colors(3)
//' image(matrix(Cols, ncol = 10), col = col)
//' segments((Obs.loc[Bd[,1],1]-1)/9, (Obs.loc[Bd[,1],2]-1)/9, (Obs.loc[Bd[,2],1]-1)/9, (Obs.loc[Bd[,2],2]-1)/9)
//' segments((Obs.loc[,1]-1)/9, (Obs.loc[,2]-1)/9, (Obs.loc[,1]-1)/9, (Obs.loc[,2]-1)/9, lwd = 33, col = "black")
//' segments((Obs.loc[,1]-1)/9, (Obs.loc[,2]-1)/9, (Obs.loc[,1]-1)/9, (Obs.loc[,2]-1)/9, lwd = 32, col = "white")
//' Comp = setComp(Bd, 100)
//' col.comp = rainbow(100)
//' text((Obs.loc[,1]-1)/div,(Obs.loc[,2]-1)/div, labels = 1:100, col = col.comp[Comp])
//' @export
// [[Rcpp::export]]
NumericVector flipComp(NumericVector CompID, NumericVector CompIDList, double ncolor) {

  int nvert = CompID.length();
  int nComp = CompIDList.length();

  NumericVector Spin(nvert);

  for(int i = 0; i < nComp; i++){
    double id = CompIDList(i);
    double spin = unif_rand() * ncolor;
    for(int j = 0; j < nvert; j++){
      double idComp = CompID(j);
      if(idComp == id){
        Spin(j) = spin;
      }
    }
  }
  return Spin;
}

//' Function to calculate the canonical statistic.
//' @param Bds matrix of bonds
//' @param Cols Vector of colors.
//' @param ncolor number of spins.
//' @examples
//' Obs.loc = as.matrix(expand.grid(1:4,1:4))
//' Bds = getBonds(Obs.loc, NN = 4)
//' Bds = Bds[which(Bds[,1]>0),]
//' Betas = c(0.8, 0.8, 0.2)
//' Cols = sample(1:3, 16, replace = TRUE)
//' par(mfrow = c(1,1))
//' col = grey.colors(3)
//' image(matrix(Cols, ncol = 4), col = col)
//' CS = getCanStatF(Bds, Cols, 3)
//' @export
// [[Rcpp::export]]
NumericMatrix getCanStatF(NumericMatrix Bds, NumericVector Cols, double ncolors) {

  int nbds = Bds.nrow();
  int nvert = Cols.length();

  NumericMatrix CanStatF(nvert,ncolors);

  for(int i = 0; i < nbds; i++){
    if(Cols(Bds(i,0)-1) == Cols(Bds(i,1)-1)){
        CanStatF(Bds(i,0)-1,Cols(Bds(i,0)-1)-1)++;
        CanStatF(Bds(i,1)-1,Cols(Bds(i,1)-1)-1)++;
    }
  }
  return CanStatF;
}

//' Function to calculate the canonical statistic.
//' @param Bds matrix of bonds
//' @param Cols Vector of colors.
//' @param ncolor number of spins.
//' @examples
//' Obs.loc = as.matrix(expand.grid(1:4,1:4))
//' Bds = getBonds(Obs.loc, NN = 4)
//' Bds = Bds[which(Bds[,1]>0),]
//' Betas = c(0.8, 0.8, 0.2)
//' Cols = sample(1:3, 16, replace = TRUE)
//' par(mfrow = c(1,1))
//' col = grey.colors(3)
//' image(matrix(Cols, ncol = 4), col = col)
//' CS = getCanStat(Bds, Cols, 3)
//' @export
// [[Rcpp::export]]
NumericVector getCanStat(NumericMatrix Bds, NumericVector Cols, double ncolors) {

  int nbds = Bds.nrow();
  int nvert = Cols.length();

  NumericVector CanStat(ncolors);

  for(int i = 0; i < nbds; i++){
    if(Cols(Bds(i,0)-1) == Cols(Bds(i,1)-1)){
      CanStat(Cols(Bds(i,0)-1)-1)++;
    }
  }
  return CanStat;
}

//' Function performing the loop in the SW algorithm.
//' @param Bds matrix of bonds
//' @param Cols Vector of colors.
//' @param ncolor number of spins.
//' @param Nrun number of runs.
//' @param Betas vector of parameters.
//' @examples
//' Obs.loc = as.matrix(expand.grid(1:10,1:10))
//' Bds = getBonds(Obs.loc, NN = 4)
//' Betas = c(0.8, 0.8, 0.2)
//' Cols = sample(1:3, 100, replace = TRUE)
//' par(mfrow = c(1,1))
//' col = grey.colors(3)
//' image(matrix(Cols, ncol = 10), col = col)
//' CS = loopSW(Bds, Cols, 3, 1000, Betas)
//' image(matrix(CS, ncol = 10), col = col)
//' @export
// [[Rcpp::export]]
NumericVector loopSW(NumericMatrix Bds, NumericVector Cols, double ncolors, int Nrun, NumericVector Betas) {

  int nbds = Bds.nrow();
  int nvert = Cols.length();

  NumericVector ColsNew = Cols;

  for(int i = 0; i < Nrun; i++){

    // Select Bonds
    NumericMatrix Comp_Bonds = selectBonds(Bds, ColsNew, Betas);

    //Set Components
    NumericVector Comp = setComp(Comp_Bonds, nvert);



    // Flipping spins for each component
    NumericVector z_temp = flipComp(Comp, unique(Comp), ncolors);

    ColsNew = ceiling(z_temp);

  }

  return ColsNew;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

// [[Rcpp::export]]
NumericMatrix RcppGibbs(int n, int thn) {
  
  int i,j;
  NumericMatrix mat(n, 2);
  
  // The rest of the code follows the R version
  double x=0, y=0;
  
  for (i=0; i<n; i++) {
    for (j=0; j<thn; j++) {
      x = R::rgamma(3.0,1.0/(y*y+4));
      y = R::rnorm(1.0/(x+1),1.0/sqrt(2*x+2));
    }
    mat(i,0) = x;
    mat(i,1) = y;
  }
  
  return mat;             // Return to R
}

// [[Rcpp::export]]
NumericMatrix RcppLoopGibbs(int n, int thn, 
                            NumericMatrix alphaH, NumericMatrix betaH,
                            NumericVector gammaH, NumericVector sigma2H,
                            NumericVector tauH, NumericVector rhoH,
                            NumericVector phiH, NumericMatrix piH, NumericMatrix zH,
                            NumericMatrix Xspt, NumericMatrix Xbasis) {
  
  int i,j;
  int nAlpha = alphaH.nrow() * alphaH.ncol();
  int nBeta = betaH.nrow() * betaH.ncol();
  int nGamma = gammaH.length();
  int nSigma = sigma2H.length();
  int nTau = tauH.length();
  int nRho = rhoH.length();
  int nPhi = phiH.length();
  int nPi =piH.nrow() * piH.ncol();
  int nZ = zH.nrow() * zH.ncol();
  int nComp = zH.ncol();
  int nSites = zH.nrow();
  int nBasis = betaH.ncol();
  
  NumericMatrix mat(n, nAlpha + nBeta + nGamma + nSigma + nTau + nRho + nPhi + nPi + nZ);
  
  // The rest of the code follows the R version
  double x=0, y=0;
  
  for (i=0; i<n; i++) {
    for (j=0; j<thn; j++) {
      x = R::rgamma(3.0,1.0/(y*y+4));
      y = R::rnorm(1.0/(x+1),1.0/sqrt(2*x+2));
    }
    mat(i,0) = x;
    mat(i,1) = y;
  }
  
  return mat;             // Return to R
}

/*** R
Obs.loc = expand.grid(1:10,1:10)
bond = getBonds(as.matrix(Obs.loc), NN = 8, th = 3)

Betas = c(1.8, 0.2, 0.2)
Cols = sample(1:3, 100, replace = TRUE)

Bd = selectBonds(bond, Cols, Betas)

par(mfrow = c(2,2))
col = grey.colors(3)
image(matrix(Cols, ncol = 10), col = col)
div = 10-1
segments((Obs.loc[Bd[,1],1]-1)/div, (Obs.loc[Bd[,1],2]-1)/div, (Obs.loc[Bd[,2],1]-1)/div, (Obs.loc[Bd[,2],2]-1)/div)
segments((Obs.loc[,1]-1)/div, (Obs.loc[,2]-1)/div, (Obs.loc[,1]-1)/div, (Obs.loc[,2]-1)/div, lwd = 13, col = "black")
segments((Obs.loc[,1]-1)/div, (Obs.loc[,2]-1)/div, (Obs.loc[,1]-1)/div, (Obs.loc[,2]-1)/div, lwd = 12, col = "white")

Comp = setComp(matrix(Bd, ncol = 2), nvert = 100)

getCompStat(Comp, unique(Comp), Cols)

col.comp = rainbow(100)
text((Obs.loc[,1]-1)/div,(Obs.loc[,2]-1)/div, labels = 1:25, col = col.comp[Comp])

Spin = ceiling(flipComp(Comp, unique(Comp), ncolor = 3))

Cols.new = ceiling(Spin)

image(matrix(Cols.new, ncol = 10), col = col)
segments((Obs.loc[Bd[,1],1]-1)/div, (Obs.loc[Bd[,1],2]-1)/div, (Obs.loc[Bd[,2],1]-1)/div, (Obs.loc[Bd[,2],2]-1)/div)
segments((Obs.loc[,1]-1)/div, (Obs.loc[,2]-1)/div, (Obs.loc[,1]-1)/div, (Obs.loc[,2]-1)/div, lwd = 13, col = "black")
segments((Obs.loc[,1]-1)/div, (Obs.loc[,2]-1)/div, (Obs.loc[,1]-1)/div, (Obs.loc[,2]-1)/div, lwd = 12, col = "white")
text((Obs.loc[,1]-1)/div,(Obs.loc[,2]-1)/div, labels = 1:25, col = col.comp[Comp])

Bd = selectBonds(bond, Cols.new, Betas)

image(matrix(Cols.new, ncol = 10), col = col)
div = 10-1
segments((Obs.loc[Bd[,1],1]-1)/div, (Obs.loc[Bd[,1],2]-1)/div, (Obs.loc[Bd[,2],1]-1)/div, (Obs.loc[Bd[,2],2]-1)/div)
segments((Obs.loc[,1]-1)/div, (Obs.loc[,2]-1)/div, (Obs.loc[,1]-1)/div, (Obs.loc[,2]-1)/div, lwd = 13, col = "black")
segments((Obs.loc[,1]-1)/div, (Obs.loc[,2]-1)/div, (Obs.loc[,1]-1)/div, (Obs.loc[,2]-1)/div, lwd = 12, col = "white")

Comp = setComp(Bd, 100)

col.comp = rainbow(100)
text((Obs.loc[,1]-1)/div,(Obs.loc[,2]-1)/div, labels = 1:100, col = col.comp[Comp])


cbind(1:100,Comp)

*/
