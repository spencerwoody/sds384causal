
#include <Rcpp.h>

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// [[Rcpp::export]]
Rcpp::NumericVector randomShuffle(Rcpp::NumericVector a) {

    // clone a into b to leave a alone
    Rcpp::NumericVector b = Rcpp::clone(a);

    std::random_shuffle(b.begin(), b.end(), randWrapper);

    return b;
}

IntegerVector RankFun (NumericVector Y) {

  double Yi;
  int RankI;

  int N = Y.size();

  IntegerVector RankVec(N);

  for (int i = 0; i < N; i++) {

    Yi = Y[i];

    RankI = 0;

    for (int j = 0; j < N; j++) {
      if (Y[j] < Y[i]) {
	RankI++;
      }
    }

    RankVec[i] = RankI;

  }
  
  return RankVec;
  
}

IntegerVector TieFun (NumericVector Y) {

  double Yi;
  int RankI;

  int N = Y.size();

  IntegerVector RankVec(N);

  for (int i = 0; i < N; i++) {

    Yi = Y[i];

    RankI = 0;

    for (int j = 0; j < N; j++) {
      if (Y[j] == Y[i]) {
	RankI++;
      }
    }

    RankVec[i] = RankI;

  }
  
  return RankVec;
  
}
