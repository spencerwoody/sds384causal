
RankFun <- Rcpp::cppFunction("IntegerVector RankFun (NumericVector Y) {

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
  
}")

TieFun <- Rcpp::cppFunction("IntegerVector TieFun (NumericVector Y) {

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
  
}")

normRank <- function(W, Y, C = 0) {

  N <- length(W)

  W <- as.numeric(W)

  if (all(sort(unique(W)) != c(0, 1))) stop("W must be a binary treatment")

  Y[W = 1] <- Y[W = 1] - C

  Rank <- rep(NA, N)

  Rank <- RankFun(Y) + (1 + TieFun(Y)) / 2 - (N + 1) / 2

  RTbar <- mean(Rank[W == 1])
  RCbar <- mean(Rank[W == 0])

  Trank <- abs(RTbar - RCbar) 

  Trank

}

