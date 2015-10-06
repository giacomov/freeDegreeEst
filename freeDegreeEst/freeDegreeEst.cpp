/*
  y = freeDegreeTIEst(x,estLength,pen,maxPolyOrder);
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/stl_iterator.hpp>


//#include "mex.h"
//#include "matrix.h"
//#include "mathHelper.c"

/* Code added by GV for compiling outside matlab */

double realmax = 1e238 ;
double realmin = 1e-307 ;
double PI = 3.14159265359;

//inline int max ( int a, int b ) { return a > b ? a : b; }
//inline int min ( int a, int b ) { return a < b ? a : b; }

using namespace std;

/* End of code added by GV */

struct tree
{
  struct tree *leftChild;
  struct tree *rightChild;
  int nObs, firstObsInInterval, maxPolyOrder, optPolyOrder;
  double **polyCoeffs, *likelihoods;
  double intervalLeft, intervalRight,optLikelihood;
};

void printTreeNode(struct tree *T, FILE *outfile)
{
  int i,k;
  fprintf(outfile, "Interval [%f,%f], polynomial order = %d\n",
	 T->intervalLeft,T->intervalRight, T->optPolyOrder);
  fprintf(outfile, "\tChebyshev polynomial coefficients:\t");
  for (k = 0; k < T->optPolyOrder; k++) {
    fprintf(outfile, "%lf\t",T->polyCoeffs[T->optPolyOrder-1][k]);
  }
  fprintf(outfile,"\n");
}

void printTree(struct tree *T, FILE *outFile)
{
  if (T == NULL) 
    return;
  if (T->optPolyOrder > -1)
    printTreeNode(T,outFile);
  else {
    printTree(T->leftChild,outFile);
    printTree(T->rightChild,outFile);
  }
}

void printTreeNodeDebug(struct tree *T, int level, double *x)
{
  int i,k;
  printf("\nLevel = %d, nObs = %d, left interval = %f, right interval = %f\n",
	 level,T->nObs,
	 T->intervalLeft,T->intervalRight);
  printf("\tmaxPolyOrder = %d, optPolyOrder = %d\n",
	 T->maxPolyOrder,T->optPolyOrder);
  for (i = 0; i < T->maxPolyOrder; i++) {
    printf("%d order poly coeffs:\t",i+1);
    for (k = 0; k <= i; k++) {
      printf("\t%lf",T->polyCoeffs[i][k]);
    }
    printf("\n");
  }
  printf("Likelihoods\n");
  for (i = 0; i < T->maxPolyOrder; i++) {
    printf("\t%lf",T->likelihoods[i]);
  }
  printf("\n");
  /*
  for (i = 0; i < T->nObs; i++)
  {
    printf("\t%lf",x[i+T->firstObsInInterval]);
  }
  if (i > 0)
    printf("\n");
  */
}
void printTreeDebug(struct tree *T, int level, double *x)
{
  if (T == NULL) 
    return;
  printTreeNodeDebug(T,level,x);
  printTreeDebug(T->leftChild,level+1,x);
  printTreeDebug(T->rightChild,level+1,x);
}

void freeTree(struct tree *T)
{
	int i;
	/*printf("fT leftChild\n");*/
  if (T->leftChild != NULL)   {
    freeTree(T->leftChild);
  }
	/*printf("fT rightChild\n");*/
  if (T->rightChild != NULL)   {
    freeTree(T->rightChild);
  }
	/*printf("fT polyCoeffs\n");*/
	if (T->polyCoeffs != NULL) {
		/*printf("fT maxPolyOrder = %d\n",T->maxPolyOrder);*/
		for (i = 0; i < T->maxPolyOrder; i++) {
			free(T->polyCoeffs[i]);
		}
		free(T->polyCoeffs);
	} 
	/*printf("fT likelihoods\n");*/
	if (T->likelihoods != NULL) {
		free(T->likelihoods);
	}
	/*printf("fT T\n");*/
  free(T);
}


void coeffsToEval2(double *x, struct tree *T, int polyOrder, double *est) {
  double *y, *da, *db, *dc;
  int i,k;
  int estLength = T->nObs;
  double a,b;
  
  a = T->intervalLeft;
  b = T->intervalRight;
  y = (double *)malloc(sizeof(double)*estLength);
  da = (double *)malloc(sizeof(double)*estLength);
  db = (double *)malloc(sizeof(double)*estLength);
  dc = (double *)malloc(sizeof(double)*estLength);
  for (i = 0; i < T->nObs; i++) {
    y[i] = (x[i+T->firstObsInInterval]-a)/(b-a)*2.0-1.0;
    da[i] = 0.0;
    db[i] = 0.0;
    dc[i] = 0.0;
  }
  for (k = polyOrder-1; k>=1; k--) {
    for (i = 0; i < T->nObs; i++) {
      dc[i] = 2*y[i]*db[i]-da[i]+T->polyCoeffs[polyOrder-1][k];
      da[i] = db[i];
      db[i] = dc[i];
    }
  }
  for (i = 0; i < T->nObs; i++) {
    est[i] = y[i]*dc[i]-da[i]+T->polyCoeffs[polyOrder-1][0];
  }
  free(y);
  free(da);
  free(db);
  free(dc);
}

void coeffsToEval(double *x, struct tree *T, int polyOrder, double *est) {
  double *cheb, *chebPrev, *chebNext, *y;
  int i,k;
  int estLength = T->nObs;
  double a,b;
  
  a = T->intervalLeft;
  b = T->intervalRight;

  y = (double *)malloc(sizeof(double)*estLength);
  cheb = (double *)malloc(sizeof(double)*estLength);
  chebNext = (double *)malloc(sizeof(double)*estLength);
  chebPrev = (double *)malloc(sizeof(double)*estLength);
  for (i = 0; i < T->nObs; i++) {
    y[i] = (x[i+T->firstObsInInterval]-a)/(b-a)*2.0-1.0;
    cheb[i] = 1.0;
    est[i] = 0.0;
  }
  for (k = 0; k < polyOrder; k++) {
    for (i = 0; i < estLength; i++) {
      est[i] += cheb[i]*T->polyCoeffs[polyOrder-1][k];
    }
    if (k==0) {
      for (i = 0; i < estLength; i++) {
	chebNext[i] = y[i];
      }
    }
    else {
      for (i = 0; i < estLength; i++) {
	chebNext[i] = 2.0*y[i]*cheb[i] - chebPrev[i];
      }
    }
    
    for (i = 0; i < estLength; i++) {
      chebPrev[i] = cheb[i];
      cheb[i] = chebNext[i];
    }
  }
  free(y);
  free(cheb);
  free(chebPrev);
  free(chebNext);
}

void coeffsToEst2(int estLength, int estStart, double *coeffs, 
		 int numCoeffs, double *est) {
  double *y, *da, *db, *dc;
  int i,k;
  
  y = (double *)malloc(sizeof(double)*estLength);
  da = (double *)malloc(sizeof(double)*estLength);
  db = (double *)malloc(sizeof(double)*estLength);
  dc = (double *)malloc(sizeof(double)*estLength);
  for (i = 0; i < estLength; i++) {
    y[i] = (double)i/(double)estLength*2.0-1.0;
    da[i] = 0.0;
    db[i] = 0.0;
    dc[i] = 0.0;
  }
  for (k = numCoeffs-1; k>=1; k--) {
    for (i = 0; i < estLength; i++) {
      dc[i] = 2*y[i]*db[i]-da[i]+coeffs[k];
      da[i] = db[i];
      db[i] = dc[i];
    }
  }
  for (i = 0; i < estLength; i++) {
    est[i+estStart] = y[i]*dc[i]-da[i]+coeffs[0];
  }
  free(y);
  free(da);
  free(db);
  free(dc);
}

void coeffsToEst(int estLength, int estStart, double *coeffs, 
		 int numCoeffs, double *est) {
  double *cheb, *chebPrev, *chebNext, *y;
  int i,k;

  cheb = (double *)malloc(sizeof(double)*estLength);
  chebNext = (double *)malloc(sizeof(double)*estLength);
  chebPrev = (double *)malloc(sizeof(double)*estLength);
  y = (double *)malloc(sizeof(double)*estLength);

  for (i = 0; i < estLength; i++) {
    cheb[i] = 1.0;
    est[i+estStart] = 0.0;
    y[i] = (double)i/(double)estLength*2.0-1.0;
  }
  for (k = 0; k < numCoeffs; k++) {
    for (i = 0; i < estLength; i++) {
      est[i+estStart] += cheb[i]*coeffs[k];
    }
    if (k==0) {
      for (i = 0; i < estLength; i++) {
	chebNext[i] = y[i];
      }
    }
    else {
      for (i = 0; i < estLength; i++) {
	chebNext[i] = 2*y[i]*cheb[i] - chebPrev[i];
      }
    }
    
    for (i = 0; i < estLength; i++) {
      chebPrev[i] = cheb[i];
      cheb[i] = chebNext[i];
    }
  }
  free(y);
  free(cheb);
  free(chebPrev);
  free(chebNext);
}

void treeToEst(struct tree *T, double *est, 
	       int estLength, int estStart, int level)
{
  int breakpoint;
  
  if (T->optPolyOrder == -1){
    /* split */
    breakpoint = (int) floor(estLength/2);
    treeToEst(T->leftChild, est, breakpoint, estStart,level+1);
    treeToEst(T->rightChild, est, estLength-breakpoint, 
	      estStart+breakpoint,level+1);
  }
  else {
    /* prune */
    coeffsToEst(estLength,estStart,T->polyCoeffs[T->optPolyOrder-1],
		T->optPolyOrder,est);
  }
}

int countObsBelowThresh(double *x, int startingPoint,
			int n, double thresh)
{
  int i,count=0;
  for (i=startingPoint; i < (n+startingPoint); i++)
    {
      if (x[i]<=thresh)
	count++;
      else
	i = n+startingPoint;
    }
  /*if ((x[startingPoint] < 83)&&(x[startingPoint] > 82))
    {
    printf("? x[sP] = %f <= %f, sP = %d, n = %d, count = %d\n",
    x[startingPoint],thresh,startingPoint,n,count);
    }*/
  return count;
}

void fitPolysOnTree(struct tree *T, int r, double *x)
{
  double *cheb, *chebPrev, *chebNext, *estEval, *est, *y;
  int i,k,estLength = 2000;
  double a,b,sum,estMin,renorm, intLength, intLeft, intRight;
	
  intLeft = T->intervalLeft;
  intRight = T->intervalRight;
  intLength = intRight-intLeft;
  T->maxPolyOrder = max(1,min(r,T->nObs));
	
  T->likelihoods = (double *)malloc(sizeof(double)*T->maxPolyOrder);
  T->polyCoeffs = (double **)malloc(sizeof(double*)*T->maxPolyOrder);
	
  if (T->maxPolyOrder == 1) {
    T->polyCoeffs[0] = (double *)malloc(sizeof(double)*1);
    T->polyCoeffs[0][0] = (double)T->nObs/intLength;
    T->likelihoods[0] = T->nObs*(-log(T->polyCoeffs[0][0]));
  }
  else {
    y = (double *)malloc(sizeof(double)*T->nObs);
    cheb = (double *)malloc(sizeof(double)*T->nObs);
    chebNext = (double *)malloc(sizeof(double)*T->nObs);
    chebPrev = (double *)malloc(sizeof(double)*T->nObs);
    estEval = (double *)malloc(sizeof(double)*T->nObs);
    est = (double *)malloc(sizeof(double)*estLength);
    for (i = 0; i < T->nObs; i++) {
      y[i] = (x[i+T->firstObsInInterval]-intLeft)/intLength*2.0-1.0;
      cheb[i] = 1.0;
      estEval[i] = 0.0;
    }
    
    /* calculate unconstrained coefficients */
    for (k = 0; k < T->maxPolyOrder; k++) {
			
      T->polyCoeffs[k] = (double *)malloc(sizeof(double)*(k+1));
      if (k > 0) {
	for (i = 0; i <= k; i++) {
	  T->polyCoeffs[k][i] = T->polyCoeffs[k-1][i];
	}
      }
      T->polyCoeffs[k][k] = 0.0;
      for (i = 0; i < T->nObs; i++) {
	if (pow(y[i],2.0) != 1.0) {
	  T->polyCoeffs[k][k] += cheb[i]*pow(1-pow(y[i],2.0),-0.5);
	}
      }
      
      if (T->nObs>0)
	T->polyCoeffs[k][k] *= 2.0/PI/T->nObs;
      else
	T->polyCoeffs[k][k] = 0.0;
      
      if (k==0) {
	T->polyCoeffs[k][k] /= 2.0;
	for (i = 0; i < T->nObs; i++) {
	  chebNext[i] = y[i];
	}
      }
      else {
	for (i = 0; i < T->nObs; i++) {
	  chebNext[i] = 2*y[i]*cheb[i] - chebPrev[i];
	}
      }
      for (i = 0; i < T->nObs; i++) {
	chebPrev[i] = cheb[i];
	cheb[i] = chebNext[i];
      }
    }
		
    /* renormalize coefficients */
    if (T->nObs > 0) {
      for (k=0; k < T->maxPolyOrder; k++) {
	/* compute min */
	coeffsToEst2(estLength,0,T->polyCoeffs[k],k+1,est);
	estMin = realmax;
	for (i=0; i<estLength; i++) {
	  estMin = min(estMin,est[i]);
	}
				
	/* change min */
	if (estMin < 0) {
	  a = 1.0/(1.0-2.0*estMin);
	  b = (1.0-a)/2.0;
	  for (i=0; i<=k; i++) {
	    T->polyCoeffs[k][i] *= a;
	  }
	  T->polyCoeffs[k][0] += b;
	}
				
				
	/* compute integral */
	/*
	  coeffsToEst(estLength,0,T->polyCoeffs[k],k+1,est);
	  sum = 0;
	  for (i=0; i<estLength; i++) {
	  sum += est[i]/((double)estLength);
	  }
				 
	*/
	sum = T->polyCoeffs[k][0];
	for (i=2; i < (k+1); i+=2) {
	  sum -= T->polyCoeffs[k][i]/((double)i+1)/((double)i-1);
	}
				
	/* change integral */
	if (sum > 0){
	  renorm = ((double)T->nObs/sum/intLength);
	  for (i=0; i<=k; i++) {
	    T->polyCoeffs[k][i] *= renorm;
	  }
	}
				
	if ((intLeft==0)&(intRight==1)) {
	  sum = T->polyCoeffs[k][0];
	  for (i=2; i < (k+1); i+=2) {
	    sum -= T->polyCoeffs[k][i]/((double)i+1)/((double)i-1);
	  }
	}
	/*coeffsToEval(x,T,k+1,estEval);*/
	coeffsToEval2(x,T,k+1,estEval);
				
	T->likelihoods[k] = 0;
	for (i = 0; i < T->nObs; i++) {
	  T->likelihoods[k] -= log(max(realmin,estEval[i]));
	}
      }
    }
    else {
      T->likelihoods[0] = 0.0;
    }
    free(cheb);
    free(chebNext);
    free(chebPrev);
    free(estEval);
    free(y);
    free(est);
  }
	
  if (T->leftChild != NULL)
    fitPolysOnTree(T->leftChild,r,x);
  if (T->rightChild != NULL)
    fitPolysOnTree(T->rightChild,r,x);
	
}



void addNode(double newIntervalLeft,
	     double newIntervalRight, int newNObs,
	     int newFirstObsInInterval, double minWidth,
	     double *x, struct tree *newNode)
{
  double breakpoint;
  int nObsL, nObsR;
	
  /*  newNode = (struct tree *) malloc(sizeof(struct tree));*/
  newNode->intervalLeft = newIntervalLeft;
  newNode->intervalRight = newIntervalRight;
  newNode->nObs = newNObs;
  newNode->firstObsInInterval = newFirstObsInInterval;
  newNode->maxPolyOrder = 0;
	
  if (((newIntervalRight-newIntervalLeft) >= minWidth)&&(newNObs>0))  {
    breakpoint = (newIntervalRight-newIntervalLeft)/2.0+newIntervalLeft;
    nObsL = countObsBelowThresh(x,newFirstObsInInterval,
				newNObs,breakpoint);
    nObsR = newNObs-nObsL;
		
    newNode->leftChild = (struct tree *) malloc(sizeof(struct tree));
    newNode->rightChild = (struct tree *) malloc(sizeof(struct tree));
    //printf("add left\n");
    addNode(newIntervalLeft,breakpoint,
	    nObsL,newFirstObsInInterval,
	    minWidth,x,newNode->leftChild);
    //printf("add right\n");
    addNode(breakpoint,newIntervalRight,
	    nObsR,newFirstObsInInterval+nObsL,
	    minWidth,x,newNode->rightChild);
  }
  else {
    newNode->leftChild = NULL;
    newNode->rightChild = NULL;
  }
  //printf("up ");
  return;
}

void DensityEstGrowTree(struct tree *T, double *x, int n, 
			int maxPolyOrder, int estLength)
{
  double minWidth = max(1.0/estLength,1.0/n);
  /*double minWidth = 0.5;*/
  
  std::cout << "Adding node.." << std::endl;
  
  addNode(0.0,1.0,n,0,minWidth,x,T);
  
  std::cout << "done" << std::endl;
  
  std::cout << "Fitting polynomial..." << std::endl;
  fitPolysOnTree(T, maxPolyOrder, x);
  std::cout << "done" << std::endl;
	
  /*printTree(T,0,x);*/
}	

void DensityEstPruneTree(struct tree *T, double *x, 
			 int n, double pen, int level)
{
  double Lsplit = 0;
  double L,Lmin = realmax;
  int k,i,kmin = 0;
	
  double *estEval;
  
  if (T == NULL) 
    return;
  if (T->leftChild != NULL) {
    DensityEstPruneTree(T->leftChild,x,n,pen,level+1);
    Lsplit += T->leftChild->optLikelihood;
  }
  else {
    Lsplit = realmax;
  }
  if (T->rightChild != NULL) {
    DensityEstPruneTree(T->rightChild,x,n,pen,level+1);
    Lsplit += T->rightChild->optLikelihood;
  }
  else {
    Lsplit = realmax;
  }
	
  for (k = 0; k<T->maxPolyOrder; k++) {
    L = T->likelihoods[k] + (pen*((double)(k+1)/2.0)+(2.0+(double)(k+1))*log(2.0));
    if (L < Lmin) {
      Lmin = L;
      kmin = k;
    }
  }
  
  if (Lsplit < Lmin) {
    T->optLikelihood = Lsplit;
    T->optPolyOrder = -1;
  }
  else {
    T->optLikelihood = Lmin;
    T->optPolyOrder = kmin+1;
  }
}	


void computeDensity(double *x, int n, double *y, int estLength, double pen, int maxPolyOrder)
{
  double renorm;
  int i;
  struct tree *T;
  
  std::cout << "Building tree..." << std::endl;
  
  T = (struct tree *) malloc(sizeof(struct tree));
  DensityEstGrowTree(T,x,n,maxPolyOrder,estLength);
  
  std::cout << "Done building tree" << std::endl;
  
  std::cout << "Pruning tree..." << std::endl;
  DensityEstPruneTree(T,x,n,pen,0);
  std::cout << "done pruning tree" << std::endl;
  
  treeToEst(T,y,estLength,0,0);
  renorm = 1.0/T->nObs;
  for (i = 0; i < estLength; i++) {
    y[i] *= renorm;
  }
  freeTree(T);
}
#ifdef MATLAB_CODE

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  double *x, *y;
  int m,n,i,nDim,estlength;
  int estLength, maxPolyOrder;
  double pen;
  
  /* check for correct # of input variables */
  if (nrhs>4){
    mexErrMsgTxt("There are at most 4 input parameters allowed!");
    return;
  }
  if (nrhs<1){
    mexErrMsgTxt("There is at least 1 input parameter required!");
    return;
  }
  x = mxGetPr(prhs[0]);
  nDim = mxGetNumberOfDimensions(prhs[0]);
	
  if (nDim > 1) {
    m = (int)(mxGetDimensions(prhs[0]))[0];
    n = (int)(mxGetDimensions(prhs[0]))[1];
    if ((m!=1)&&(n!=1)) {
      mexErrMsgTxt("The input array cannot have more than one dimension.");
      return;
    }
  }
	
  n = mxGetNumberOfElements(prhs[0]);
	
  
  if ((nrhs < 2)||(mxGetN(prhs[1])==0)) {
    estLength = (double) pow(2.0,floor(log((double)n)/log(2.0)));
  }
  else {
    estLength = (double) *mxGetPr(prhs[1]);
  }
	
	
	
  /* set penalty
     default = log(n)/2 
  */
  if ((nrhs < 3)||(mxGetN(prhs[2])==0)) {
    pen = (double) log(n)/2.0;
  }
  else {
    pen = (double) *mxGetPr(prhs[2]);
  }
  if (pen < 0) {
    mexErrMsgTxt("The penalty must be positive.");
  }
	
  if ((nrhs < 4)||(mxGetN(prhs[3])==0)) {
    maxPolyOrder = n;
  }
  else {
    maxPolyOrder = (double) *mxGetPr(prhs[3]);
  }
	
  plhs[0] = mxCreateDoubleMatrix(estLength,1,mxREAL);
  y = mxGetPr(plhs[0]);
  computeDensity(x,n,y,estLength,pen,maxPolyOrder);
}

#endif

/* Wrapper by GV */

typedef std::vector<double> doubleArray;

inline std::vector< double > to_std_vector( const boost::python::api::object& iterable )
{
    return std::vector< double >( boost::python::stl_input_iterator< double >( iterable ),
                             boost::python::stl_input_iterator< double >( ) );
}

doubleArray getDensity( boost::python::object& iterable, int n, int estLength, double pen, int maxPolyOrder)
{
  
  doubleArray x = to_std_vector( iterable );
  
  double *xx = &x[0];
    
  doubleArray y( estLength, 0.0 );
  
  double *yy = &y[0];
  
  int i;
  
  for(i=0;i < n;++i)
  {
      printf("%f", xx[i]);
  }
  printf("\n");
  
  printf("%i\n", n);
  
  for(i=0;i < estLength;++i)
  {
      printf("%f", yy[i]);
  }
  printf("\n");
  
  printf("%i\n", estLength);
  printf("%f\n", pen);
  printf("%i\n", maxPolyOrder); 
  
  
  
  computeDensity(xx, n, yy, estLength, pen, maxPolyOrder);
  
  return y;

}


BOOST_PYTHON_MODULE(freeDegreeEst)
{
    using namespace boost::python;
    
    class_<doubleArray>("doubleArray")
        .def(vector_indexing_suite<doubleArray>() );
    
    def("getDensity", getDensity);

}

