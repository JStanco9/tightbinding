// John Stanco 11/1/18

#ifndef TB_OPTIMIZE_H
#define TB_OPTIMIZE_H

// Generic implementation of Nelder-Mead Simplex Algorithm over some arbitrary vector space
// for arbitrary function


typedef std::vector<arma::vec> Simplex;
double size( const Simplex &s );

// This method optimizes maps from vector spaces to a set X
// with respect to the class Comparator.
template<class Func, class T, class Comparator = std::less<T>>
class FnOptimizer {
  Func f;
  Comparator cmpr;
	arma::vec centroid;
  arma::vec xworst, xbest, xr, xe, xc;
  Simplex smplx;
  std::vector<T> vals;
  double fbest, fworst, fr, fe, fc;
  size_t iter;
public:
	FnOptimizer( Func func, Comparator cmpr );
	arma::vec optimize( arma::vec x, double etol = 1e-8 );
private:
  void init( const arma::vec &x );
	void evaluate();
	void reflect();
	void expand();
	void contract();
	void shrink();
  void sort();
};


template<class Func, class T, class Comparator>
FnOptimizer<Func, T, Comparator>::FnOptimizer( Func func, Comparator c ) :
  f{ func }, cmpr{ c } {}


Simplex constructSimplex( arma::vec v, double r ) {
  int n = v.size();
  Simplex smplx( n + 1 );
  for( auto& pt : smplx ) { pt.zeros(); }
  int tmp = r;
  for( size_t i = n; i > 0; i-- ) {
    smplx[i - 1](i) = tmp;
    for( int j = i; j > 0; j-- ) {
      smplx[i - 1](j - 1) += -tmp / i;
    }
    tmp *= sqrt( i * i - 1 ) / i;
  }
  for( int i = 0; i < n; i++ ) {
    smplx[i] += v;
  }
  return smplx;
}


template<class Func, class T, class Comparator>
void FnOptimizer<Func, T, Comparator>::init( const arma::vec &x ) {
  smplx = constructSimplex( x, 0.005 );
  vals.resize( smplx.size() );

  centroid = arma::zeros<arma::vec>( x.size() );
  xworst = arma::zeros<arma::vec>( x.size() );
  xbest = arma::zeros<arma::vec>( x.size() );
  xr = arma::zeros<arma::vec>( x.size() );
  xe = arma::zeros<arma::vec>( x.size() );
  xc = arma::zeros<arma::vec>( x.size() );

  fworst = 0; fbest = 0; 
  fr = 0; fe = 0; fc = 0;
  iter = 0;

  evaluate();
  sort();
}


// construct factorial table
size_t factorial( size_t n ) {
  size_t tmp = 1;
  for( size_t i = 0; i < n; ++i ) { tmp *= i; }
  return tmp;
}


double volume( const Simplex &smplx ) {
  arma::mat hull( smplx.size(), smplx.size() );
  for( size_t i = 0; i < smplx.size(); ++i ) {
    hull.col( i ) = smplx[i + 1] - smplx[0];
  }
  return det( hull ) / factorial( smplx.size() );
}


template<class Func, class T, class Comparator>
void FnOptimizer<Func, T, Comparator>::evaluate() {
  for( size_t i = 0; i < smplx.size(); ++i ) {
    vals[i] = Func( smplx[i] );
  }
}


template<class Func, class T, class Comparator>
void FnOptimizer<Func, T, Comparator>::sort(){
  std::sort( smplx.begin(), smplx.end(), [this]( const arma::vec &lhs, const arma::vec &rhs ) {
    return cmpr( f( lhs ), f( rhs ) );
  } );

  fbest = f( smplx.front() ); 
  fworst = f( smplx.back() );
}


/*
void FnOptimizer::output()
{
  std::cout << "Iteration #:   " << iter << std::endl;
  std::cout << "Simplex Size:  " << std::setprecision(16)
  << size << std::endl;
  std::cout << "Minimum Gap:   " << std::setprecision(16)
  << f(0, 0) << "\n" << std::endl;
  std::cout << "Simplex Vertices:    " << "\n";
  trans(x).print();
  std::cout << "Function at Vertices " << "\n";
  trans(f.col(0)).print();
  std::cout << "\n" << std::endl;
}
*/


template<class Func, class T, class Comparator>
void FnOptimizer<Func, T, Comparator>::reflect() {
  centroid.zeros();
  for( int i = 0; i < smplx.size(); ++i ) {
    centroid += smplx[i];
  }
  centroid /= ( double )smplx.size();
  xr = centroid + ( centroid - smplx.back() );
  fr = f( xr );
}


template<class Func, class T, class Comparator>
void FnOptimizer<Func, T, Comparator>::expand() {
  xe = centroid + ( 2 * ( xr - centroid ) );
  fe = f( xe );
  if( fe < fr ) {
    smplx.back() = xe;
    vals.back() = fe;
  } else {
    smplx.back() = xr;
    vals.back() = fr;
  }
}


template<class Func, class T, class Comparator>
void FnOptimizer<Func, T, Comparator>::contract() {
  if( fr < vals.back() ) {
    xworst = smplx.back();
    fworst = vals.back();
  } else {
    xworst = xr;
    fworst = fr;
  }
  xc = centroid + ( .5 * ( smplx.back() - centroid ) );
  fc = f( xc );
  
}


template<class Func, class T, class Comparator>
void FnOptimizer<Func, T, Comparator>::shrink() {
  for( size_t i = 0; i < smplx.size(); i++ ) {
    smplx[i] = smplx.front() + ( .5 * ( smplx[i] - smplx.front() ) );
    vals[i] = f( smplx[i] );
  }
}

const size_t maxIter = 1000;

template<class Func, class T, class Comparator>
arma::vec FnOptimizer<Func, T, Comparator>::optimize( arma::vec x, double etol ) {
  init( x );
  while( ( volume( smplx ) > etol || ( fworst - fbest ) > etol ) && iter < maxIter ) {
    //if( !( iter % 10 ) ) { output(); }
    ++iter;
    // 2) Reflection Step
    reflect();
    if( vals.front() <= fr && fr < vals[ vals.size() - 2 ] ) {
      smplx.back() = xr;
      vals.back() = fr;
    } else if( fr < vals.front() ) {
      // 3) Expansion step 
      expand();
    } else if ( fr >= vals[ vals.size() - 2 ] ) {
      // 4) Contraction Step
      contract();
      if( fc < fworst ) {
        smplx.back() = xc;
        vals.back() = fc;
      } else {
        shrink();
      }
    }
    sort();
  }
  return smplx[0];
  //std::cout << "Band minimum " << f(0, 0) << "." << std::endl;
  //std::cout << "Iterations taken: " << iter << "\n" << std::endl;
}

#endif /* TB_OPTIMIZE_H */