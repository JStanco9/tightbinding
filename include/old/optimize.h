// John Stanco 11/4/18



// This method optimizes maps from vector spaces to a set X
// with respect to the class Comparator.
template<class Func, class T = double, class Comparator = std::less<T>>
class FnOptimizer {
  Func f;
  Comparator cmpr;
	arma::vec centroid;
  arma::vec xr, xe, xc;
  std::vector<std::pair<arma::vec, T>> smplx;
  double fr, fe, fc;
public:
	FnOptimizer( Comparator cmpr = std::less<T>{} );
	FnOptimizer( Func func, Comparator cmpr = std::less<T>{} );
	std::pair<arma::vec, T> optimize( Func func, Comparator cmpr, arma::vec x, double etol = 1e-14 );
	std::pair<arma::vec, T> optimize( Func func, arma::vec x, double etol = 1e-14 );
	std::pair<arma::vec, T> optimize( arma::vec x, double etol = 1e-14 );
private:
	void reflect();
	void expand();
	void contract();
	void shrink();
  void sortSmplxByVal();
};


template<class Func, class T, class Comparator>
FnOptimizer<Func, T, Comparator>::FnOptimizer( Comparator c ) :
	cmpr{ c } {}


template<class Func, class T, class Comparator>
FnOptimizer<Func, T, Comparator>::FnOptimizer( Func func, Comparator c ) :
  f{ func }, cmpr{ c } {}



template<class Func, class T, class Comparator>
void FnOptimizer<Func, T, Comparator>::sortSmplxByVal(){
  std::sort( smplx.begin(), smplx.end(), [this]( const std::pair<arma::vec, T> &lhs, const std::pair<arma::vec, T> &rhs ) {
    return cmpr( lhs.second, rhs.second );
  } );
}


template<class Func, class T, class Comparator>
void FnOptimizer<Func, T, Comparator>::reflect() {
  xr = centroid + ( centroid - smplx.back().first );
  fr = f( xr );
}


template<class Func, class T, class Comparator>
void FnOptimizer<Func, T, Comparator>::expand() {
  xe = centroid + 2 * ( xr - centroid );
  fe = f( xe );
  if( cmpr( fe, fr ) ) {
    smplx.back().first = xe;
    smplx.back().second = fe;
  } else {
    smplx.back().first = xr;
    smplx.back().second = fr;
  }
}


template<class Func, class T, class Comparator>
void FnOptimizer<Func, T, Comparator>::contract() {
  xc = centroid + ( .5 * ( smplx.back().first - centroid ) );
  fc = f( xc );
}


template<class Func, class T, class Comparator>
void FnOptimizer<Func, T, Comparator>::shrink() {
  for( size_t i = 1; i < smplx.size(); ++i ) {
  	auto xs = smplx.front().first + ( .5 * ( smplx[i].first - smplx.front().first ) );
    smplx[i] = { xs, f( xs ) };
  }
}

const size_t maxIter = 200;
template<class Func, class T, class Comparator>
std::pair<arma::vec, T> FnOptimizer<Func, T, Comparator>::optimize( arma::vec x, double etol ) {
	if( x.size() == 0 ) { throw std::logic_error( "invalid vector entered" ); }
  smplx = constructSimplex<Func, T>( f, x, 1 );
  for( size_t i = 0; i < maxIter; ++i ) {
	 sortSmplxByVal();
	 if( std_dev( smplx ) < etol ) { break; }
	 centroid = computeCentroid( smplx );

    // 3) Reflection step
    reflect();
    if(( cmpr( smplx.front().second, fr ) || smplx.front().second == fr )
    		 && cmpr( fr, smplx[ smplx.size() - 2 ].second ) ) {
      smplx.back().first = xr;
      smplx.back().second = fr;
      continue;
    }

    // 4) Expansion step 
    if( cmpr( fr, smplx.front().second ) ) {
      expand();
      continue;
    } 

    // 5) Contraction step
    contract();
    if( cmpr( fc, smplx.back().second ) ) {
			smplx.back().first = xc;
			smplx.back().second = fc;
			continue;
	  }
    // 6) Shrink 
    shrink();
  }
  sortSmplxByVal();
  return smplx.front();
}


template<class Func, class T, class Comparator>
std::pair<arma::vec, T> FnOptimizer<Func, T, Comparator>::optimize( Func func, Comparator c, arma::vec x, double etol ) {
	if( x.size() == 0 ) { throw std::logic_error( "invalid vector entered" ); }
	f = func;
	cmpr = c;
	return optimize( x, etol );
}


template<class Func, class T, class Comparator>
std::pair<arma::vec, T> FnOptimizer<Func, T, Comparator>::optimize( Func func, arma::vec x, double etol ) {
	if( x.size() == 0 ) { throw std::logic_error( "invalid vector entered" ); }
	f = func;
	return optimize( x, etol );
}