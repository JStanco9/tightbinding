// John Stanco 10/23/18

/*
			Software Package designed to compute and print the band structure of an arbitry crystal solid
			using hopping energy values calculated by the Wannier90 software package.
*/

#include <armadillo>
#include <cmath>
#include <iomanip>

#ifndef TIGHTBINDING_H
#define TIGHTBINDING_H

class BandStructure;
class W90FileReader;
class BZIntegrator;
class Lattice;
class Operator;
class GapMinimizer;


typedef std::complex<double> cx_double;

#ifndef TB_HAM_H
#define TB_HAM_H


class Operator : public arma::cx_mat {
public:
	Operator();
	explicit Operator( size_t n );
	Operator( const arma::cx_mat &A );
	Operator( size_t m, size_t n ) = delete;
	const arma::cx_mat &eigVecs();
	const arma::vec &eigVals();
private:
	void computeEigs();
	mutable bool isSet;
	arma::cx_mat vecs;
	arma::vec vals;
};


#endif /* TB_HAM_H */

typedef Operator TBHam;

struct LatticeSite {
	arma::Col<int> coords;
	arma::cx_mat hoppingEnergies;
	arma::umat wsWeights;
	std::vector<arma::Col<int>> wsvec;
	size_t hrWeight;
};


class TBSpec {
public:
	size_t nBands;
	size_t nSites;
	size_t nHopping;
	std::vector<size_t> hrWeights;
	std::vector<LatticeSite> sites;
	size_t latDim;
	int dataPos;
	bool use_ws;
	std::string hr;
	std::string wsvec;
	void getFileInfo( std::ifstream &hrFile );
	void fillSites( std::ifstream &hrFile );
	void fillSites( std::ifstream &hrFile, std::ifstream &wsFile );
private:
	LatticeSite getNextSite( std::ifstream &hrFile );
	LatticeSite getNextSite( std::ifstream &hrFile, std::ifstream &wsFile );
};


class TBModel {
public:
	TBModel( std::string hr, std::string wsvec, size_t nThr = 1 );
	TBModel( std::string hr, size_t nThr = 1 );
	const TBHam &Ham( arma::vec k ) const;
	const arma::cx_mat &eigVecs( arma::vec k ) const;
	const arma::vec &eigVals( arma::vec k ) const;
	void setNThreads( size_t nThr );
	size_t nThreads() const;
private:
	TBSpec getTBSpecs( const std::string &hr );
	TBSpec getTBSpecs( const std::string &hr, const std::string &wsvec );
	void replaceHam( const arma::vec &k ) const;
	mutable arma::vec kCurr;
	mutable TBHam hamCurr;
	TBSpec specs;
	size_t _nThr;
};


template<class Model>
class InjectionComponent {
	InjectionComponent( const Model &m ) : tb{ m } {}
	InjectionComponent( const InjectionComponent &other ) = delete;
	InjectionComponent &operator=( const InjectionComponent &other ) = delete;
private:
	const Model &tb;
};


// consider adding
template<class Func, class PointGrid>
double integrate( Func f, const PointGrid &grid ) {
	double sum = 0;
	for( const auto& k : grid ) { sum += f( k ); }
	return sum * grid.volume() / grid.size();
}


template<class Func, class PointGrid, class Stream>
void plot( Func f, const PointGrid &grid, Stream &ofs, size_t prec = 4 ) {
	size_t wid = prec + 6;
	for( const auto &k : grid ) {
		for( const auto &ki : k ) {
			ofs << std::setw( wid ) << std::setprecision( prec ) << k;
		}
		ofs << std::setprecision( prec ) << f( k ) << "\n";
	}
}


// func could be band-gap
template<class Func, class PointGrid>
void plot( Func f, const PointGrid &grid, std::basic_ostream<char> &ofs = std::cout ) {
	for( const auto &k : grid ) {
		for( const auto &ki : k ) {
			ofs << std::setw( 10 ) << std::setprecision( 4 ) << k;
		}
		ofs << std::setprecision( 4 ) << f( k ) << "\n";
	}
}


#ifndef TB_LATTICE_H
#define TB_LATTICE_H


class Lattice {
public:
	Lattice( const std::vector<arma::vec> &vecs );
	Lattice( const arma::mat &vecs );
	const arma::mat &vectors() const;
	double dimension() const;
	double volume() const;
	Lattice dual() const;
private:
	size_t dim;
	arma::mat latVecs;
};

#endif /* TB_LATTICE_H */
#ifndef TB_MP_GRID_H
#define TB_MP_GRID_H


class MPGrid : public std::vector<arma::vec> {
public:
	explicit MPGrid( const Lattice &lat, size_t q );
	explicit MPGrid( const Lattice &lat, std::vector<size_t> q );
	MPGrid( arma::mat vecs, size_t q );
	MPGrid( arma::mat vecs, std::vector<size_t> q );
	MPGrid( const MPGrid &other ) = delete;
	MPGrid &operator=( const MPGrid &other ) = delete;
	double volume() const;
private:
	arma::vec toVecIndices( size_t idx ) const;
	void buildGrid( const arma::mat &vecs );
	void computeStride();
	arma::uvec dims;
	arma::uvec strideLen;
	double vol;
};

#endif /* TB_MP_GRID_H */
#ifndef TB_VALENCE_H
#define TB_VALENCE_H
	
	typedef std::vector<int> element;
	typedef std::vector<element> molecule;
	size_t nValence( element elem, bool spinOrb = false );
	size_t nValence( molecule mole, bool spinOrb = false );

#endif /* TB_VALENCE_H */
#ifndef TB_OPTIMIZE_H
#define TB_OPTIMIZE_H

/*
template<class T> class Vertex : public std::pair<arma::vec, T>{};
template<class T> class Simplex : public std::vector<std::pair<arma::vec, T>>{
	public: Simplex( size_t n ) : std::vector<std::pair<arma::vec, T>>( n ) {}
};
*/

template<class Func, class T>
std::vector<std::pair<arma::vec, T>> constructSimplex( Func f, arma::vec v, double r ) {
  int n = v.size();
  std::vector<std::pair<arma::vec, T>> smplx( n + 1 );
  for( auto& pt : smplx ) { pt.first.resize( n ); pt.first.zeros(); }
  double tmp = r;
  for( size_t i = n; i > 0; i-- ) {
    smplx[i].first(i - 1) = tmp;
    for( int j = i; j > 0; j-- ) {
      smplx[j - 1].first(i - 1) += -tmp / i;
    }
    tmp *= sqrt( i * i - 1 ) / i;
  }
  for( auto &pt : smplx ) { 
  	pt.first += v;
  	pt.second = f( pt.first );
  }
  return smplx;
}


inline size_t factorial( size_t n ) {
  size_t tmp = 1;
  for( size_t i = 1; i <= n; ++i ) { tmp *= i; }
  return tmp;
}


template<class T>
double volume( const std::vector<std::pair<arma::vec, T>> &smplx ) {
  arma::mat hull( smplx.size() - 1, smplx.size() - 1 );
  for( size_t i = 0; i < hull.n_cols; ++i ) {
    hull.col( i ) = ( smplx[i + 1].first - smplx.front().first );
  }
  return det( hull ) / factorial( hull.n_cols );
}


template<class T>
arma::vec computeCentroid( const std::vector<std::pair<arma::vec, T>> &smplx ) {
	arma::vec centroid = arma::zeros<arma::vec>( smplx.size() - 1 );
	for( size_t i = 0; i < smplx.size() - 1; ++i ) {
		centroid += smplx[i].first;
	}
	return centroid / ( smplx.size() - 1 );
}


//compute std::deviation of container
template<class T>
double std_dev( const std::vector<std::pair<arma::vec, T>> &smplx ) {
	//compute mean
	if( smplx.size() == 0 ) { throw std::logic_error( "invalid container" ); }
	if( smplx.size() == 1 ) { return 0; }
	double m = 0, v = 0, dx;
	for( const auto &vtx : smplx ) { m += vtx.second; } 
	m /= smplx.size();
	for( const auto &vtx : smplx ) {
		dx = vtx.second - m;
		v += dx * dx;
	}
	v /= ( smplx.size() - 1 );
	return sqrt( v );
}


// function version
template<class Func, class T = double, class Comparator>
std::pair<arma::vec, T> optimize( Func f, arma::vec k, Comparator cmpr, size_t maxIter, double etol ) {
  if( k.size() == 0 ) { throw std::logic_error( "invalid vector entered" ); }

  // initialize
  auto smplx = constructSimplex<Func, T>( f, k, 1 );
  arma::vec centroid, xr, xe, xc;
  double fr, fe, fc;
  size_t iter = 0;

  while( true ) {
    std::sort( smplx.begin(), smplx.end(), [&]( const std::pair<arma::vec, T> &lhs, const std::pair<arma::vec, T> &rhs ) {
    	return cmpr( lhs.second, rhs.second );
  	} );
 
    if( iter++ > maxIter || std_dev( smplx ) < etol ) { break; }
    centroid = computeCentroid( smplx );

    // 3) Reflection step
    xr = centroid + ( centroid - smplx.back().first );
    fr = f( xr );
    if(( cmpr( smplx.front().second, fr ) || smplx.front().second == fr )
         && cmpr( fr, smplx[ smplx.size() - 2 ].second )) {
      smplx.back() = { xr, fr };
      continue;
    }

    // 4) Expansion step 
    if( cmpr( fr, smplx.front().second ) ) {
      xe = centroid + 2 * ( xr - centroid );
      fe = f( xe );
      if( cmpr( fe, fr ) ) {
        smplx.back() = { xe, fe };
      } else {
        smplx.back() = { xr, fr };
      }
      continue;
    } 

    // 5) Contraction step
    xc = centroid + ( .5 * ( smplx.back().first - centroid ) );
    fc = f( xc );
    if( cmpr( fc, smplx.back().second ) ) {
      smplx.back() = { xc, fc };
      continue;
    }

    // 6) Shrink 
    arma::vec xs;
    for( size_t i = 1; i < smplx.size(); ++i ) {
    	xs = smplx.front().first + ( .5 * ( smplx[i].first - smplx.front().first ) );
      smplx[i] = { xs, f( xs ) };
    }
  }
  return smplx.front();
}


template<class Func, class T = double, class Comparator = std::less<T>>
std::pair<arma::vec, T> optimize( Func f, arma::vec k, size_t maxIter, double etol ) {
  return optimize( f, k, Comparator{}, maxIter, etol );
}


/* TODO:

	1. Band Gap
	2. Band structure
	3. Fermi's Golden Rule
	4. Function Optimizer ( DONE )

*/


#endif /* TB_OPTIMIZE_H */
#endif /* TIGHTBINDING_H */
