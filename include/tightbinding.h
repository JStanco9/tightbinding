// John Stanco 10/23/18

/*
			Software Package designed to compute and print the band structure of an arbitry crystal solid
			using hopping energy values calculated by the Wannier90 software package.
*/

#include <armadillo>
#include <cmath>

#ifndef TIGHTBINDING_H
#define TIGHTBINDING_H

class BandStructure;
class W90FileReader;
class BZIntegrator;
class Lattice;
class Operator;
class GapMinimizer;


typedef std::complex<double> cx_double;


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
	MPGrid( const MPGrid &other ) = delete;
	MPGrid &operator=( const MPGrid &other ) = delete;
	double volume() const;
private:
	arma::vec toVecIndices( size_t idx );
	void buildGrid( const Lattice &lat );
	void computeStride();
	arma::uvec dims;
	arma::uvec strideLen;
	double vol;
};

#endif /* TB_MP_GRID_H */
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
	//double operator()( const arma::vec &k ) {
		//do something with tb.Ham( k );
		// auto eigs = tb.eigVals( k );
		// need to access fermi energy
		// This way we can determine the gap...
	//}
private:
	const Model &tb;
};


// consider adding
template<class PointGrid, class Func>
double BZIntegrate( const PointGrid &grid, Func f ) {
	double sum = 0;
	for( const auto& k : grid ) {
		sum += f( k );
	}
	return sum * grid.volume() / grid.size();
}

/*
template<class Model>
double injectionCurrent( Model &tb, const MPGrid &grid ) {
	return BZIntegrate( grid, InjectionComponent{ tb } );
}

void print( const BandStructure &bands, std::ofstream& ofs );
*/
/*
	TBModel tb{ "TaAs" };
	// must define unit-cell of lattice

	auto current = BZIntegrate(  )


*/
//gives BZIntegrate( grid, Func{ TBModel{ "TaAs" } } );

#endif /* TIGHTBINDING_H */
