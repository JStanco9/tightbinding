// John Stanco 10/25/18

#include "tightbinding.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <thread>

#define MAX_THREAD 10

const double pi = std::acos( -1 );
const cx_double I = { 0, 1 };
const cx_double twoPiI = 2 * pi * I;


int getLineSize( const std::string &hrLine ) {
	double x; int lineSize = 0;
	std::istringstream iss( hrLine );
	while( iss >> x ) {
		++lineSize;
	}
	return lineSize;
}


void TBSpec::getFileInfo( std::ifstream &hrFile ) {
	hrFile.seekg( 0, std::ios::beg );
	std::string hrLine;
	getline( hrFile, hrLine );

	hrFile >> nBands;
	hrFile >> nSites;
	nHopping = nSites * nBands * nBands;
	hrWeights.resize( nSites );

	for( auto &wt : hrWeights ) { hrFile >> wt; }

	dataPos = (int)hrFile.tellg() + 1;
	hrFile.seekg( dataPos );
	getline( hrFile, hrLine );
	latDim = getLineSize( hrLine ) - 4;
}


cx_double round( cx_double x, double cutoff ) {
  double rtmp, itmp;
  if ( norm(x) > cutoff ) {
    rtmp = x.real() / cutoff;
    itmp = x.imag() / cutoff;
    return cx_double{ round( rtmp ), round( itmp ) } * cutoff;
  }
  return {0, 0};
}


LatticeSite TBSpec::getNextSite( std::ifstream &hrFile ) {
	LatticeSite site;
	size_t nB2 = nBands * nBands;
	site.hoppingEnergies.resize( nBands, nBands );
	int r, row, col; double a, b;
	site.coords.resize( latDim );
	for( size_t m = 0; m < nB2; ++m ) {
		for( size_t l = 0; l < latDim; ++l ) {
			hrFile >> r;
			site.coords(l) = r;
		}
		hrFile >> row; hrFile >> col;
		hrFile >> a; hrFile >> b;
		site.hoppingEnergies(row - 1, col - 1) = round( { a, b }, .0005 );
	}
	return site;
}


LatticeSite TBSpec::getNextSite( std::ifstream &hrFile, std::ifstream &wsFile ) {
	auto site = getNextSite( hrFile );
	site.wsWeights.resize( nBands, nBands );
	size_t nB2 = nBands * nBands;
	int x, row, col; arma::Col<int> v( latDim );
	for( size_t m = 0; m < nB2; ++m ) {
		for( size_t i = 0; i < latDim; ++i ) { wsFile >> x; }
		wsFile >> row; wsFile >> col;
		wsFile >> site.wsWeights(row-1, col-1);
		for( size_t i = 0; i < site.wsWeights(row-1, col-1); ++i ) {
			for( size_t j = 0; j < latDim; ++j ) {
				wsFile >> x;
				v( j ) = x;
			}
			site.wsvec.push_back( v );
		}
	}
	return site;
}


void TBSpec::fillSites( std::ifstream &hrFile, std::ifstream &wsFile ) {
	hrFile.seekg( dataPos );
	//get rid of first line of wsvec
	std::string firstLine;
	getline( wsFile, firstLine );
	sites.resize( nSites );
	for( size_t i = 0; i < nSites; ++i ) {
		sites[i] = getNextSite( hrFile, wsFile );
		sites[i].hrWeight = hrWeights[i];
	}
}


void TBSpec::fillSites( std::ifstream &hrFile ) {
	sites.resize( nSites );
	for( size_t i = 0; i < nSites; ++i ) {
		sites[i] = getNextSite( hrFile );
		sites[i].hrWeight = hrWeights[i];
	}
}


class TBHamFactory {
	// private singleton constructor;
	TBHamFactory() {}
public:
	TBHamFactory( const TBHamFactory &other ) = delete;
	TBHamFactory &operator=( const TBHamFactory &other ) = delete;
	TBHam buildHam( const arma::vec &k, const TBSpec &specs, size_t nThr );
	static TBHamFactory &Instance() { static TBHamFactory instance; return instance; }
private:
	TBHam localHam_ws( const LatticeSite &site, const arma::vec &k );
	TBHam localHam_no_ws ( const LatticeSite &site, const arma::vec &k );
	TBHam build_ws( const arma::vec &k, const TBSpec &specs, size_t nThr );
	TBHam build_no_ws( const arma::vec &k, const TBSpec &specs, size_t nThr );
	void fillHam_ws( TBHam *H, const arma::vec *k, const TBSpec *specs, size_t start, size_t end );
	void fillHam_no_ws( TBHam *H, const arma::vec *k, const TBSpec *specs, size_t start, size_t end );

};


TBHam TBHamFactory::localHam_ws( const LatticeSite &site, const arma::vec &k ) {
	TBHam H( site.hoppingEnergies.n_cols );  H.zeros();
	cx_double w; int m = 0;
	for( size_t i = 0; i < site.hoppingEnergies.n_rows; ++i ) {
		for( size_t j = 0; j < site.hoppingEnergies.n_cols; ++j ) {
			w = { 0, 0 };
			for( size_t l = 0; l < site.wsWeights(i, j); ++l ) {
				w += exp( twoPiI * dot( site.coords - site.wsvec[m++], k ) );
			}
			H(i, j) += w * site.hoppingEnergies(i, j) / (double)( site.hrWeight * site.wsWeights(i, j) );
		}
	}
	return H;
}


TBHam TBHamFactory::localHam_no_ws( const LatticeSite &site, const arma::vec &k ) {
	return (arma::cx_mat)( site.hoppingEnergies * exp( twoPiI * dot( site.coords, k ) ) / site.hrWeight );
}


void TBHamFactory::fillHam_ws( TBHam *H, const arma::vec *k, const TBSpec *specs, size_t start, size_t end ) {
	for( size_t i = start; i < end; ++i ) {
		*H += localHam_ws( specs->sites[i], *k );
	}
}


void TBHamFactory::fillHam_no_ws( TBHam *H, const arma::vec *k, const TBSpec *specs, size_t start, size_t end ) {
	for( size_t i = start; i < end; ++i ) {
		*H += localHam_no_ws( specs->sites[i], *k );
	}
}


// supports multi-threading
TBHam TBHamFactory::build_ws( const arma::vec &k, const TBSpec &specs, size_t nThr ) {
	TBHam H( specs.nBands ); H.zeros();
	if( nThr == 1 ) {
		for( auto& site : specs.sites ) {
			H += localHam_ws( site, k );
		}
		return H;
	}
	size_t start, end, nSitesPerProcess = specs.nSites / nThr;
	std::vector<std::thread> threads( nThr );
	for( size_t i = 0; i < nThr; ++i ) {
		start = i * nSitesPerProcess;
		end = ( i + 1 ) * nSitesPerProcess;
		if( specs.nSites - end < nSitesPerProcess ) {
				end = specs.nSites;
		}
		threads[i] = std::thread( &TBHamFactory::fillHam_ws, this, &H, &k, &specs, start, end );
	}

	for( auto &thr : threads ) { thr.join(); }
	return H;
}


// supports multi-threading
TBHam TBHamFactory::build_no_ws( const arma::vec &k, const TBSpec &specs, size_t nThr ) {
	TBHam H( specs.nBands ); H.zeros();
	if( nThr == 1 ) {
		for( auto& site : specs.sites ) {
			H += localHam_no_ws( site, k );
		}
		return H;
	}
	size_t start, end, nSitesPerProcess = specs.nSites / nThr;
	std::vector<std::thread> threads( nThr );
	for( size_t i = 0; i < nThr; ++i ) {
		start = i * nSitesPerProcess;
		end = ( i + 1 ) * nSitesPerProcess;
		if( specs.nSites - end < nSitesPerProcess ) {
				end = specs.nSites;
		}
		threads[i] = std::thread( &TBHamFactory::fillHam_no_ws, this, &H, &k, &specs, start, end );
	}

	for( auto &thr : threads ) { thr.join(); }
	return H;
}


TBHam TBHamFactory::buildHam( const arma::vec &k, const TBSpec &specs, size_t nThr ) {
	return specs.use_ws?
		build_ws( k, specs, nThr ) :
		build_no_ws( k, specs, nThr );
}


TBSpec TBModel::getTBSpecs( const std::string &hr ) {
	std::ifstream hrFile{ hr, std::ifstream::in | std::ios::binary };
	if( !hrFile.is_open() ) { throw "Could not open file " + hr; }
	TBSpec spec;
	specs.hr = hr;
	specs.use_ws = false;
	specs.getFileInfo( hrFile );
	specs.fillSites( hrFile );
	return spec;
}


TBSpec TBModel::getTBSpecs( const std::string &hr, const std::string &wsvec ) {
	std::ifstream hrFile{ hr, std::ifstream::in | std::ios::binary };
	std::ifstream wsFile{ wsvec, std::ifstream::in | std::ios::binary };
	if( !hrFile.is_open() ) { throw "Could not open file " + hr; }
	if( !wsFile.is_open() ) { throw "Could not open file " + wsvec; }
	TBSpec spec;
	spec.hr = hr;
	spec.wsvec = wsvec;
	spec.use_ws = true;
	spec.getFileInfo( hrFile );
	spec.fillSites( hrFile, wsFile );
	return spec;
}


TBModel::TBModel( std::string hr, size_t nThr ) :
	specs{ getTBSpecs( hr ) }, 
	_nThr{ nThr } {
		if( _nThr == 0 || _nThr > MAX_THREAD ) { throw "Invalid number of threads entered"; }

	}


TBModel::TBModel( std::string hr, std::string wsvec, size_t nThr ) :
	specs{ getTBSpecs( hr, wsvec ) }, 
	_nThr{ nThr } {
		if( _nThr == 0 || _nThr > MAX_THREAD ) { throw "Invalid number of threads entered"; }
	}


void TBModel::replaceHam( const arma::vec &k ) const {
	kCurr = k;
	hamCurr = TBHamFactory::Instance().buildHam( k, specs, _nThr );
}


const TBHam &TBModel::Ham( arma::vec k ) const {
	if( approx_equal( k, kCurr, "reldiff", 1e-10 ) ) { return hamCurr; }
	replaceHam( k );
	return hamCurr;
}


const arma::cx_mat &TBModel::eigVecs( arma::vec k ) const {
	if( approx_equal( k, kCurr, "reldiff", 1e-10 ) ) { return hamCurr.eigVecs(); }
	replaceHam( k );
	return hamCurr.eigVecs();
}


const arma::vec &TBModel::eigVals(  arma::vec k ) const {
	if( approx_equal( k, kCurr, "reldiff", 1e-10 ) ) { return hamCurr.eigVals(); }
	replaceHam( k );
	return hamCurr.eigVals();
}


void TBModel::setNThreads( size_t nThr ) {
	if( _nThr == 0 || _nThr > MAX_THREAD ) {
		throw "Invalid number of threads entered";
	} _nThr = nThr;
}


size_t TBModel::nThreads() const {
	return _nThr;
}
