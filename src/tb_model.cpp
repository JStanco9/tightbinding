// John Stanco 10/25/18

#include "tightbinding.h"
#include <cmath>
#include <fstream>
#include <sstream>

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


void TBFileSpec::fill( std::ifstream &hrFile ) {
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


HamFactory::HamFactory( const TBFileSpec &spec ) :
	specs{ spec } {
		specs.use_ws_distance?
			fillSites_ws() :
			fillSites_no_ws();
	}


Ham HamFactory::localHam_ws( const latticeSite &site, const arma::vec &k ) {
	Ham H( specs.nBands );  H.zeros();
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


Ham HamFactory::localHam_no_ws( const latticeSite &site, const arma::vec &k ) {
	return (arma::cx_mat)( site.hoppingEnergies * exp( twoPiI * dot( site.coords, k ) ) / site.hrWeight );
}


Ham HamFactory::build_ws( const arma::vec &k ) {
	Ham H( specs.nBands ); H.zeros();
	for( const auto &site : sites ) {
		H += localHam_ws( site, k );
	}
	return H;
}


Ham HamFactory::build_no_ws( const arma::vec &k ) {
	Ham H( specs.nBands ); H.zeros();
	for( const auto& site : sites ) {
		H += localHam_no_ws( site, k );
	}
	return H;
}


Ham HamFactory::buildHam( const arma::vec &k ) {
	return specs.use_ws_distance?
		build_ws( k ) :
		build_no_ws( k );
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


latticeSite HamFactory::getNextSite( std::ifstream &hrFile ) {
	latticeSite site;
	size_t nB2 = specs.nBands * specs.nBands;
	site.hoppingEnergies.resize( specs.nBands, specs.nBands );
	int r, row, col; double a, b;
	site.coords.resize( specs.latDim );
	for( size_t m = 0; m < nB2; ++m ) {
		for( size_t l = 0; l < specs.latDim; ++l ) {
			hrFile >> r;
			site.coords(l) = r;
		}
		hrFile >> row; hrFile >> col;
		hrFile >> a; hrFile >> b;
		site.hoppingEnergies(row - 1, col - 1) = round( { a, b }, .0005 );
	}
	return site;
}


latticeSite HamFactory::getNextSite( std::ifstream &hrFile, std::ifstream &wsFile ) {
	auto site = getNextSite( hrFile );

	site.wsWeights.resize( specs.nBands, specs.nBands );
	size_t nB2 = specs.nBands * specs.nBands;

	int x, row, col; arma::Col<int> v( specs.latDim );
	for( size_t m = 0; m < nB2; ++m ) {
		for( size_t i = 0; i < specs.latDim; ++i ) { wsFile >> x; }

		wsFile >> row; wsFile >> col;
		wsFile >> site.wsWeights(row-1, col-1);

		for( size_t i = 0; i < site.wsWeights(row-1, col-1); ++i ) {
			for( size_t j = 0; j < specs.latDim; ++j ) {
				wsFile >> x;
				v( j ) = x;
			}
			site.wsvec.push_back( v );
		}
	}
	return site;
}


void HamFactory::fillSites_ws() {
	std::ifstream hrFile( specs.hr, std::ifstream::in | std::ios::binary );
	std::ifstream wsFile( specs.wsvec, std::ifstream::in | std::ios::binary );
	if( !hrFile.is_open() ) { throw "failed to open file " + specs.hr; }
	if( !wsFile.is_open() ) { throw "failed to open file " + specs.wsvec; }
	hrFile.seekg( specs.dataPos );
	//get rid of first line of wsvec
	std::string firstLine;
	getline( wsFile, firstLine );
	sites.resize( specs.nSites );
	for( size_t i = 0; i < specs.nSites; ++i ) {
		sites[i] = getNextSite( hrFile, wsFile );
		sites[i].hrWeight = specs.hrWeights[i];
	}
}


void HamFactory::fillSites_no_ws() {
	std::ifstream hrFile( specs.hr, std::ifstream::in | std::ios::binary );
	if( !hrFile.is_open() ) { throw "failed to open file " + specs.hr; }
	sites.resize( specs.nSites );
	for( size_t i = 0; i < specs.nSites; ++i ) {
		sites[i] = getNextSite( hrFile );
		sites[i].hrWeight = specs.hrWeights[i];
	}
}


TBModel::TBModel( std::string hr ) :
	factory{ getTBSpecs( hr ) } {}


TBModel::TBModel( std::string hr, std::string wsvec ) :
	factory{ getTBSpecs( hr, wsvec ) } {}


HamFactory TBModel::getTBSpecs( const std::string &hr ) {
	TBFileSpec specs;
	std::ifstream hrFile{ hr, std::ifstream::in | std::ios::binary };
	if( !hrFile.is_open() ) { throw "Could not open file " + hr; }
	specs.hr = hr;
	specs.use_ws_distance = false;
	specs.fill( hrFile );
	return HamFactory{ specs };
}


HamFactory TBModel::getTBSpecs( const std::string &hr, const std::string &wsvec ) {
	TBFileSpec specs;
	std::ifstream hrFile{ hr, std::ifstream::in | std::ios::binary };
	if( !hrFile.is_open() ) { throw "Could not open file " + hr; }
	specs.hr = hr;
	specs.wsvec = wsvec;
	specs.use_ws_distance = true;
	specs.fill( hrFile );
	return HamFactory{ specs };
}


void TBModel::replaceHam( const arma::vec &k ) const {
	kCurr = k;
	hamCurr = factory.buildHam( k );
}


const Ham &TBModel::Ham( arma::vec k ) const {
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
