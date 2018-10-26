// John Stanco 10/23/18

#include "mp_grid.h"
#include "lattice.h"


void MPGrid::computeStride() {
	strideLen(dims.size() - 1) = 1;
	for( int i = dims.size() - 2; i >= 0; --i ) {
		strideLen(i) = strideLen(i + 1) * dims(i + 1);
	}
}


arma::vec MPGrid::toVecIndices( size_t index ) const {
		size_t y_i = 0;
		auto xs = arma::vec( dims.size() );
		for( int i = xs.size() - 1; i > 0; --i ){
			xs(i) = ((index - y_i) / strideLen(i)) % dims(i);
			y_i += xs(i) * strideLen(i);
		}
		xs(0) = ((index - y_i) / strideLen(0));
		return xs;
}


void MPGrid::buildGrid( const Lattice &lat ) {
	size_t i = 0;
	for( auto& k : *this ) {
		k = toVecIndices( i );

		// convert to Monkhorst-Pack indices
		for( size_t j = 0; j < k.size(); ++j ) {
			k(j) = (2*k(j) - dims(j) + 1) / (2*dims(j));
		}

		k = lat.vectors() * k;
		++i;
	}
}


MPGrid::MPGrid( const Lattice &lat, size_t q ) :
	std::vector<arma::vec>( pow( q, lat.dimension() ) ),
	dims{ arma::uvec( lat.dimension() ) },
	strideLen{ arma::uvec( lat.dimension() ) },
	vol{ lat.volume() } {
	dims.fill( q );
	computeStride();
	buildGrid( lat );
}

size_t prod( std::vector<size_t> v ) {
	size_t acc = 1;
	for( auto& x : v ) {
		acc *= x;
	}
	return acc;
}

MPGrid::MPGrid( const Lattice &lat, std::vector<size_t> q ) :
	std::vector<arma::vec>( prod(q) ),
	dims{ arma::uvec( lat.dimension() ) },
	strideLen{ arma::uvec( lat.dimension() ) },
	vol{ lat.volume() } {
	if( q.size() != lat.dimension() ) {
		throw "Improper lattice dimensions";
	}

	for( size_t i = 0; i < q.size(); ++i ) { dims(i) = q[i]; }
	computeStride();
	buildGrid( lat );
}

double MPGrid::volume() const {
	return vol;
}
