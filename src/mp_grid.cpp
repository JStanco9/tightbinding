// John Stanco 10/23/18

#include "tightbinding.h"
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


void MPGrid::buildGrid( const arma::mat &vecs ) {
	size_t i = 0;
	for( auto& k : *this ) {
		k = toVecIndices( i );

		// convert to Monkhorst-Pack indices
		for( size_t j = 0; j < k.size(); ++j ) {
			k(j) = (2*k(j) - dims(j) + 1) / (2*dims(j));
		}

		k = vecs * k;
		++i;
	}
}

//
//	API for determining the number of valence electrons in a crystal system
//	Could simply map the entire periodic table...
//	Function n_valence .. takes in vector of crystals...



// vol of region enclosed by vectors in 1, 2, 3 dimensions.
double volume( const arma::mat &span ) {
  if( span.n_cols == 0 || span.n_rows == 0 ||
      span.n_rows > 3 || span.n_cols > span.n_rows ) { throw "Invalid dimensions entered"; }
  if( span.n_cols == span.n_rows ) { return det( span ); }
  if( span.n_cols == 1 ) { return norm( span.col(0) ); }
  //must be 3x2
  return norm( cross( span.col(0), span.col(1) ) );
}


MPGrid::MPGrid( arma::mat vecs, size_t q ) :
	std::vector<arma::vec>( pow(q, vecs.n_cols) ),
	dims{ arma::uvec( vecs.n_cols ) },
	strideLen{ arma::uvec( vecs.n_cols ) },
	vol{ ::volume( vecs ) } {
		dims.fill( q );
		computeStride();
		buildGrid( vecs );
}


size_t prod( std::vector<size_t> v ) {
	size_t acc = 1;
	for( auto& x : v ) {
		acc *= x;
	}
	return acc;
}


MPGrid::MPGrid( arma::mat vecs, std::vector<size_t> q ) :
	std::vector<arma::vec>( prod(q) ),
	dims{ arma::uvec( vecs.n_cols ) },
	strideLen{ arma::uvec( vecs.n_cols ) },
	vol{ ::volume( vecs ) } {
	if( q.size() != vecs.n_cols ) {
		throw "Improper lattice dimensions";
	}
		for( size_t i = 0; i < q.size(); ++i ) { dims(i) = q[i]; }
		computeStride();
		buildGrid( vecs );
}


MPGrid::MPGrid( const Lattice &lat, size_t q ) :
	MPGrid( lat.vectors(), q ) {}


MPGrid::MPGrid( const Lattice &lat, std::vector<size_t> q ) :
	MPGrid( lat.vectors(), q ) {}


double MPGrid::volume() const {
	return vol;
}
