// John Stanco 10/23/18

#include <cmath>

#include "lattice.h"

const double pi = std::acos( -1 );

Lattice::Lattice( const std::vector<arma::vec> &vecs ) :
dim{ vecs.size() }, latVecs{ arma::mat( vecs.size(), vecs.size() ) } {
  //check dimensions
  for( size_t i = 0; i < vecs.size(); ++i ) {
    if( vecs[i].size() != vecs.size() ) { throw "improper lattice dimensions"; }
    latVecs.col( i ) = vecs[i];
  }
}


Lattice::Lattice( const arma::mat &vecs ) :
dim{ vecs.n_cols }, latVecs{ vecs } {
  if( vecs.n_cols != vecs.n_rows ) { throw "improper lattice dimensions"; }
}


const arma::mat &Lattice::vectors() const {
  return latVecs;
}


double Lattice::dimension() const {
  return dim;
}


double Lattice::volume() const {
  return det( latVecs );
}


Lattice Lattice::dual() const {
  return Lattice{ 2 * pi * latVecs.i().t() };
}
