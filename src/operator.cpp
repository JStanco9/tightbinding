// John Stanco 10/24/18

#include "tightbinding.h"


Operator::Operator() : arma::cx_mat(), isSet{ false } {}

Operator::Operator( const arma::cx_mat &A ) : arma::cx_mat{ A },
isSet{ false } {}


Operator::Operator( size_t n ) :
	arma::cx_mat( n, n ),
	isSet{ false },
	vecs{ arma::cx_mat( n, n ) },
	vals{ arma::vec( n ) } {}


void Operator::computeEigs() {
	eig_sym( vals, vecs, *this );
}


const arma::cx_mat &Operator::eigVecs() {
	if( isSet ) { return vecs; }
	computeEigs();
	isSet = true;
	return vecs;
}


const arma::vec &Operator::eigVals() {
	if( isSet ) { return vals; }
	computeEigs();
	isSet = true;
	return vals;
}
