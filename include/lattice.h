// John Stanco 10/23/18

#include <armadillo>
#include <vector>

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
