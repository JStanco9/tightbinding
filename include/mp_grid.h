// John Stanco 10/23/18

#include <vector>
#include <armadillo>

/*
	Monkhorst-Pack grid storage:
	https://journals.aps.org/prb/abstract/10.1103/PhysRevB.13.
 	https://tu-freiberg.de/sites/default/files/media/institut-fuer-theoretische-physik-10451/Lehre/Dichtefunktionaltheorie/physrevb.13.5188.pdf
*/

#ifndef TB_MP_GRID_H
#define TB_MP_GRID_H

class Lattice;

class MPGrid : public std::vector<arma::vec> {
public:
	explicit MPGrid( const Lattice &lat, size_t q );
	explicit MPGrid( const Lattice &lat, std::vector<size_t> q );
	MPGrid( const MPGrid &other ) = delete;
	MPGrid &operator=( const MPGrid &other ) = delete;
	double volume() const;
private:
	arma::vec toVecIndices( size_t idx ) const;
	void buildGrid( const Lattice &lat );
	void computeStride();
	arma::uvec dims;
	arma::uvec strideLen;
	double vol;
};

#endif /* TB_MP_GRID_H */
