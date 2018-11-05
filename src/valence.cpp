// John Stanco 10/29/18

#include "tightbinding.h"
/*

		API for determining the number of valence orbitals of an atom
		Which can be used to obtain the number of valence bands in a 
		crystal system consisting of many atoms.

*/


/*

	In the future, we will use the total number of valence electrons given 
	in the file from PwSCF output.

*/


std::vector<std::vector<int>> levelSizes = {
		{ 2, 0 },
		{ 8, 2, 0 },
		{ 18, 8, 2, 0 },
		{ 32, 18, 8, 2, 0 },
		{ 32, 18, 8, 2, 0 },
		{ 18, 8, 2, 0 },
		{ 8, 2, 0 }
	};


size_t nValence( element elem, bool spinOrb ) {
	if( elem.size() == 0 || elem.size() > 7 ) { throw "Invalid element entered"; }
	//determine size of orbital
	size_t nV = 0;
	for( size_t i = 0; i < elem.size() - 1; ++i ) {
		for( auto &levelSize : levelSizes[i] ) {
			if( elem[i] >= levelSize ) {
				nV += elem[i] - levelSize;
				break;
			}
		}
	}
	nV += elem[elem.size()-1];
	return spinOrb?
		2 * nV :
		nV;
}


size_t nValence( molecule mole, bool spinOrb ) {
	size_t nV = 0;
	for( const auto& elem : mole ) {
		nV += nValence( elem, spinOrb );
	}
	return nV;
}