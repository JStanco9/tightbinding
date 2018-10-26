# tightbinding
Interface for constructing tight-binding Hamiltonian from W90 output data.  Entirely encapsulates W90 implementation to ensure scalability.

Makes use of the following libraries:
  Armadillo C++: for fast Linear Algebra and Eigen-Decomposition.
  MPI ( Message Passing Interface ): for multi-threaded construction of Hamiltonian from hopping energies.
  
To Install: ( installs library to /usr/local/lib, header to /usr/local/include )
  cd ~/.../tightbinding &&
  make install
  
To Use: be sure that $PATH variable contains /usr/local/
  #include <tightbinding>
  compile with -ltightbinding, 
  
