/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>

#ifndef __lsn_statistics__
#define __lsn_statistics__

double covariance(std::vector<double> x, std::vector<double> y);

double autocorrelation(std::vector<double> x, long unsigned int lag);

template <typename T>
class SquareMatrix {
    std::vector<T> inner;
    long unsigned int dim;

  public:
    SquareMatrix(long unsigned int d);

    T &operator()(long unsigned int row, long unsigned int column);
};

template <typename T>
SquareMatrix<T>::SquareMatrix(long unsigned int d) : dim(d) {
    inner.resize(dim * dim);
}

template <typename T>
T &SquareMatrix<T>::operator()(long unsigned int row, long unsigned int column) {
    assert(row < dim && column < dim && "out of bounds\n");
    return inner.at(dim * row + column);
}

#endif // __Random_Sim__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
