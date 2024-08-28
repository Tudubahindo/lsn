/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Random_Sim__
#define __Random_Sim__

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

// namespace Sim
//{
class Random {

  private:
    const uint64_t m_tot{34522712143931ull};
    uint64_t l_tot, n_tot;

  protected:
  public:
    // constructors
    Random();
    // destructor
    ~Random();
    // methods
    void SetRandom(int *, int, int);
    void SaveSeed();
    double Rannyu(void);
    double Rannyu(double min, double max);
    double Gauss(double mean, double sigma);
};

//}
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
