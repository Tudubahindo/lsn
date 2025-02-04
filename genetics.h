/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "random.h"
#include "stats.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#define _USE_MATH_DEFINES

#ifndef __genetics__
#define __genetics__

typedef std::tuple<double, double> city;

double d2(city a, city b);

double d1(city a, city b);

struct map {
    std::vector<city> cities;
    SquareMatrix<double> distances;

    map(std::vector<city> c);
    std::vector<long unsigned int> bruteforce();
};

class chromosome {

  private:
    std::vector<long unsigned int> _genes;
    std::shared_ptr<map> _c;
    double _fitness;

  public:
    bool check();

    void swap(int a, int b);

    void scramble(Random &rnd);

    bool is_tangled();

    void untangle(Random &rnd);

    chromosome(map m, bool test, Random &rnd);

    chromosome(map m, std::vector<long unsigned int> in);

    chromosome(chromosome p, chromosome m, int stop, Random &rnd);

    void mutation_random_swap(Random &rnd);

    void mutation_shift(Random &rnd);

    void mutation_permutation(Random &rnd);

    void mutation_inversion(Random &rnd);

    double L();

    double fit_getter() const;

    std::string print() const;

    std::vector<long unsigned int> intprint() const;
};

inline bool operator<(const chromosome &lhs, const chromosome &rhs) { return lhs.fit_getter() < rhs.fit_getter(); }
inline bool operator>(const chromosome &lhs, const chromosome &rhs) { return lhs.fit_getter() > rhs.fit_getter(); }
inline bool operator<=(const chromosome &lhs, const chromosome &rhs) { return lhs.fit_getter() <= rhs.fit_getter(); }
inline bool operator>=(const chromosome &lhs, const chromosome &rhs) { return lhs.fit_getter() >= rhs.fit_getter(); }

int selection(std::vector<chromosome> &population, double num);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
