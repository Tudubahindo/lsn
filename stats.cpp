/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "stats.h"

double covariance(std::vector<double> x, std::vector<double> y) {
    assert(x.size() == y.size() && "to compute covariance the two input vectors need to be the same size\n");
    double avg_xy{}, avg_x{}, avg_y{};
    for (long unsigned int i = 0; i < x.size(); ++i) {
        avg_x += x.at(i);
        avg_y += y.at(i);
        avg_xy += x.at(i) * y.at(i);
    }
    avg_x /= (double)x.size();
    avg_y /= (double)x.size();
    avg_xy /= (double)x.size();

    return avg_xy - (avg_x * avg_y);
}

double autocorrelation(std::vector<double> x, long unsigned int lag) {
    assert(lag < x.size() && "lag needs to be smaller than vector size\n");
    double avg_fb{}, avg_f{}, avg_b{}, avg_ff{}, avg_bb{};
    for (long unsigned int i = 0; i < x.size() - lag; ++i) {
        avg_f += x.at(i);
        avg_ff += std::pow(x.at(i), 2);
        avg_b += x.at(i + lag);
        avg_bb += std::pow(x.at(i + lag), 2);
        avg_fb += x.at(i) * x.at(i + lag);
    }
    avg_f /= (double)(x.size() - lag);
    avg_b /= (double)(x.size() - lag);
    avg_ff /= (double)(x.size() - lag);
    avg_bb /= (double)(x.size() - lag);
    avg_fb /= (double)(x.size() - lag);

    return (avg_fb - (avg_f * avg_b)) / (std::sqrt(avg_ff - std::pow(avg_f, 2)) * std::sqrt(avg_bb - std::pow(avg_b, 2)));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
