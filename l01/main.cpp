/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "../random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char *argv[]) {

    Random rnd;
    int seed[4];
    int p1, p2;
    std::ifstream Primes("Primes");
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    } else
        std::cerr << "PROBLEM: Unable to open Primes" << std::endl;
    Primes.close();

    std::ifstream input("seed.in");
    std::string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else
        std::cerr << "PROBLEM: Unable to open seed.in" << std::endl;

    //----------------------------------exercise 01.1---------------------------------------

    static constexpr int total = 1000000;              // one million throws
    static constexpr int blocksize = 1000;             // one-thousand-throws sized blocks
    static constexpr int blocknum = total / blocksize; // one thousand blocks

    std::stringstream buffer; // output buffer

    double running_avg{};
    double running_avg2{};

    double running_variance{};
    double running_variance2{};

    for (int i = 0; i < blocknum; ++i) {
        double avg{};
        double variance{};
        for (int j = 0; j < blocksize; ++j) {
            double num = rnd.Rannyu();
            avg += num;
            variance += std::pow(
                num - 0.5,
                2); // hardcoded variables are never good, however Im lazy
        }
        avg /= (double)blocksize;
        variance /= (double)blocksize;
        // std::cerr << avg << std::endl;
        double avg2 = std::pow(avg, 2);
        double variance2 = std::pow(variance, 2);

        running_avg = running_avg * ((double)i / (double)(i + 1)) + avg / (double)(i + 1);
        running_avg2 = running_avg2 * ((double)i / (double)(i + 1)) + avg2 / (double)(i + 1);
        double sigma = std::sqrt((running_avg2 - std::pow(running_avg, 2)) / i);
        if (i == 0)
            sigma = 0;

        running_variance = running_variance * ((double)i / (double)(i + 1)) +
                           variance / (double)(i + 1);
        running_variance2 = running_variance2 * ((double)i / (double)(i + 1)) +
                            variance2 / (double)(i + 1);
        double variance_sigma =
            std::sqrt((running_variance2 - std::pow(running_variance, 2)) / i);
        if (i == 0)
            variance_sigma = 0;

        buffer << running_avg << "\t" << sigma << "\t";
        buffer << running_variance << "\t" << variance_sigma << "\n";
    }

    std::ofstream out;
    out.open("output.dat");
    out << buffer.str();
    out.close();

    static constexpr int TOT = 10000;
    static constexpr int M = 100; // 100 intervals
    static constexpr int throws = 100;
    static constexpr int N = TOT / M;

    std::vector<double> chi2;

    std::stringstream chi_buffer;

    for (int j = 0; j < throws; ++j) {
        std::vector<int> n(M, 0); // M dimensional vector of 0s, each element is
                                  // the number of occurrences in its interval
        for (int k = 0; k < TOT; ++k) {
            double num = rnd.Rannyu();
            int index = static_cast<int>(
                std::floor(num * M)); // pick an interval at random
            // std::cerr<<index<<std::endl;
            n.at(index)++; // add an occurence to the lucky interval
        }

        double chi = 0; // compute chi squared
        for (int i = 0; i < M; ++i) {
            chi += std::pow(n.at(i) - N, 2) / N;
        }

        chi2.push_back(chi);
        chi_buffer << chi << "\n";
    }

    std::ofstream chi_out;
    chi_out.open("chi.dat");
    chi_out << chi_buffer.str();
    chi_out.close();

    //----------------------------------exercise 01.2---------------------------------------

    static constexpr int NUM[4] = {1, 2, 10, 100};

    std::stringstream d6_buffer, de_buffer, dl_buffer;

    for (int j = 0; j < TOT; ++j) {
        for (int k = 0; k < 4; ++k) {
            double d6{}, de{}, dl{};
            for (int i = 0; i < NUM[k]; ++i) {
                d6 += std::floor(6 * rnd.Rannyu()) + 1;
                de += rnd.Exponential(1.0);
                dl += rnd.Lorentzian(0.0, 1.0);
            }

            d6 /= NUM[k];
            d6_buffer << d6 << "\t";

            de /= NUM[k];
            de_buffer << de << "\t";

            dl /= NUM[k];
            dl_buffer << dl << "\t";
        }

        d6_buffer << "\n";
        de_buffer << "\n";
        dl_buffer << "\n";
    }

    std::ofstream d6_out, de_out, dl_out;

    d6_out.open("d6.dat");
    d6_out << d6_buffer.str();
    d6_out.close();

    de_out.open("de.dat");
    de_out << de_buffer.str();
    de_out.close();

    dl_out.open("dl.dat");
    dl_out << dl_buffer.str();
    dl_out.close();

    //----------------------------------exercise 01.3---------------------------------------

    static constexpr int d = 6;
    static constexpr int L = 5;

    double running_pi = 0;
    double running_pi2 = 0;

    std::stringstream pi_buffer;

    for (int i = 0; i < blocknum; ++i) {
        int N_hit = 0;
        int N_thr = 0;
        for (int j = 0; j < blocksize; ++j) {
            double x = rnd.Rannyu(); // 0 < x < L (because of symmetry, I only
                                     // consider a quarter of the 2L x 2L box)
            double y = rnd.Rannyu(); // 0 < y < L (symmetry)
            double z = rnd.Rannyu(); // 0 < z < d (only 1 coordinate 4
                                     // transaltion symmetry)
            double r = std::pow(x, 2) + std::pow(y, 2);
            if (r <= 1) { // if r > L, reject (the point isn't in the circle)
                N_thr++;
                r = std::sqrt(r); // I do the sqrt AFTER checking the previous
                                  // condition to save unnecessary sqrts
                if (z * d * r + y * L > d * r) {
                    N_hit++;
                }
            }
        }
        double pi = 2 * L * (double)N_thr / (d * (double)N_hit);
        double pi2 = std::pow(pi, 2);

        running_pi =
            running_pi * ((double)i / (double)(i + 1)) + pi / (double)(i + 1);
        running_pi2 =
            running_pi2 * ((double)i / (double)(i + 1)) + pi2 / (double)(i + 1);
        double sigma = std::sqrt((running_pi2 - std::pow(running_pi, 2)) / i);
        if (i == 0)
            sigma = 0;

        pi_buffer << running_pi << "\t" << sigma << "\n";
    }

    std::ofstream pi_out;
    pi_out.open("pi.dat");
    pi_out << pi_buffer.str();
    pi_out.close();

    rnd.SaveSeed();
    return 0;
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
