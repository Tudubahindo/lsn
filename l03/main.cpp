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
#include <random>
#include <sstream>
#include <string>
#include <vector>

double N(double x) {
    return 0.5 * (1 + std::erf(x / std::sqrt(2)));
}

double GBM(double m, double s, double t) {
    return std::exp((m - 0.5 * std::pow(s, 2)) * t + s * N(t));
}

double d1(double T, double t, double r, double s, double K, double S) {
    return (std::log(S /* *GBM(r, s, t)*/ / K) + (r + 0.5 * std::pow(s, 2) * (T - t))) / (s * std::sqrt(T - t));
}

double d2(double T, double t, double r, double s, double K, double S) {
    return d1(T, t, r, s, K, S) - s * std::sqrt(T - t);
}

double call(double T, double t, double r, double s, double K, double S) {
    return S /* *GBM(r, s, t)*/ * N(d1(T, t, r, s, K, S)) - K * std::exp(-r * (T - t)) * N(d2(T, t, r, s, K, S));
}

double put(double T, double t, double r, double s, double K, double S) {
    return S /* *GBM(r, s, t)*/ * (N(d1(T, t, r, s, K, S)) - 1) - K * std::exp(-r * (T - t)) * (N(d2(T, t, r, s, K, S)) - 1);
}

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

    //----------------------------------exercise 03.1---------------------------------------

    double S = 100;
    double K = 100;
    static constexpr double T = 1;
    double r = 0.1;
    double s = 0.25;

    std::cerr << "Black-Scholes\n"
              << call(T, 0, r, s, K, S) << "\n"
              << put(T, 0, r, s, K, S) << "\n";

    static constexpr int total = 1000000;              // one million throws
    static constexpr int blocksize = 1000;             // one-thousand-throws sized blocks
    static constexpr int blocknum = total / blocksize; // one thousand blocks

    std::stringstream buffer;

    double running_call{};
    double running_call2{};

    double running_put{};
    double running_put2{};

    for (int i = 0; i < blocknum; ++i) {
        double avg_c{};
        double avg_p{};
        for (int j = 0; j < blocksize; ++j) {
            double num = S * std::exp((r - std::pow(s, 2) / 2) * T - s * rnd.Gauss(0, T)) - K;
            avg_c += std::max(0.0, num) * std::exp(-r * T);
            avg_p += std::max(0.0, -num) * std::exp(-r * T);
        }
        avg_c /= (double)blocksize;
        avg_p /= (double)blocksize;
        double avg_c2 = std::pow(avg_c, 2);
        double avg_p2 = std::pow(avg_p, 2);

        running_call = running_call * ((double)i / (double)(i + 1)) + avg_c / (double)(i + 1);
        running_call2 = running_call2 * ((double)i / (double)(i + 1)) + avg_c2 / (double)(i + 1);
        double sigma_c = std::sqrt((running_call2 - std::pow(running_call, 2)) / i);
        if (i == 0)
            sigma_c = 0;

        buffer << running_call << "\t" << sigma_c << "\t";

        running_put = running_put * ((double)i / (double)(i + 1)) + avg_p / (double)(i + 1);
        running_put2 = running_put2 * ((double)i / (double)(i + 1)) + avg_p2 / (double)(i + 1);
        double sigma_p = std::sqrt((running_put2 - std::pow(running_put, 2)) / i);
        if (i == 0)
            sigma_p = 0;

        buffer << running_put << "\t" << sigma_p << "\n";
    }

    std::ofstream out;
    out.open("output.dat");
    out << buffer.str();
    out.close();

    //----------------------------------exercise 03.2---------------------------------------

    rnd.SetRandom(seed, p1, p2); // reset rnd

    running_call = 0;
    running_call2 = 0;

    running_put = 0;
    running_put2 = 0;

    buffer.str(std::string());

    static constexpr int instants = 100;
    static constexpr double delta_t = T / instants;

    for (int i = 0; i < blocknum; ++i) {
        double avg_c = 0;
        double avg_p = 0;
        for (int j = 0; j < blocksize; ++j) {
            double S_t = S;
            for (int t = 0; t < instants; ++t) {
                S_t = S_t * std::exp((r - std::pow(s, 2) / 2) * delta_t - s * rnd.Gauss(0, 1) * std::sqrt(delta_t));
            }
            double num = S_t - K;
            avg_c += std::max(0.0, num) * std::exp(-r * T);
            avg_p += std::max(0.0, -num) * std::exp(-r * T);
        }
        avg_c /= (double)blocksize;
        avg_p /= (double)blocksize;
        // std::cerr << avg << std::endl;
        double avg_c2 = std::pow(avg_c, 2);
        double avg_p2 = std::pow(avg_p, 2);

        running_call = running_call * ((double)i / (double)(i + 1)) + avg_c / (double)(i + 1);
        running_call2 = running_call2 * ((double)i / (double)(i + 1)) + avg_c2 / (double)(i + 1);
        double sigma_c = std::sqrt((running_call2 - std::pow(running_call, 2)) / i);
        if (i == 0)
            sigma_c = 0;

        buffer << running_call << "\t" << sigma_c << "\t";

        running_put = running_put * ((double)i / (double)(i + 1)) + avg_p / (double)(i + 1);
        running_put2 = running_put2 * ((double)i / (double)(i + 1)) + avg_p2 / (double)(i + 1);
        double sigma_p = std::sqrt((running_put2 - std::pow(running_put, 2)) / i);
        if (i == 0)
            sigma_p = 0;

        buffer << running_put << "\t" << sigma_p << "\n";
    }

    out.open("output100.dat");
    out << buffer.str();
    out.close();

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
