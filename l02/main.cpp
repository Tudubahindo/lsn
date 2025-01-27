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
#include <tuple>
#include <vector>
#define _USE_MATH_DEFINES

static constexpr double D = 24.0 / std::pow(M_PI, 2) - 1.0;

double f(double x) {
    return 0.5 * M_PI * std::cos(0.5 * M_PI * x);
}

double p(double x) {
	return (0.25*M_PI*x + 1.0) * (1.0 - x) / (M_PI/24.0 + 0.5);
}

double g(double x) {
    return f(x) / p(x);
}

typedef std::tuple<int, int, int> discrete_vector;
typedef std::tuple<double, double, double> continuous_vector;

discrete_vector discrete_step(discrete_vector position, int direction, bool plus) {
    switch (direction) {
    case 0:
        std::get<0>(position) += 2 * plus - 1;
        break;

    case 1:
        std::get<1>(position) += 2 * plus - 1;
        break;

    case 2:
        std::get<2>(position) += 2 * plus - 1;
        break;

    default:
        break;
    }
    return position;
}

continuous_vector continuous_step(continuous_vector position, double theta, double phi) {
    std::get<0>(position) += std::sin(theta) * std::cos(phi);
    std::get<1>(position) += std::sin(theta) * std::sin(phi);
    std::get<2>(position) += std::cos(theta);

    return position;
}

double discrete_module(discrete_vector v) {
    return std::pow(std::get<0>(v), 2) + std::pow(std::get<1>(v), 2) + std::pow(std::get<2>(v), 2);
}

double continuous_module(continuous_vector v) {
    return std::pow(std::get<0>(v), 2) + std::pow(std::get<1>(v), 2) + std::pow(std::get<2>(v), 2);
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

    //----------------------------------exercise 02.1---------------------------------------

    static constexpr int total = 1000000;              // one million throws
    static constexpr int blocksize = 1000;             // one-thousand-throws sized blocks
    static constexpr int blocknum = total / blocksize; // one thousand blocks

    std::stringstream buffer;

    double running_avg{};
    double running_avg2{};

    double running_importance{};
    double running_importance2{};

    for (int i = 0; i < blocknum; ++i) {
        double avg{};
        double importance{};
        for (int j = 0; j < blocksize; ++j) {
            double num = rnd.Rannyu();
            avg += f(num);
			importance += g(rnd.AcceptReject(0.0, 1.0, 1.6, p));
        }
        avg /= (double)blocksize;
        importance /= (double)blocksize;

        // std::cerr << avg << std::endl;
        double avg2 = std::pow(avg, 2);
        double importance2 = std::pow(importance, 2);

        running_avg = running_avg * ((double)i / (double)(i + 1)) + avg / (double)(i + 1);
        running_avg2 = running_avg2 * ((double)i / (double)(i + 1)) + avg2 / (double)(i + 1);
        double sigma = std::sqrt((running_avg2 - std::pow(running_avg, 2)) / i);
        if (i == 0)
            sigma = 0;

        running_importance = running_importance * ((double)i / (double)(i + 1)) + importance / (double)(i + 1);
        running_importance2 = running_importance2 * ((double)i / (double)(i + 1)) + importance2 / (double)(i + 1);
        double importance_sigma = std::sqrt((running_importance2 - std::pow(running_importance, 2)) / i);
        if (i == 0)
            importance_sigma = 0;

        buffer << running_avg << "\t" << sigma << "\t" << running_importance << "\t" << importance_sigma << "\n";
    }

    std::ofstream out;
    out.open("output.dat");
    out << buffer.str();
    out.close();

    //----------------------------------exercise 02.2---------------------------------------

	rnd.SetRandom(seed, p1, p2);       // reset rnd
    static constexpr int length = 100; // length of the RW

    std::stringstream d_buffer;
    std::stringstream c_buffer;

    std::vector<discrete_vector> d_pos(total, {0, 0, 0});
    std::vector<continuous_vector> c_pos(total, {0., 0., 0.});

    for (int i = 0; i < length; ++i) {

        double running_d_mod{};
        double running_d_mod2{};
        double running_c_mod{};
        double running_c_mod2{};

        for (int k = 0; k < total; ++k) {
            int direction = static_cast<int>(std::floor(6 * rnd.Rannyu()));
            bool verse = direction % 2;
            direction /= 2;
            d_pos.at(k) = discrete_step(d_pos.at(k), direction, verse);

            double theta = M_PI * rnd.Rannyu();
            double phi = 2 * M_PI * rnd.Rannyu();
            c_pos.at(k) = continuous_step(c_pos.at(k), theta, phi);
        }

        for (int k = 0; k < blocknum; ++k) {
            double d_mod{}, c_mod{};
            for (int j = 0; j < blocksize; ++j) {
                d_mod += discrete_module(d_pos.at(j + blocksize * k));
                c_mod += continuous_module(c_pos.at(j + blocksize * k));
            }
            d_mod /= (double)blocksize;
            c_mod /= (double)blocksize;
            double d_mod2 = std::pow(d_mod, 2);
            double c_mod2 = std::pow(c_mod, 2);

            running_d_mod = running_d_mod * ((double)k / (double)(k + 1)) + d_mod / (double)(k + 1);
            running_d_mod2 = running_d_mod2 * ((double)k / (double)(k + 1)) + d_mod2 / (double)(k + 1);
            double d_sigma = std::sqrt((running_d_mod2 - std::pow(running_d_mod, 2)) / k);
            if (k == 0)
                d_sigma = 0;

            running_c_mod = running_c_mod * ((double)k / (double)(k + 1)) + c_mod / (double)(k + 1);
            running_c_mod2 = running_c_mod2 * ((double)k / (double)(k + 1)) + c_mod2 / (double)(k + 1);
            double c_sigma = std::sqrt((running_c_mod2 - std::pow(running_c_mod, 2)) / k);
            if (k == 0)
                c_sigma = 0;

            d_buffer << i << "\t" << std::sqrt(running_d_mod) << "\t" << 0.5 * d_sigma / std::sqrt(running_d_mod) << "\n";
            c_buffer << i << "\t" << std::sqrt(running_c_mod) << "\t" << 0.5 * c_sigma / std::sqrt(running_c_mod) << "\n";
        }
    }

    std::ofstream d_out;
    d_out.open("discrete_output.dat");
    d_out << d_buffer.str();
    d_out.close();

    std::ofstream c_out;
    c_out.open("continuous_output.dat");
    c_out << c_buffer.str();
    c_out.close();

    rnd.SaveSeed("seed2.out");
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
