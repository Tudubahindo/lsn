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
#include "../stats.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <tuple>
#include <vector>
#define _USE_MATH_DEFINES

static constexpr double initial_mu = 1.0;
static constexpr double initial_sigma = 1.0;

struct walker {

    double x;

    void step(double a, double xi) { x += a * xi; }
};

double trial(walker w, double mu, double s2) {
    double psi = std::exp(-std::pow(w.x - mu, 2) / (2 * s2)) +
                 std::exp(-std::pow(w.x + mu, 2) / (2 * s2));
    return std::pow(psi, 2);
}

double Acceptance(walker v, walker w, double mu, double s2) {
    return std::min(1.0, trial(v, mu, s2) / trial(w, mu, s2));
}

//---------------------------------------------------------------------------------------------------

std::tuple<double, double> Energy(double step_length, double mu, double sigma,
                                  Random &rnd, bool calibration,
                                  bool verbose) {

    static constexpr int total = 10000;                // ten thousand rw steps
    static constexpr int blocksize = 100;              // one hundred blocks
    static constexpr int blocknum = total / blocksize; // one hundred steps per block

	std::stringstream buffer, hist_buffer;
    walker sampler{};

    double running_H{};
    double running_H2{};
    double sigma_H{};
	std::vector<double> H_ac;
    long rate{};

    for (int i = 0; i < blocknum; ++i) {
        double H{};
        for (int j = 0; j < blocksize; ++j) {

            walker slave = sampler;
            double xi = 2 * rnd.Rannyu() - 1;
            slave.step(step_length, xi);
            double num = rnd.Rannyu();
            if (num < Acceptance(slave, sampler, mu, std::pow(sigma, 2))) {
                sampler = slave;
                rate++;
            }
            H += std::pow(sampler.x, 4) - 2.5 * std::pow(sampler.x, 2) -
                 0.5 *
                     ((std::pow(sampler.x, 2) + std::pow(mu, 2)) -
                      2 * sampler.x * mu *
                          std::tanh(sampler.x * mu / std::pow(sigma, 2)) -
                      std::pow(sigma, 2)) /
                     std::pow(sigma, 4);
			if (verbose == true && calibration == false) 
				hist_buffer << sampler.x << "\n";
			if (verbose == true && calibration == true) 
				H_ac.push_back(sampler.x);
        }

        H /= (double)blocksize;
        double H2 = std::pow(H, 2);

        running_H = running_H * ((double)i / (double)(i + 1)) + H / (double)(i + 1);
        running_H2 = running_H2 * ((double)i / (double)(i + 1)) + H2 / (double)(i + 1);
        sigma_H = std::sqrt((running_H2 - std::pow(running_H, 2)) / i);
        if (i == 0)
            sigma_H = 0;
		if (verbose == true)
			buffer << running_H << "\t" << sigma_H << "\n";
    }

	if (verbose == true) {
	    std::ofstream out;
		std::string name = "output_";
		if (calibration == true)
			name += "calibration.dat";
		else
			name += "final.dat";
		out.open(name);
	    out << buffer.str();
		out.close();

		if (calibration == false) {
			out.open("output_hist.dat");
			out << hist_buffer.str();
			out.close();
		}

		if (calibration == true) {
			out.open("output_autoc.dat");
			for (int i=0; i<blocksize; ++i){
				out << autocorrelation(H_ac, i) << "\n";
			}
			out.close();
		}
	}

    if (calibration == false) {
        return {running_H, sigma_H};
    } else {
        return {running_H, (double)rate / (double)total};
    }
}

struct param_space_walker {

    double s;
    double m;
    double L;
    double l;

    void step(double a, double xi, double ypsilon, double step_length,
              Random &rnd) {
        s += a * xi;
        m += a * ypsilon;
        std::tie(L, l) = Energy(step_length, m, s, rnd, false, false);
    }

    void recompute(double step_length, Random &rnd) {
        std::tie(L, l) = Energy(step_length, m, s, rnd, false, false);
    }

    param_space_walker(double step_length, Random &rnd) {
        s = initial_sigma;
        m = initial_mu;
        std::tie(L, l) = Energy(step_length, m, s, rnd, false, false);
    }
};

double param_space_acceptance(param_space_walker v, param_space_walker w,
                              double beta) {
    return std::min(1.0, std::exp(-beta * (v.L - w.L)));
}

//---------------------------------------------------------------------------------------------------

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

    //---------------------------------------------------------------------------------------------------

    double step_length = 0.1;
    double best_acceptance{};

    std::vector<walker> wavefunction{};
	std::stringstream buffer;

    for (int k = 1; k < 51; ++k) {
        double rate = std::get<1>(
            Energy(0.1 * k, initial_mu, initial_sigma, rnd, true, false));
        if (k == 1) {
            best_acceptance = rate;
        } else {
            if (std::abs(rate - 0.5) < std::abs(best_acceptance - 0.5)) {
                best_acceptance = rate;
                step_length = 0.1 * k;
            }
        }
    }

	auto trash = Energy(step_length, initial_mu, initial_sigma, rnd, true, true);
    std::cerr << step_length << "\n";

    param_space_walker parameters_sampler =
        param_space_walker(step_length, rnd);

    double beta = 0.001;
    static constexpr int sameTsteps = 100;
    static constexpr int doublings = 20;
    static constexpr double a = 0.2; // manual fine tuning

    for (int i = 0; i < doublings; ++i) {
        for (int j = 0; j < sameTsteps; ++j) {
            param_space_walker slave = parameters_sampler;
            parameters_sampler.recompute(step_length, rnd);
            double xi = 2 * rnd.Rannyu() - 1;
            double ypsilon = 2 * rnd.Rannyu() - 1;
            slave.step(a, xi, ypsilon, step_length, rnd);
            double num = rnd.Rannyu();
            if (num < param_space_acceptance(slave, parameters_sampler, beta)) {
                parameters_sampler = slave;
            }
            buffer << parameters_sampler.m << "\t" << parameters_sampler.s
				   << "\t" << parameters_sampler.L << "\t" << parameters_sampler.l << "\n";
        }
        beta *= 2;
    }

	trash = Energy(step_length, parameters_sampler.m, parameters_sampler.s, rnd, false, true);

    std::ofstream out;
    out.open("output.dat");
    out << buffer.str();
    out.close();

	/*for (int i=0; i<100; ++i)
		std::cerr << std::get<0>(Energy(step_length, -1.42806, 0.0145151, rnd, false, false)) << "\t" << std::get<1>(Energy(step_length, -1.42806, 0.0145151, rnd, false, false)) << "\n";
	*/

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
