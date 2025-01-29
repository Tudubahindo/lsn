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
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#define _USE_MATH_DEFINES

struct walker {

    double x;
    double y;
    double z;

    walker(double x_0, double y_0, double z_0) : x{x_0}, y{y_0}, z{z_0} {}

    walker() = default;

    ~walker() = default;

    double r() {
        return std::sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2));
    }

    void step(double a, double xi, double ypsilon, double zeta) {
        x += a * xi;
        y += a * ypsilon;
        z += a * zeta;
    }
};

double psi_100(walker w) { return std::exp(-2 * w.r()); }

double psi_210(walker w) {
    return std::exp(-w.r()) * std::pow(w.z, 2);
}

double Acceptance(double (*psi)(walker), walker v, walker w) {
    return std::min(1.0, psi(v) / psi(w));
}

int main(int argc, char *argv[]) {

    double x_0 = 100.3, y_0 = 73.0, z_0 = -37.901;

    if (argc > 2) {
        x_0 = std::atof(argv[1]);
        y_0 = std::atof(argv[2]);
        z_0 = std::atof(argv[3]);
    }

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

    //----------------------------------exercise 05.1---------------------------------------

    static constexpr int total = 1000000;              // one million throws
    static constexpr int blocksize = 1000;             // one-thousand-throws sized blocks
    static constexpr int blocknum = total / blocksize; // one thousand blocks

    std::stringstream buffer;
    std::stringstream buffer_100;
    std::stringstream buffer_210;
    std::stringstream buffer_autoc;

    walker s_sampler{};
    walker p_sampler{};

    double best_s_acceptance{};
    double best_p_acceptance{};

    double step_length_s = 0.1;
    double step_length_p = 0.1;

    for (int k = 1; k < 51; ++k) {
        s_sampler = walker();
        p_sampler = walker();

        long rate_s{};
        long rate_p{};
        for (int i = 0; i < blocksize; ++i) {
            double xi = 2 * rnd.Rannyu() - 1; // smart idea: use the same rng for both orbitals to halve runtime
            double ypsilon = 2 * rnd.Rannyu() - 1;
            double zeta = 2 * rnd.Rannyu() - 1;
            double num = rnd.Rannyu();

            walker slave = s_sampler;
            slave.step(0.1 * k, xi, ypsilon, zeta);
            if (num < Acceptance(psi_100, slave, s_sampler)) {
                s_sampler = slave;
                rate_s++;
            }

            slave = p_sampler;
            slave.step(0.1 * k, xi, ypsilon, zeta);
            if (num < Acceptance(psi_210, slave, p_sampler)) {
                p_sampler = slave;
                rate_p++;
            }
        }
        if (k == 1) {
            best_s_acceptance = (double)rate_s / (double)blocksize;
            best_p_acceptance = (double)rate_p / (double)blocksize;
        } else {
            if (std::abs((double)rate_s / (double)blocksize - 0.5) < std::abs(best_s_acceptance - 0.5)) {
                best_s_acceptance = (double)rate_s / (double)blocksize;
                step_length_s = 0.1 * k;
            }
            if (std::abs((double)rate_p / (double)blocksize - 0.5) < std::abs(best_p_acceptance - 0.5)) {
                best_p_acceptance = (double)rate_p / (double)blocksize;
                step_length_p = 0.1 * k;
            }
        }
    }

    std::cerr << step_length_s << "\t" << step_length_p << "\n";
    std::cerr << best_s_acceptance << "\t" << best_p_acceptance << "\n";

    s_sampler = walker();
    p_sampler = walker();

    double r_s{};
    double r_s2{};
    double r_p{};
    double r_p2{};

    std::vector<double> s_orbital;
    std::vector<double> p_orbital;

    for (int i = 0; i < blocknum; ++i) {
        double avg_s{};
        double avg_p{};

        for (int j = 0; j < blocksize; ++j) {
            double xi = 2 * rnd.Rannyu() - 1; // smart idea: use the same rng for both orbitals to halve runtime
            double ypsilon = 2 * rnd.Rannyu() - 1;
            double zeta = 2 * rnd.Rannyu() - 1;
            double num = rnd.Rannyu();

            walker slave = s_sampler;
            slave.step(step_length_s, xi, ypsilon, zeta);
            if (num < Acceptance(psi_100, slave, s_sampler)) {
                s_sampler = slave;
            }

            slave = p_sampler;
            slave.step(step_length_p, xi, ypsilon, zeta);
            if (num < Acceptance(psi_210, slave, p_sampler)) {
                p_sampler = slave;
            }

            buffer_100 << s_sampler.x << "\t" << s_sampler.y << "\t" << s_sampler.z << "\t" << s_sampler.r() << "\n";
            buffer_210 << p_sampler.x << "\t" << p_sampler.y << "\t" << p_sampler.z << "\t" << p_sampler.r() << "\n";

            s_orbital.push_back(s_sampler.r());
            p_orbital.push_back(p_sampler.r());

            avg_s += s_sampler.r();
            avg_p += p_sampler.r();
        }

        avg_s /= (double)blocksize;
        avg_p /= (double)blocksize;

        r_s = r_s * ((double)i / (double)(i + 1)) + avg_s / (double)(i + 1);
        r_p = r_p * ((double)i / (double)(i + 1)) + avg_p / (double)(i + 1);

        r_s2 = r_s2 * ((double)i / (double)(i + 1)) +
               std::pow(avg_s, 2) / (double)(i + 1);
        r_p2 = r_p2 * ((double)i / (double)(i + 1)) +
               std::pow(avg_p, 2) / (double)(i + 1);

        double sigma_s = std::sqrt((r_s2 - std::pow(r_s, 2)) / i);
        double sigma_p = std::sqrt((r_p2 - std::pow(r_p, 2)) / i);
        if (i == 0) {
            sigma_s = 0;
            sigma_p = 0;
        }

        buffer << r_s << "\t" << sigma_s << "\t" << r_p << "\t" << sigma_p << "\n";
    }

    for (int i = 0; i < 10 * blocksize; ++i) {
        buffer_autoc << autocorrelation(s_orbital, i) << "\t" << autocorrelation(p_orbital, i) << "\n";
    }

    std::ofstream out;
    out.open("output.dat");
    out << buffer.str();
    out.close();

    out.open("output_100.dat");
    out << buffer_100.str();
    out.close();

    out.open("output_210.dat");
    out << buffer_210.str();
    out.close();

    out.open("output_autoc.dat");
    out << buffer_autoc.str();
    out.close();

    buffer_100.str(std::string());
    buffer_210.str(std::string());

    rnd.SetRandom(seed, p1, p2); // reset rnd

    s_sampler = walker(x_0, y_0, z_0); // now I start far away from the origin (initial coords provided as arguments)
    p_sampler = walker(x_0, y_0, z_0);

    for (int i = 0; i < total; ++i) {
        double xi = 2 * rnd.Rannyu() - 1; // smart idea: use the same rng for both orbitals to halve runtime
        double ypsilon = 2 * rnd.Rannyu() - 1;
        double zeta = 2 * rnd.Rannyu() - 1;
        double num = rnd.Rannyu();

        walker slave = s_sampler;
        slave.step(step_length_s, xi, ypsilon, zeta);
        if (num < Acceptance(psi_100, slave, s_sampler)) {
            s_sampler = slave;
        }

        slave = p_sampler;
        slave.step(step_length_p, xi, ypsilon, zeta);
        if (num < Acceptance(psi_210, slave, p_sampler)) {
            p_sampler = slave;
        }

        buffer_100 << s_sampler.x << "\t" << s_sampler.y << "\t" << s_sampler.z << "\t" << s_sampler.r() << "\n";
        buffer_210 << p_sampler.x << "\t" << p_sampler.y << "\t" << p_sampler.z << "\t" << p_sampler.r() << "\n";
    }

    out.open("output_100_far.dat");
    out << buffer_100.str();
    out.close();

    out.open("output_210_far.dat");
    out << buffer_210.str();
    out.close();

    buffer.str(std::string());
    buffer_100.str(std::string());
    buffer_210.str(std::string());

    rnd.SetRandom(seed, p1, p2); // reset rnd

    best_s_acceptance = 0;
    best_p_acceptance = 0;
    step_length_s = 0;
    step_length_p = 0;

    for (int k = 1; k < 51; ++k) {
        s_sampler = walker();
        p_sampler = walker();

        long rate_s{};
        long rate_p{};
        for (int i = 0; i < blocksize; ++i) {
            double xi = rnd.Gauss(0.0, 1.0); // smart idea: use the same rng for both orbitals to halve runtime
            double ypsilon = rnd.Gauss(0.0, 1.0);
            double zeta = rnd.Gauss(0.0, 1.0);
            double num = rnd.Rannyu();

            walker slave = s_sampler;
            slave.step(0.1 * k, xi, ypsilon, zeta);
            if (num < Acceptance(psi_100, slave, s_sampler)) {
                s_sampler = slave;
                rate_s++;
            }

            slave = p_sampler;
            slave.step(0.1 * k, xi, ypsilon, zeta);
            if (num < Acceptance(psi_210, slave, p_sampler)) {
                p_sampler = slave;
                rate_p++;
            }
        }
        if (k == 1) {
            best_s_acceptance = (double)rate_s / (double)blocksize;
            best_p_acceptance = (double)rate_p / (double)blocksize;
        } else {
            if (std::abs((double)rate_s / (double)blocksize - 0.5) < std::abs(best_s_acceptance - 0.5)) {
                best_s_acceptance = (double)rate_s / (double)blocksize;
                step_length_s = 0.1 * k;
            }
            if (std::abs((double)rate_p / (double)blocksize - 0.5) < std::abs(best_p_acceptance - 0.5)) {
                best_p_acceptance = (double)rate_p / (double)blocksize;
                step_length_p = 0.1 * k;
            }
        }
    }

    std::cerr << step_length_s << "\t" << step_length_p << "\n";
    std::cerr << best_s_acceptance << "\t" << best_p_acceptance << "\n";

    r_s = 0;
    r_s2 = 0;
    r_p = 0;
    r_p2 = 0;

    s_sampler = walker();
    p_sampler = walker();

    for (int i = 0; i < blocknum; ++i) {
        double avg_s{};
        double avg_p{};

        for (int j = 0; j < blocksize; ++j) {
            double xi = rnd.Gauss(0.0, 1.0); // smart idea: use the same rng for both orbitals to halve runtime
            double ypsilon = rnd.Gauss(0.0, 1.0);
            double zeta = rnd.Gauss(0.0, 1.0);
            double num = rnd.Rannyu();

            walker slave = s_sampler;
            slave.step(step_length_s, xi, ypsilon, zeta);
            if (num < Acceptance(psi_100, slave, s_sampler)) {
                s_sampler = slave;
            }

            slave = p_sampler;
            slave.step(step_length_p, xi, ypsilon, zeta);
            if (num < Acceptance(psi_210, slave, p_sampler)) {
                p_sampler = slave;
            }

            buffer_100 << s_sampler.x << "\t" << s_sampler.y << "\t" << s_sampler.z << "\t" << s_sampler.r() << "\n";
            buffer_210 << p_sampler.x << "\t" << p_sampler.y << "\t" << p_sampler.z << "\t" << p_sampler.r() << "\n";

            avg_s += s_sampler.r();
            avg_p += p_sampler.r();
        }

        avg_s /= (double)blocksize;
        avg_p /= (double)blocksize;

        r_s = r_s * ((double)i / (double)(i + 1)) + avg_s / (double)(i + 1);
        r_p = r_p * ((double)i / (double)(i + 1)) + avg_p / (double)(i + 1);

        r_s2 = r_s2 * ((double)i / (double)(i + 1)) +
               std::pow(avg_s, 2) / (double)(i + 1);
        r_p2 = r_p2 * ((double)i / (double)(i + 1)) +
               std::pow(avg_p, 2) / (double)(i + 1);

        double sigma_s = std::sqrt((r_s2 - std::pow(r_s, 2)) / i);
        double sigma_p = std::sqrt((r_p2 - std::pow(r_p, 2)) / i);
        if (i == 0) {
            sigma_s = 0;
            sigma_p = 0;
        }

        buffer << r_s << "\t" << sigma_s << "\t" << r_p << "\t" << sigma_p << "\n";
    }

    out.open("output_gauss.dat");
    out << buffer.str();
    out.close();

    out.open("output_100_gauss.dat");
    out << buffer_100.str();
    out.close();

    out.open("output_210_gauss.dat");
    out << buffer_210.str();
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
