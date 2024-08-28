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

// namespace Sim
// {

void Random::Initializer(const std::string &seed_file,
                         const std::string &prime_file,
                         const std::size_t rows_to_skip) {
    int seed[4];
    int p1, p2;
    std::ifstream Primes(prime_file);
    if (Primes.is_open()) {
        for (std::size_t i = 0; i < rows_to_skip && !Primes.eof(); ++i) {
            Primes.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        if (!Primes.eof()) {
            Primes >> p1 >> p2;
        } else {
            std::cerr << "ERROR: Failed to go to line" << rows_to_skip
                      << " Aborting\n";
            exit(-1);
        }
    } else
        std::cerr << "PROBLEM: Unable to open Primes" << std::endl;
    Primes.close();

    std::ifstream input(seed_file);
    std::string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else
        std::cerr << "PROBLEM: Unable to open seed.in" << std::endl;
}

void Random ::SaveSeed(const std::string &filename) const {
    std::ofstream WriteSeed;
    WriteSeed.open(filename);
    if (WriteSeed.is_open()) {
        WriteSeed << ((l_tot >> 36) & 4095) << " " << ((l_tot >> 24) & 4095)
                  << " " << ((l_tot >> 12) & 4095) << " " << (l_tot & 4095)
                  << std::endl;
    } else
        std::cerr << "PROBLEM: Unable to open " << filename << std::endl;
    WriteSeed.close();
    return;
}

double Random ::Gauss(const double mean, const double sigma) {
    double s = Rannyu();
    double t = Rannyu();
    double x = sqrt(-2. * std::log(1. - s)) * std::cos(2. * M_PI * t);
    return mean + x * sigma;
}

double Random ::Rannyu(const double min, const double max) {
    assert(max > min && "Max should be greather than min\n");
    return min + (max - min) * Rannyu();
}

double Random ::Rannyu(void) {
    constexpr double twom48 = 1. / (1ull << 48);
    l_tot = l_tot * m_tot + n_tot;
    l_tot &= ((1ull << 48) - 1);
    double r = twom48 * l_tot;
    return r;
}

void Random ::SetRandom(int *s, int p1, int p2) {
    l_tot = (static_cast<uint64_t>(s[0]) << (12 * 3)) +
            (static_cast<uint64_t>(s[1]) << (12 * 2)) +
            (static_cast<uint64_t>(s[2]) << 12) + static_cast<uint64_t>(s[3]);

    uint16_t n1 = 0;
    uint16_t n2 = 0;
    uint16_t n3 = p1;
    uint16_t n4 = p2;

    n_tot = (static_cast<uint64_t>(n1) << (12 * 3)) +
            (static_cast<uint64_t>(n2) << (12 * 2)) +
            (static_cast<uint64_t>(n3) << (12 * 1)) + static_cast<uint64_t>(n4);
}
double Random::Exponential(const double lambda) {
    assert(lambda > 0 && "Lambda Parameter should be greater than zero\n");
    double y = Rannyu();
    double r = -std::log(1 - y) / lambda;
    return r;
}

double Random::Lorentzian(const double mu, const double gamma) {
    assert(gamma > 0 && "Gamma should be greater than zero\n");
    double y = Rannyu();
    double r = gamma * std::tan(M_PI * (y - 0.5)) + mu;
    return r;
}

double Random::AcceptReject(const double a, const double b, const double max,
                            std::function<double(double)> &PDF) {
    double x = 0, y = 0;

    do {
        x = Rannyu(a, b);
        y = Rannyu(0, max);
    } while (PDF(x) < y);

    return x;
}

double Random::AcceptReject(const double a, const double b, const double max,
                            const std::function<double(double)> &PDF) {
    double x = 0, y = 0;

    do {
        x = Rannyu(a, b);
        y = Rannyu(0, max);
    } while (PDF(x) < y);

    return x;
}

double Random::ExternalInvCum(std::function<double(double)> &ICDF) {
    return ICDF(Rannyu());
}

uint64_t Random::Ranint() {
    l_tot = l_tot * m_tot + n_tot;
    l_tot &= ((1ull << 48) - 1);
    uint64_t r = l_tot;
    return r;
}

uint64_t Random::Ranint(const uint64_t min, const uint64_t max) {
    assert((max > min) && "supplied wrong value range");
    auto range = max - min + 1;
    auto res = Ranint() % range + min;
    return res;
}

// } // namespace Sim

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
