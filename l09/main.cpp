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
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <vector>
#define _USE_MATH_DEFINES

typedef std::tuple<double, double> city;

double d2(city a, city b) {
    return std::pow(std::get<0>(a) - std::get<0>(b), 2) + std::pow(std::get<1>(a) - std::get<1>(b), 2);
}

double d1(city a, city b) {
    return std::sqrt(std::pow(std::get<0>(a) - std::get<0>(b), 2) + std::pow(std::get<1>(a) - std::get<1>(b), 2));
}

struct map {
    std::vector<city> cities;
    arma::mat distances;

    map(std::vector<city> c) {
        cities = c;
        arma::mat slave(c.size(), c.size());
        for (long unsigned int i = 0; i < c.size(); ++i) {
            for (long unsigned int j = 0; j < c.size(); ++j) {
                slave(i, j) = d1(c.at(i), c.at(j));
            }
        }
        distances = slave;
    }
};

struct chromosome {

    std::vector<int> genes;
    std::shared_ptr<map> c;
    double fitness;

    bool check() {
        bool flag = true;
        if (genes.at(0) != 1)
            flag = false;
        if (genes.size() != c->cities.size())
            flag = false;
        for (long unsigned int i = 1; i <= genes.size(); ++i) {
            bool found = false;
            for (long unsigned int j = 0; j < genes.size(); ++j) {
                if (genes.at(j) == i) {
                    found = true;
                }
            }
            if (found == false)
                flag = false;
        }
        return flag;
    }

    void swap(int a, int b) {
        int slave = genes.at(a);
        genes.at(a) = genes.at(b);
        genes.at(b) = slave;
    }

    void scramble(Random &rnd) {
        int N = genes.size() - 2;
        if (N > 0) {
            for (int i = 0; i < 10 * N; ++i) {
                int first = 1 + static_cast<int>(std::floor((N + 1) * rnd.Rannyu()));
                int second = 1 + static_cast<int>(std::floor(N * rnd.Rannyu()));
                if (second == first)
                    second = N + 1;
                swap(first, second);
            }
        }
    }

    chromosome(map m, Random &rnd) {
        int N = m.cities.size();
        c = std::make_shared<map>(m);
        for (int i = 1; i <= N; ++i) {
            genes.push_back(i);
        }
        scramble(rnd);
        fitness = std::pow(100.0 / L(), 4);
    }

    chromosome(chromosome p, chromosome m, int stop) {
        c = p.c;
        for (int i = 0; i < stop; ++i) {
            genes.push_back(p.genes.at(i));
        }
        for (long unsigned int i = 0; i < m.genes.size(); ++i) {
            bool found = false;
            for (long unsigned int j = 0; j < genes.size(); ++j) {
                if (genes.at(j) == m.genes.at(i)) {
                    found = true;
                }
            }
            if (found == false) {
                genes.push_back(m.genes.at(i));
            }
        }
        fitness = std::pow(100.0 / L(), 4);
    }

    void mutation_random_swap(Random &rnd) {
        int N = genes.size() - 2;
        if (N > 0) {
            int first = 1 + static_cast<int>(std::floor((N + 1) * rnd.Rannyu()));
            int second = 1 + static_cast<int>(std::floor(N * rnd.Rannyu()));
            if (second == first)
                second = N + 1;
            swap(first, second);
        }
        fitness = std::pow(100.0 / L(), 4);
    }

    void mutation_shift(Random &rnd) {
        int N = genes.size() - 2;
        if (N > 0) {
            int num = 1 + static_cast<int>(std::floor((N / 2) * rnd.Rannyu()));     // dimension of the subarray to switch forward (between 1 and N)
            int steps = 1 + static_cast<int>(std::floor((N / 2) * rnd.Rannyu()));   // number of steps forward (between 1 and N)
            int start = 1 + static_cast<int>(std::floor((N - num) * rnd.Rannyu())); // index of the start of said subarray (between 1 and N-num)
            std::vector<int> genes_extended = genes;
            for (long unsigned int i = 1; i < genes.size(); ++i) {
                genes_extended.push_back(genes.at(i));
            }
            for (int k = 0; k < steps; ++k) { // k max is steps-1 (k max is N-1)
                int j = start + k + num;
                int slave = genes_extended.at(j);
                while (j > start + k) {
                    genes_extended.at(j) = genes_extended.at(j - 1);
                    --j;
                };
                genes_extended.at(start + k) = slave;
            }
            for (long unsigned int i = 1; i < genes.size(); ++i) {
                genes.at(i) = genes_extended.at(i);
            }

            for (int i = 0; i < start - 1; ++i) {
                genes.at(i + 1) = genes_extended.at(genes.size() + i);
            }
        }
        fitness = std::pow(100.0 / L(), 4);
    }

    void mutation_permutation(Random &rnd) {
        int N = genes.size() - 2;
        int Nhalf = N / 2;
        if (Nhalf > 0) {
            int num = 1 + static_cast<int>(std::floor(Nhalf * rnd.Rannyu()));                     // dimension of the 2 subarrays to permute
            int space = static_cast<int>(std::floor((N - 2 * num) * rnd.Rannyu()));               // space between the 2 subarrays
            int start = 1 + static_cast<int>(std::floor((N - (2 * num + space)) * rnd.Rannyu())); // index of the first element of the first subarray
            for (int i = 0; i < num; ++i) {
                swap(start + i, start + num + space + i);
            }
        }
        fitness = std::pow(100.0 / L(), 4);
    }

    void mutation_inversion(Random &rnd) {
        int N = genes.size() - 2;
        if (N > 0) {
            int num = 1 + static_cast<int>(std::floor(N * rnd.Rannyu()));           // dimension of the subarray to invert
            int start = 1 + static_cast<int>(std::floor((N - num) * rnd.Rannyu())); // index of the start of said subarray
            int end = start + num;                                                  // index of the end of said subarray
            while (end > start) {
                swap(start, end);
                ++start;
                --end;
            };
        }
        fitness = std::pow(100.0 / L(), 4);
    }

    double L() {
        double length{};
        for (long unsigned int i = 1; i < genes.size(); ++i) {
            length += c->distances(genes.at(i - 1) - 1, genes.at(i) - 1);
        }
        return length;
    }
};

int selection(std::vector<chromosome> &population, double num) {

    double fit_tot{};
    for (long unsigned int i = 0; i < population.size(); ++i) {
        fit_tot += population.at(i).fitness;
    }

    num *= fit_tot;
    int count = -1;
    while (num > 0) {
        ++count;
        num -= population.at(count).fitness;
    };

    return count;
}

//---------------------------------------------------------------------------------------------------------------------------------------

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

    static constexpr int Ncities = 34;
    static constexpr double radius = 1.0;
    static constexpr int popsize = 100;
    static constexpr double Pmutation = 0.02; // 0.98^4 ~ 0.92
    static constexpr double Pcrossover = 0.6;
    std::vector<city> circle;

    for (int i = 0; i < Ncities; ++i) {
        double theta = 2 * M_PI * rnd.Rannyu();
        double x = radius * std::cos(theta);
        double y = radius * std::sin(theta);
        city newcity = {x, y};
        circle.push_back(newcity);
    }

    map circle_map(circle);

    std::vector<chromosome> trial;
    for (int i = 0; i < popsize; ++i) {
        chromosome slave = chromosome(circle, rnd);
        trial.push_back(slave);
    }

    bool ok = true;
    for (int generation = 0; generation < 100; ++generation) {

        std::vector<chromosome> newgen;

        for (int i = 0; i < popsize / 2; ++i) { // crossover
            int pater = selection(trial, rnd.Rannyu());
            int mater = selection(trial, rnd.Rannyu());
            double cross = rnd.Rannyu();
            chromosome son = trial.at(pater);
            chromosome daughter = trial.at(mater);

            if (cross < Pcrossover) {
                son = chromosome(trial.at(pater), trial.at(mater), Ncities / 2);
                daughter = chromosome(trial.at(mater), trial.at(pater), Ncities / 2);
            }
            newgen.push_back(son);
            newgen.push_back(daughter);
        }

        while (newgen.size() < trial.size()) {
            int index = selection(trial, rnd.Rannyu());
            newgen.push_back(trial.at(index));
        }

        trial = newgen;

        for (int i = 0; i < popsize; ++i) { // mutations
            double num_swap = rnd.Rannyu();
            double num_shift = rnd.Rannyu();
            double num_perm = rnd.Rannyu();
            double num_inv = rnd.Rannyu();

            if (num_swap < Pmutation) {
                trial.at(i).mutation_random_swap(rnd);
            }
            if (num_shift < Pmutation) {
                trial.at(i).mutation_shift(rnd);
            }
            if (num_perm < Pmutation) {
                trial.at(i).mutation_permutation(rnd);
            }
            if (num_inv < Pmutation) {
                trial.at(i).mutation_inversion(rnd);
            }

            if (trial.at(i).check() == false) {
                ok = false;
            }
        }
        std::cout << generation << "\t" << trial.at(0).check() << "\t" << trial.at(0).L() << "\t" << trial.at(0).fitness << "\n";
    }
    std::cerr << "\n"
              << ok << "\n";

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
