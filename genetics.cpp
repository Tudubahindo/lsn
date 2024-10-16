/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "genetics.h"
#include "stats.h"

double d2(city a, city b) {
    return std::pow(std::get<0>(a) - std::get<0>(b), 2) + std::pow(std::get<1>(a) - std::get<1>(b), 2);
}

double d1(city a, city b) {
    return std::sqrt(std::pow(std::get<0>(a) - std::get<0>(b), 2) + std::pow(std::get<1>(a) - std::get<1>(b), 2));
}

map::map(std::vector<city> c) : cities(c), distances(SquareMatrix<double>(c.size())) {
    for (long unsigned int i = 0; i < c.size(); ++i) {
        for (long unsigned int j = 0; j < c.size(); ++j) {
            distances(i, j) = d1(c.at(i), c.at(j));
        }
    }
}

bool chromosome::check() {
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

void chromosome::swap(int a, int b) {
    int slave = genes.at(a);
    genes.at(a) = genes.at(b);
    genes.at(b) = slave;
}

void chromosome::scramble(Random &rnd) {
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

bool chromosome::is_tangled() {
    bool flag = false;
    for (long unsigned int j = 1; j < genes.size() - 1; ++j) {
        if (c->distances(genes.back() - 1, genes.at(0) - 1) + c->distances(genes.at(j) - 1, genes.at(j + 1) - 1) > c->distances(genes.back() - 1, genes.at(j) - 1) + c->distances(genes.at(0) - 1, genes.at(j + 1) - 1)) {
            flag = true;
        }
    }

    for (long unsigned int i = 0; i < genes.size() - 3; ++i) {
        for (long unsigned int j = i + 2; j < genes.size() - 1; ++j) {
            if (c->distances(genes.at(i) - 1, genes.at(i + 1) - 1) + c->distances(genes.at(j) - 1, genes.at(j + 1) - 1) > c->distances(genes.at(i) - 1, genes.at(j) - 1) + c->distances(genes.at(i + 1) - 1, genes.at(j + 1) - 1)) {
                flag = true;
            }
        }
    }

    return flag;
}

void chromosome::untangle() {
    bool flag = true;
    while (flag == true) {
        flag = false;
        for (long unsigned int j = 1; j < genes.size() - 1; ++j) {
            if (c->distances(genes.back() - 1, genes.at(0) - 1) + c->distances(genes.at(j) - 1, genes.at(j + 1) - 1) > c->distances(genes.back() - 1, genes.at(j) - 1) + c->distances(genes.at(0) - 1, genes.at(j + 1) - 1)) {
                flag = true;
                for (long unsigned int k = 0; k < 0.5 * (j + 1); ++k) {
                    swap(k, j - k);
                }
            }
        }

        for (long unsigned int i = 0; i < genes.size() - 3; ++i) {
            for (long unsigned int j = i + 2; j < genes.size() - 1; ++j) {
                if (c->distances(genes.at(i) - 1, genes.at(i + 1) - 1) + c->distances(genes.at(j) - 1, genes.at(j + 1) - 1) > c->distances(genes.at(i) - 1, genes.at(j) - 1) + c->distances(genes.at(i + 1) - 1, genes.at(j + 1) - 1)) {
                    flag = true;
                    for (long unsigned int k = 0; k < 0.5 * (j - i); ++k) {
                        swap(i + 1 + k, j - k);
                    }
                }
            }
        }
    }

    while (genes.at(0) != 1) {
        int slave = genes.at(0);
        for (long unsigned int i = 1; i < genes.size(); ++i) {
            genes.at(i - 1) = genes.at(i);
        }
        genes.back() = slave;
    }
    fitness = std::pow(100.0 / L(), 4);
}

chromosome::chromosome(map m, Random &rnd) {
    int N = m.cities.size();
    c = std::make_shared<map>(m);
    for (int i = 1; i <= N; ++i) {
        genes.push_back(i);
    }
    scramble(rnd);
    untangle();
    fitness = std::pow(100.0 / L(), 4);
}

chromosome::chromosome(map m, std::vector<long unsigned int> in) {
    genes = in;
    c = std::make_shared<map>(m);
    fitness = std::pow(100.0 / L(), 4);
}

chromosome::chromosome(chromosome p, chromosome m, int stop) {
    /*c = p.c;
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
    fitness = std::pow(100.0 / L(), 4);*/

    c = p.c;
    for (int i = 0; i < stop; ++i) {
        genes.push_back(p.genes.at(i));
    }

    long unsigned int i = 0;
    while (genes.size() < m.genes.size()) {
        bool found = false;
        for (long unsigned int j = 0; j < genes.size(); ++j) {
            if (genes.at(j) == m.genes.at(i)) {
                found = true;
            }
        }
        if (found == false) {
            bool tangle = false;
            for (long unsigned int j = 0; j < genes.size() - 2; ++j) {
                if (c->distances(m.genes.at(i) - 1, genes.back() - 1) + c->distances(genes.at(j) - 1, genes.at(j + 1) - 1) > c->distances(genes.at(j + 1) - 1, genes.back() - 1) + c->distances(genes.at(j) - 1, m.genes.at(i) - 1)) {
                    tangle = true;
                }
            }

            if (tangle == false) {
                genes.push_back(m.genes.at(i));
                i = -1;
            }
        }
        ++i;

        if (i >= m.genes.size()) {
            for (long unsigned int i = 0; i < m.genes.size(); ++i) {
                bool found = false;
                for (long unsigned int j = 0; j < genes.size(); ++j) {
                    if (genes.at(j) == m.genes.at(i)) {
                        found = true;
                    }
                }
                if (found == false) {
                    genes.push_back(m.genes.at(i));
                    i = 0;
                }
            }
        }
    }
    untangle();
    fitness = std::pow(100.0 / L(), 4);
}

void chromosome::mutation_random_swap(Random &rnd) {
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

void chromosome::mutation_shift(Random &rnd) {
    int N = genes.size() - 2;
    if (N > 0) {
        int num = 1 + static_cast<int>(std::floor((N / 2) * rnd.Rannyu()));     // dimension of the subarray to switch forward (between 1 and N)
        int steps = 1 + static_cast<int>(std::floor((N / 2) * rnd.Rannyu()));   // number of steps forward (between 1 and N)
        int start = 1 + static_cast<int>(std::floor((N - num) * rnd.Rannyu())); // index of the start of said subarray (between 1 and N-num)
        std::vector<long unsigned int> genes_extended = genes;
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

void chromosome::mutation_permutation(Random &rnd) {
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

void chromosome::mutation_inversion(Random &rnd) {
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

double chromosome::L() {
    double length{};
    for (long unsigned int i = 1; i < genes.size(); ++i) {
        length += c->distances(genes.at(i - 1) - 1, genes.at(i) - 1);
    }
    length += c->distances(genes.back() - 1, genes.at(0) - 1);
    return length;
}

double chromosome::fit_getter() const {
    return fitness;
}

std::string chromosome::print() const {
    std::stringstream buffer;
    for (auto i : genes) {
        buffer << i;
        if (i != genes.back()) {
            buffer << "\t";
        }
    }
    return buffer.str();
}

std::vector<long unsigned int> chromosome::intprint() const {
    return genes;
}

int selection(std::vector<chromosome> &population, double num) {

    double fit_tot{};
    for (long unsigned int i = 0; i < population.size(); ++i) {
        fit_tot += population.at(i).fit_getter();
    }

    num *= fit_tot;
    int count = -1;
    while (num > 0) {
        ++count;
        num -= population.at(count).fit_getter();
    };

    return count;
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
