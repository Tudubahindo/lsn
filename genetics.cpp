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

std::vector<long unsigned int> map::bruteforce() {
    std::vector<long unsigned int> best;
    std::vector<long unsigned int> slave;
    double shortest_path{};
    for (long unsigned int i = 1; i < cities.size(); ++i) {
        best.push_back(i);
        slave.push_back(i);
        shortest_path += distances(i, i - 1);
    }
    shortest_path += distances(0, cities.size() - 1);

    while (std::next_permutation(slave.begin(), slave.end())) {
        double path{};
        for (long unsigned int i = 1; i < slave.size(); ++i) {
            path += distances(slave.at(i), slave.at(i - 1));
        }
        path += distances(slave.at(0), 0);
        path += distances(slave.size() - 1, 0);

        if (path < shortest_path) {
            shortest_path = path;
            best = slave;
        }
    }

    best.insert(best.begin(), 0);
    for (long unsigned int i = 0; i < best.size(); ++i) {
        ++best.at(i);
    }
    return best;
}

bool chromosome::check() {
    bool flag = true;
    if (_genes.at(0) != 1)
        flag = false;
    if (_genes.size() != _c->cities.size())
        flag = false;
    for (long unsigned int i = 1; i <= _genes.size(); ++i) {
        bool found = false;
        for (long unsigned int j = 0; j < _genes.size(); ++j) {
            if (_genes.at(j) == i) {
                found = true;
            }
        }
        if (found == false)
            flag = false;
    }
    return flag;
}

void chromosome::swap(int a, int b) {
    int slave = _genes.at(a);
    _genes.at(a) = _genes.at(b);
    _genes.at(b) = slave;
}

void chromosome::scramble(Random &rnd) {
    int N = _genes.size() - 2;
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
    for (long unsigned int j = 1; j < _genes.size() - 1; ++j) {
        if (_c->distances(_genes.back() - 1, _genes.at(0) - 1) + _c->distances(_genes.at(j) - 1, _genes.at(j + 1) - 1) > _c->distances(_genes.back() - 1, _genes.at(j) - 1) + _c->distances(_genes.at(0) - 1, _genes.at(j + 1) - 1)) {
            flag = true;
        }
    }

    for (long unsigned int i = 0; i < _genes.size() - 3; ++i) {
        for (long unsigned int j = i + 2; j < _genes.size() - 1; ++j) {
            if (_c->distances(_genes.at(i) - 1, _genes.at(i + 1) - 1) + _c->distances(_genes.at(j) - 1, _genes.at(j + 1) - 1) > _c->distances(_genes.at(i) - 1, _genes.at(j) - 1) + _c->distances(_genes.at(i + 1) - 1, _genes.at(j + 1) - 1)) {
                flag = true;
            }
        }
    }

    return flag;
}

void chromosome::untangle(Random &rnd) {
    long unsigned int start = static_cast<int>(std::floor(_genes.size() * rnd.Rannyu())); // start untagling from a random point
    for (long unsigned int j = 0; j < start; ++j) {
        int slave = _genes.at(0);
        for (long unsigned int i = 1; i < _genes.size(); ++i) {
            _genes.at(i - 1) = _genes.at(i);
        }
        _genes.back() = slave;
    }

    bool flag = true;
    while (flag == true) {
        flag = false;
        for (long unsigned int j = 1; j < _genes.size() - 1; ++j) {
            if (_c->distances(_genes.back() - 1, _genes.at(0) - 1) + _c->distances(_genes.at(j) - 1, _genes.at(j + 1) - 1) > _c->distances(_genes.back() - 1, _genes.at(j) - 1) + _c->distances(_genes.at(0) - 1, _genes.at(j + 1) - 1)) {
                flag = true;
                for (long unsigned int k = 0; k < 0.5 * (j + 1); ++k) {
                    swap(k, j - k);
                }
            }
        }

        for (long unsigned int i = 0; i < _genes.size() - 3; ++i) {
            for (long unsigned int j = i + 2; j < _genes.size() - 1; ++j) {
                if (_c->distances(_genes.at(i) - 1, _genes.at(i + 1) - 1) + _c->distances(_genes.at(j) - 1, _genes.at(j + 1) - 1) > _c->distances(_genes.at(i) - 1, _genes.at(j) - 1) + _c->distances(_genes.at(i + 1) - 1, _genes.at(j + 1) - 1)) {
                    flag = true;
                    for (long unsigned int k = 0; k < 0.5 * (j - i); ++k) {
                        swap(i + 1 + k, j - k);
                    }
                }
            }
        }
    }

    while (_genes.at(0) != 1) {
        int slave = _genes.at(0);
        for (long unsigned int i = 1; i < _genes.size(); ++i) {
            _genes.at(i - 1) = _genes.at(i);
        }
        _genes.back() = slave;
    }
    _fitness = std::pow(1.0 / L(), 4);
}

chromosome::chromosome(map m, bool test, Random &rnd) {
    int N = m.cities.size();
    _c = std::make_shared<map>(m);
    for (int i = 1; i <= N; ++i) {
        _genes.push_back(i);
    }
    scramble(rnd);
    if (test == false) {
        untangle(rnd);
    }
    _fitness = std::pow(1.0 / L(), 4);
}

chromosome::chromosome(map m, std::vector<long unsigned int> in) {
    _genes = in;
    _c = std::make_shared<map>(m);
    _fitness = std::pow(1.0 / L(), 4);
}

chromosome::chromosome(chromosome p, chromosome m, int stop, Random &rnd) {
    /*c = p.c;
    for (int i = 0; i < stop; ++i) {
        _genes.push_back(p._genes.at(i));
    }

    for (long unsigned int i = 0; i < m._genes.size(); ++i) {
        bool found = false;
        for (long unsigned int j = 0; j < _genes.size(); ++j) {
            if (_genes.at(j) == m._genes.at(i)) {
                found = true;
            }
        }
        if (found == false) {
            _genes.push_back(m._genes.at(i));
        }
    }
    _fitness = std::pow(100.0 / L(), 4);*/

    _c = p._c;
    for (int i = 0; i < stop; ++i) {
        _genes.push_back(p._genes.at(i));
    }

    long unsigned int i = 0;
    while (_genes.size() < m._genes.size()) {
        bool found = false;
        for (long unsigned int j = 0; j < _genes.size(); ++j) {
            if (_genes.at(j) == m._genes.at(i)) {
                found = true;
            }
        }
        if (found == false) {
            bool tangle = false;
            for (long unsigned int j = 0; j < _genes.size() - 2; ++j) {
                if (_c->distances(m._genes.at(i) - 1, _genes.back() - 1) + _c->distances(_genes.at(j) - 1, _genes.at(j + 1) - 1) > _c->distances(_genes.at(j + 1) - 1, _genes.back() - 1) + _c->distances(_genes.at(j) - 1, m._genes.at(i) - 1)) {
                    tangle = true;
                }
            }

            if (tangle == false) {
                _genes.push_back(m._genes.at(i));
                i = -1;
            }
        }
        ++i;

        if (i >= m._genes.size()) {
            for (long unsigned int i = 0; i < m._genes.size(); ++i) {
                bool found = false;
                for (long unsigned int j = 0; j < _genes.size(); ++j) {
                    if (_genes.at(j) == m._genes.at(i)) {
                        found = true;
                    }
                }
                if (found == false) {
                    _genes.push_back(m._genes.at(i));
                    i = 0;
                }
            }
        }
    }
    untangle(rnd);
    _fitness = std::pow(1.0 / L(), 4);
}

void chromosome::mutation_random_swap(Random &rnd) {
    int N = _genes.size() - 2;
    if (N > 0) {
        int first = 1 + static_cast<int>(std::floor((N + 1) * rnd.Rannyu()));
        int second = 1 + static_cast<int>(std::floor(N * rnd.Rannyu()));
        if (second == first)
            second = N + 1;
        swap(first, second);
    }
    _fitness = std::pow(1.0 / L(), 4);
}

void chromosome::mutation_shift(Random &rnd) {
    int N = _genes.size() - 2;
    if (N > 0) {
        int num = 1 + static_cast<int>(std::floor((N / 2) * rnd.Rannyu()));     // dimension of the subarray to switch forward (between 1 and N)
        int steps = 1 + static_cast<int>(std::floor((N / 2) * rnd.Rannyu()));   // number of steps forward (between 1 and N)
        int start = 1 + static_cast<int>(std::floor((N - num) * rnd.Rannyu())); // index of the start of said subarray (between 1 and N-num)
        std::vector<long unsigned int> _genes_extended = _genes;
        for (long unsigned int i = 1; i < _genes.size(); ++i) {
            _genes_extended.push_back(_genes.at(i));
        }
        for (int k = 0; k < steps; ++k) { // k max is steps-1 (k max is N-1)
            int j = start + k + num;
            int slave = _genes_extended.at(j);
            while (j > start + k) {
                _genes_extended.at(j) = _genes_extended.at(j - 1);
                --j;
            };
            _genes_extended.at(start + k) = slave;
        }
        for (long unsigned int i = 1; i < _genes.size(); ++i) {
            _genes.at(i) = _genes_extended.at(i);
        }

        for (int i = 0; i < start - 1; ++i) {
            _genes.at(i + 1) = _genes_extended.at(_genes.size() + i);
        }
    }
    _fitness = std::pow(1.0 / L(), 4);
}

void chromosome::mutation_permutation(Random &rnd) {
    int N = _genes.size() - 2;
    int Nhalf = N / 2;
    if (Nhalf > 0) {
        int num = 1 + static_cast<int>(std::floor(Nhalf * rnd.Rannyu()));                     // dimension of the 2 subarrays to permute
        int space = static_cast<int>(std::floor((N - 2 * num) * rnd.Rannyu()));               // space between the 2 subarrays
        int start = 1 + static_cast<int>(std::floor((N - (2 * num + space)) * rnd.Rannyu())); // index of the first element of the first subarray
        for (int i = 0; i < num; ++i) {
            swap(start + i, start + num + space + i);
        }
    }
    _fitness = std::pow(1.0 / L(), 4);
}

void chromosome::mutation_inversion(Random &rnd) {
    int N = _genes.size() - 2;
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
    _fitness = std::pow(1.0 / L(), 4);
}

double chromosome::L() {
    double length{};
    for (long unsigned int i = 1; i < _genes.size(); ++i) {
        length += _c->distances(_genes.at(i - 1) - 1, _genes.at(i) - 1);
    }
    length += _c->distances(_genes.back() - 1, _genes.at(0) - 1);
    return length;
}

double chromosome::fit_getter() const {
    return _fitness;
}

std::string chromosome::print() const {
    std::stringstream buffer;
    for (auto i : _genes) {
        buffer << i;
        if (i != _genes.back()) {
            buffer << "\t";
        }
    }
    return buffer.str();
}

std::vector<long unsigned int> chromosome::intprint() const {
    return _genes;
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
