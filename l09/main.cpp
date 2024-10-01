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
#include "utils.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#define _USE_MATH_DEFINES

//---------------------------------------------------------------------------------------------------------------------------------------

bool genetic_algo(int popsize, int total, map in_map, std::string filename, Random &rnd) {
    static constexpr double Pmutation = 0.02; // 0.98^4 ~ 0.92
    static constexpr double Pcrossover = 0.6;
    std::stringstream buffer;

    std::vector<chromosome> people;
    for (int i = 0; i < popsize; ++i) {
        chromosome slave = chromosome(in_map, rnd);
        people.push_back(slave);
    }

    // std::sort(people.rbegin(), people.rend());

    bool ok = true;
    for (int generation = 0; generation < total; ++generation) {

        std::vector<chromosome> newgen;

        for (int i = 0; i < popsize / 2; ++i) { // crossover
            int father = selection(people, rnd.Rannyu());
            int mother = selection(people, rnd.Rannyu());
            double cross = rnd.Rannyu();
            chromosome son = people.at(father);
            chromosome daughter = people.at(mother);

            if (cross < Pcrossover) {
                son = chromosome(people.at(father), people.at(mother), in_map.cities.size() / 2);
                daughter = chromosome(people.at(mother), people.at(father), in_map.cities.size() / 2);
            }
            son.untangle();
            daughter.untangle();

            newgen.push_back(son);
            newgen.push_back(daughter);
        }

        while (newgen.size() < people.size()) {
            int index = selection(people, rnd.Rannyu());
            newgen.push_back(people.at(index));
        }

        people = newgen;

        for (int i = 0; i < popsize; ++i) { // mutations
            double num_swap = rnd.Rannyu();
            double num_shift = rnd.Rannyu();
            double num_perm = rnd.Rannyu();
            double num_inv = rnd.Rannyu();

            if (num_swap < Pmutation) {
                people.at(i).mutation_random_swap(rnd);
            }
            if (num_shift < Pmutation) {
                people.at(i).mutation_shift(rnd);
            }
            if (num_perm < Pmutation) {
                people.at(i).mutation_permutation(rnd);
            }
            if (num_inv < Pmutation) {
                people.at(i).mutation_inversion(rnd);
            }

            people.at(i).untangle();

            if (people.at(i).check() == false) {
                ok = false;
            }
        }

        std::sort(people.rbegin(), people.rend());
		double avg_L{};
		int half = popsize/2;
		for (int i = 0; i < half; ++i){
			avg_L += people.at(i).L();
		}
		avg_L /= half;

        buffer << generation << "\t" << people.at(0).check() << "\t" << people.at(0).L() << "\t" << avg_L << "\t" << people.at(0).fit_getter() << "\t" << "\n";
    }

    std::ofstream out;
    out.open("verbose_"+filename);
    for (auto i : people) {
        out << i.print() << "\n";
    }
    out.close();

	out.open(filename);
    out << buffer.str();
    out.close();
    
	return ok;
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

    //---------------------------------------------------------------------------------------------------

    static constexpr int Ncities = 34;
    static constexpr double radius = 1.0;
    static constexpr int popsize = 100;
    std::vector<city> circle;
    std::vector<city> square;
    std::stringstream buffer;

    for (int i = 0; i < Ncities; ++i) {
        double theta = 2 * M_PI * rnd.Rannyu();
        double x = radius * std::cos(theta);
        double y = radius * std::sin(theta);
        city newcity = {x, y};
        circle.push_back(newcity);
        buffer << x << "\t" << y << "\n";
    }

    std::ofstream out;
    out.open("map_circle.dat");
    out << buffer.str();
    out.close();

    buffer.str(std::string());

    for (int i = 0; i < Ncities; ++i) {
        double x = 2 * rnd.Rannyu() - 1;
        double y = 2 * rnd.Rannyu() - 1;
        city newcity = {x, y};
        square.push_back(newcity);
        buffer << x << "\t" << y << "\n";
    }

    out.open("map_square.dat");
    out << buffer.str();
    out.close();

    map square_map(square);
    map circle_map(circle);

    chromosome test = chromosome(circle_map, rnd);

    out.open("output_test.dat");
    out << test.print() << "\n";
    test.untangle();
    out << test.print();
    out.close();

    std::cerr << "\n"
              << genetic_algo(popsize, 300, square_map, "output_square.dat", rnd) << "\n";

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
