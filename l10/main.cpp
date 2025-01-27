/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "../genetics.h"
#include "../random.h"
#include "../stats.h"
#include "mpi.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#define _USE_MATH_DEFINES

int main(int argc, char *argv[]) {

    // first thing first, let's initialize MPI
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // read input file
    std::string filename;
	std::string option;
	bool flag = false;
    if (argc < 2) {
        std::cout << "Program usage: [-x] <file_name>\nType in the name of the file you want to use: ";
        std::cin >> filename;
        std::cout << "\n";
    } else {
		if (argc > 2) {
			option = argv[1];
			if (option[0] == '-') {
				filename = argv[2];
				flag = true;
			} else {
				filename = option;
			}
		} else {
			filename = argv[1];
		}
    }

    std::ifstream in(filename);

    if (!in) { // checks for errors while opening file
        std::error_code err_code(errno, std::system_category());
        std::cerr << "Error opening \"" << filename << "\" with error: "
                  << err_code.message() << std::endl;
        return 2;
    }

    std::vector<city> cities;
    while (in.good()) {
        double x, y;
        in >> x >> y;
        cities.push_back({x, y});
    }
    in.close();
    cities.pop_back(); // weird error

    Random rnd;
    int seed[4];
    int p1, p2;
    std::ifstream Primes("Primes");
    if (Primes.is_open()) {
		for (int i=0; i<rank; ++i) {
			int trash;
			Primes >> trash >> trash;
		}
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
                seed[3] += 10 + rank; //I tried some seeds and this is the best
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else {
        std::cerr << "PROBLEM: Unable to open seed.in" << std::endl;
    }

    //---------------------------------------------------------------------------------------------------

    static constexpr double Pmutation = 0.1;
    static constexpr double Pcrossover = 0.4;
    static constexpr int popsize = 500;
    static constexpr int cross_size = 50;
    static constexpr int total = 1000;
    std::stringstream buffer;

    map city_map(cities);

    std::vector<chromosome> people;
    for (int i = 0; i < popsize; ++i) {
        chromosome slave = chromosome(city_map, 0, rnd);
        people.push_back(slave);
    }

	std::sort(people.rbegin(), people.rend());
	chromosome best = people.at(0);

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
                son = chromosome(people.at(father), people.at(mother), city_map.cities.size() / 2, rnd);
                daughter = chromosome(people.at(mother), people.at(father), city_map.cities.size() / 2, rnd);
            }
            son.untangle(rnd);
            daughter.untangle(rnd);

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
            //double num_inv = rnd.Rannyu();

            if (num_swap < Pmutation) {
                people.at(i).mutation_random_swap(rnd);
            }
            if (num_shift < Pmutation) {
                people.at(i).mutation_shift(rnd);
            }
            if (num_perm < Pmutation) {
                people.at(i).mutation_permutation(rnd);
            }
            /*if (num_inv < Pmutation) {
                people.at(i).mutation_inversion(rnd);
            }*/

            people.at(i).untangle(rnd);

            if (people.at(i).check() == false) {
                ok = false;
            }
        }

        std::sort(people.rbegin(), people.rend());
        double avg_L{};
        int half = popsize / 2;
        for (int i = 0; i < half; ++i) {
            avg_L += people.at(i).L();
        }
        avg_L /= half;

		if (people.at(0).L() < best.L()) {
			best = people.at(0);
		}

        if (rank == 0)
            std::cout << generation << "\t" << people.at(0).check() << "\t" << people.at(0).L() << "\t" << avg_L << "\t" << people.at(0).fit_getter() << "\t"
                      << "\n";

        /*if (generation % 10 == 9 && flag==true) {
            long unsigned int *receiver = new long unsigned int[size * cities.size()];
            std::vector<long unsigned int> sender = people.at(0).intprint();
            std::vector<long unsigned int> inserter{};
            MPI_Allgather(sender.data(), cities.size(), MPI_UNSIGNED_LONG, receiver, cities.size(), MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
            int num = static_cast<int>(std::floor(size * rnd.Rannyu()));
            for (long unsigned int j = 0; j < cities.size(); ++j) {
                inserter.push_back(receiver[num * cities.size() + j]);
            }
            chromosome slave = chromosome(city_map, inserter);
            people.back() = slave;
            delete[] receiver;
        }*/

		if (generation % 10 == 9 && flag==true) {
            long unsigned int *receiver = new long unsigned int[size * cross_size * cities.size()];
            std::vector<long unsigned int> sender = people.at(0).intprint();
			for (int k = 1; k < cross_size; ++k) {
				std::vector<long unsigned int> slave = people.at(k).intprint();
				sender.insert(sender.end(), slave.begin(), slave.end());
			}
            MPI_Allgather(sender.data(), cross_size * cities.size(), MPI_UNSIGNED_LONG, receiver, cross_size * cities.size(), MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
            //int num = static_cast<int>(std::floor(size * rnd.Rannyu()));
			//num = rank+1;
			//if (num == size) num = 0;
			std::vector<chromosome> slave{};
			for (int k = 0; k < cross_size; ++k) {
				int number = static_cast<int>(std::floor(size * cross_size * rnd.Rannyu()));
				std::vector<long unsigned int> inserter{};
	            for (long unsigned int j = 0; j < cities.size(); ++j) {
					//inserter.push_back(receiver[num * cities.size() * cross_size + k * cities.size() + j]);
					inserter.push_back(receiver[number * cities.size() + j]);
				}
				slave.push_back(chromosome(city_map, inserter));
			}

            people.erase(people.end()-cross_size, people.end());
            //people.erase(people.begin(), people.begin() + cross_size);
			people.insert(people.end(), slave.begin(), slave.end());
            delete[] receiver;
        }
		//std::cerr << people.size() << "\n";
    }

    long unsigned int *receiver = new long unsigned int[size * cities.size()];
    std::vector<long unsigned int> sender = best.intprint();

    MPI_Gather(sender.data(), cities.size(), MPI_UNSIGNED_LONG, receiver, cities.size(), MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        for (int i = 0; i < size; ++i) {
            for (long unsigned int j = 0; j < cities.size(); ++j) {
                buffer << receiver[i * cities.size() + j] << "\t";
            }
            buffer << "\n";
        }
    }

	std::string outname = "output";
	if (flag) outname += "_x";
	outname += ".dat";
    std::ofstream out;
    out.open(outname);
    out << buffer.str();
    out.close();

    delete[] receiver;

    MPI_Finalize();

    rnd.SaveSeed();
    return (!ok);
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
