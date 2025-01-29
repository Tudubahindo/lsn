# lsn
My git repo for the Laboratorio di Simulazione Numerica class of Prof. D. E. Galli, academic year 2023/24, unimi

The course material is divided in 12 lectures, and for each lecture the students are required to complete an exercise. For each lecture I created a specific directory containing the lecture material, the assignment and my solution.

The code has been tested on tolab machines, so it should run without issues there. Jupyter notebooks require tensorflow (for l11 & l12) and geopandas (for l10). There's no need to re-run the Julia scripts as the results are already saved in the relevant output files, but just in case I will list the necessary packages: TravelingSalesmanExact, SCIP, CSV, DataFrames, DelimitedFiles.

Three exercises (l04, l06, l07) require the use of a simulator. The source code for this simulator can be found in the `simulator` directory, together with the scripts `l04.sh`, `l06.sh`, `l07.sh` and `all.sh` to produce all the necessary output files. The usage for the simulator is:

`./simulator.exe {n} {directory_path}`

where `{n}` is a natural number and `{directory_path}` is the path of the input directory. The input directory HAS to include an `input.dat` file, a `properties.dat` file, the `seed.in` and `Primes` files, and all the others necessary input files (depending on the chosen simulator mode). The directory also has to include an `out/` subdirectory where the simulator will write the output files. The `n` parameter is the number of cycles: the simulator will perform the simulation `n` times, each time increasing the starting temperature by `0.1`, and each time with different primes for the RNG initialization. All the cycles write to the same output file, in order.

All the other exercises are contained in a single c++ program named `main.cpp` located in the relevant directory. Use `make` to compile and `./main.exe` to execute (c++17 required). l10 requires MPI to compile and `mpiexec` to run, so I suggest to just run the `all.sh` script instead.

Some programs may take a LONG time to run, so I uploaded the output files as well. It should suffice to just run the jupyter notebooks.

Special thanks to my colleagues Alessio Serraino and Stefano Pilosio for the invaluable feedback they provided on some technical issues.
