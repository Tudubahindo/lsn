time ./simulator.exe 21 ../l04/md_solid_calibration/
time ./simulator.exe 21 ../l04/md_liquid_calibration/
time ./simulator.exe 21 ../l04/md_gas_calibration/
cp ../l04/md_solid/input_eq.dat ../l04/md_solid/input.dat
cp ../l04/md_liquid/input_eq.dat ../l04/md_liquid/input.dat
cp ../l04/md_gas/input_eq.dat ../l04/md_gas/input.dat
time ./simulator.exe 1 ../l04/md_solid/
time ./simulator.exe 1 ../l04/md_liquid/
time ./simulator.exe 1 ../l04/md_gas/
cp ../l04/md_solid/input_sim.dat ../l04/md_solid/input.dat
cp ../l04/md_liquid/input_sim.dat ../l04/md_liquid/input.dat
cp ../l04/md_gas/input_sim.dat ../l04/md_gas/input.dat
time ./simulator.exe 1 ../l04/md_solid/
time ./simulator.exe 1 ../l04/md_liquid/
time ./simulator.exe 1 ../l04/md_gas/
time ./simulator.exe 16 ../l06/ising_gibbs_calibration_0/
time ./simulator.exe 16 ../l06/ising_gibbs_calibration_infinity/
time ./simulator.exe 16 ../l06/ising_mrt2_calibration_0/
time ./simulator.exe 16 ../l06/ising_mrt2_calibration_infinity/
time ./simulator.exe 16 ../l06/ising_gibbs/
time ./simulator.exe 16 ../l06/ising_gibbs_m/
time ./simulator.exe 16 ../l06/ising_mrt2/
time ./simulator.exe 16 ../l06/ising_mrt2_m/
time ./simulator.exe 1 ../l07/mc_solid_calibration/
time ./simulator.exe 1 ../l07/mc_liquid_calibration/
time ./simulator.exe 1 ../l07/mc_gas_calibration/
time ./simulator.exe 1 ../l07/mc_solid/
time ./simulator.exe 1 ../l07/mc_liquid/
time ./simulator.exe 1 ../l07/mc_gas/
