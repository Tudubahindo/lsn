time ./simulator.exe 21 ../l4/md_solid_calibration/
time ./simulator.exe 21 ../l4/md_liquid_calibration/
time ./simulator.exe 21 ../l4/md_gas_calibration/
cp ../l4/md_solid/input_eq.dat ../l4/md_solid/input.dat
cp ../l4/md_liquid/input_eq.dat ../l4/md_liquid/input.dat
cp ../l4/md_gas/input_eq.dat ../l4/md_gas/input.dat
time ./simulator.exe 1 ../l4/md_solid/
time ./simulator.exe 1 ../l4/md_liquid/
time ./simulator.exe 1 ../l4/md_gas/
cp ../l4/md_solid/input_sim.dat ../l4/md_solid/input.dat
cp ../l4/md_liquid/input_sim.dat ../l4/md_liquid/input.dat
cp ../l4/md_gas/input_sim.dat ../l4/md_gas/input.dat
time ./simulator.exe 1 ../l4/md_solid/
time ./simulator.exe 1 ../l4/md_liquid/
time ./simulator.exe 1 ../l4/md_gas/
time ./simulator.exe 16 ../l6/ising_gibbs_calibration_0/
time ./simulator.exe 16 ../l6/ising_gibbs_calibration_infinity/
time ./simulator.exe 16 ../l6/ising_mrt2_calibration_0/
time ./simulator.exe 16 ../l6/ising_mrt2_calibration_infinity/
time ./simulator.exe 16 ../l6/ising_gibbs/
time ./simulator.exe 16 ../l6/ising_gibbs_m/
time ./simulator.exe 16 ../l6/ising_mrt2/
time ./simulator.exe 16 ../l6/ising_mrt2_m/
time ./simulator.exe 1 ../l7/mc_solid_calibration/
time ./simulator.exe 1 ../l7/mc_liquid_calibration/
time ./simulator.exe 1 ../l7/mc_gas_calibration/
time ./simulator.exe 1 ../l7/mc_solid/
time ./simulator.exe 1 ../l7/mc_liquid/
time ./simulator.exe 1 ../l7/mc_gas/
