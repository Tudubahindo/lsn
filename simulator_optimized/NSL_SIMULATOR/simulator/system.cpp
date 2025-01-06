/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "system.h"
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#define _USE_MATH_DEFINES 

using namespace std;
using namespace arma;

void System ::step() { // Perform a simulation step
    if (_sim_type == 0)
        this->Verlet(); // Perform a MD step
    else
        for (int i = 0; i < _npart; i++)
            this->move(int(_rnd.Rannyu() * _npart)); // Perform a MC step on a randomly choosen particle
    _nattempts += _npart;                            // update number of attempts performed on the system
    return;
}

void System ::Verlet() {
    double xnew, ynew, znew;
    for (int i = 0; i < _npart; i++) { // Force acting on particle i
        _fx(i) = this->Force(i, 0);
        _fy(i) = this->Force(i, 1);
        _fz(i) = this->Force(i, 2);
    }
    for (int i = 0; i < _npart; i++) { // Verlet integration scheme
        xnew = this->pbc(2.0 * _particle(i).getposition(0, true) - _particle(i).getposition(0, false) + _fx(i) * pow(_delta, 2), 0);
        ynew = this->pbc(2.0 * _particle(i).getposition(1, true) - _particle(i).getposition(1, false) + _fy(i) * pow(_delta, 2), 1);
        znew = this->pbc(2.0 * _particle(i).getposition(2, true) - _particle(i).getposition(2, false) + _fz(i) * pow(_delta, 2), 2);
        _particle(i).setvelocity(0, this->pbc(xnew - _particle(i).getposition(0, false), 0) / (2.0 * _delta));
        _particle(i).setvelocity(1, this->pbc(ynew - _particle(i).getposition(1, false), 1) / (2.0 * _delta));
        _particle(i).setvelocity(2, this->pbc(znew - _particle(i).getposition(2, false), 2) / (2.0 * _delta));
        _particle(i).acceptmove(); // xold = xnew
        _particle(i).setposition(0, xnew);
        _particle(i).setposition(1, ynew);
        _particle(i).setposition(2, znew);
    }
    _naccepted += _npart;
    return;
}

double System ::Force(int i, int dim) {
    double f = 0.0, dr;
    vec distance;
    distance.resize(_ndim);

    //#pragma omp parallel for
    for (int j = 0; j < _npart; j++) {
        if (i != j) {
            distance(0) = this->pbc(_particle(i).getposition(0, true) - _particle(j).getposition(0, true), 0);
            distance(1) = this->pbc(_particle(i).getposition(1, true) - _particle(j).getposition(1, true), 1);
            distance(2) = this->pbc(_particle(i).getposition(2, true) - _particle(j).getposition(2, true), 2);
            dr = sqrt(dot(distance, distance));
            if (dr < _r_cut) {
                f += distance(dim) * (48.0 / pow(dr, 14) - 24.0 / pow(dr, 8));
            }
        }
    }
    return f;
}

void System ::move(int i) { // Propose a MC move for particle i

    if (_sim_type == 3) {   // Gibbs sampler for Ising
        double spin_sum = _particle(pbc(i+1)).getspin() + _particle(pbc(i-1)).getspin();       
        double probability = 1 / (1 + exp(-2 * _beta * (_J * spin_sum + _H)));
        if (_rnd.Rannyu() < probability){
            _particle(i).setspin(1);
        }
        else{
            _particle(i).setspin(-1);
        }
    }

    else {                  // M(RT)^2
        if (_sim_type == 1) { // LJ system
            vec shift(_ndim); // Will store the proposed translation
            for (int j = 0; j < _ndim; j++) {
                shift(j) = _rnd.Rannyu(-1.0, 1.0) * _delta; // uniform distribution in [-_delta;_delta)
            }
            _particle(i).translate(shift, _side); // Call the function Particle::translate
            if (this->metro(i)) {                 // Metropolis acceptance evaluation
                _particle(i).acceptmove();
                _naccepted++;
            } else
                _particle(i).moveback(); // If translation is rejected, restore the old configuration
        } else {                         // Ising 1D
            if (this->metro(i)) {        // Metropolis acceptance evaluation for a spin flip involving spin i
                _particle(i).flip();     // If accepted, the spin i is flipped
                _naccepted++;
            }
        }
    }
    return;
}

bool System ::metro(int i) { // Metropolis algorithm
    bool decision = false;
    double delta_E, acceptance;
    if (_sim_type == 1)
        delta_E = this->Boltzmann(i, true) - this->Boltzmann(i, false);
    else
        delta_E = 2.0 * _particle(i).getspin() *
                  (_J * (_particle(this->pbc(i - 1)).getspin() + _particle(this->pbc(i + 1)).getspin()) + _H);
    acceptance = exp(-_beta * delta_E);
    if (_rnd.Rannyu() < acceptance)
        decision = true; // Metropolis acceptance step
    return decision;
}

double System ::Boltzmann(int i, bool xnew) {
    double energy_i = 0.0;
    double dx, dy, dz, dr;
    for (int j = 0; j < _npart; j++) {
        if (j != i) {
            dx = this->pbc(_particle(i).getposition(0, xnew) - _particle(j).getposition(0, 1), 0);
            dy = this->pbc(_particle(i).getposition(1, xnew) - _particle(j).getposition(1, 1), 1);
            dz = this->pbc(_particle(i).getposition(2, xnew) - _particle(j).getposition(2, 1), 2);
            dr = dx * dx + dy * dy + dz * dz;
            dr = sqrt(dr);
            if (dr < _r_cut) {
                energy_i += 1.0 / pow(dr, 12) - 1.0 / pow(dr, 6);
            }
        }
    }
    return 4.0 * energy_i;
}

double System ::pbc(double position, int i) { // Enforce periodic boundary conditions
    return position - _side(i) * rint(position / _side(i));
}

int System ::pbc(int i) { // Enforce periodic boundary conditions for spins
    if (i >= _npart)
        i = i - _npart;
    else if (i < 0)
        i = i + _npart;
    return i;
}

void System ::initialize(int index, string name) { // Initialize the System object according to the content of the input files in the ../INPUT/ directory

    _filename = name;

    int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG
    ifstream Primes(_filename + "Primes");
    for (int i=0; i<index; ++i){
        int trash;
        Primes >> trash >> trash;
    }
    Primes >> p1 >> p2;
    Primes.close();
    int seed[4]; // Read the seed of the RNG
    ifstream Seed(_filename + "seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    seed[3] += name.size();
    _rnd.SetRandom(seed, p1, p2);

    // ofstream couta("../OUTPUT/acceptance.dat"); // Set the heading line in file ../OUTPUT/acceptance.dat
    // stringstream couta;
    _a_buffer << "#   N_BLOCK:  ACCEPTANCE:"
             << "\n";
    // _a_buffer += couta.str();
    // couta.close();
 
    ifstream input(_filename + "input.dat"); // Start reading ../INPUT/input.dat
    // cerr<<"file open with success\n";
    // stringstream coutf;
    // coutf.open("../OUTPUT/output.dat");
    string property;
    double delta;
    while (!input.eof()) {
        input >> property;
        if (property == "SIMULATION_TYPE") {
            input >> _sim_type;
            if (_sim_type > 1) {
                input >> _J;
                input >> _H;
            }
            if (_sim_type > 3) {
                cerr << "PROBLEM: unknown simulation type"
                     << "\n";
                exit(EXIT_FAILURE);
            }
            if (_sim_type == 0)
                _out_buffer << "LJ MOLECULAR DYNAMICS (NVE) SIMULATION"
                           << "\n";
            else if (_sim_type == 1)
                _out_buffer << "LJ MONTE CARLO (NVT) SIMULATION"
                           << "\n";
            else if (_sim_type == 2)
                _out_buffer << "ISING 1D MONTE CARLO (MRT^2) SIMULATION"
                           << "\n";
            else if (_sim_type == 3)
                _out_buffer << "ISING 1D MONTE CARLO (GIBBS) SIMULATION"
                           << "\n";
        } else if (property == "RESTART") {
            input >> _restart;
        } else if (property == "TEMP") {
            input >> _temp;
            _temp += 0.1 * index;
            _beta = 1.0 / _temp;
            _out_buffer << "TEMPERATURE= " << _temp << "\n";
        } else if (property == "NPART") {
            input >> _npart;
            _fx.resize(_npart);
            _fy.resize(_npart);
            _fz.resize(_npart);
            _particle.set_size(_npart);
            for (int i = 0; i < _npart; i++) {
                _particle(i).initialize();
                if (_rnd.Rannyu() > 0.5)
                    _particle(i).flip(); // to randomize the spin configuration
            }
            _out_buffer << "NPART= " << _npart << "\n";
        } else if (property == "RHO") {
            input >> _rho;
            _volume = _npart / _rho;
            _side.resize(_ndim);
            _halfside.resize(_ndim);
            double side = pow(_volume, 1.0 / 3.0);
            for (int i = 0; i < _ndim; i++)
                _side(i) = side;
            _halfside = 0.5 * _side;
            _out_buffer << "SIDE5 ";
            for (int i = 0; i < _ndim; i++) {
                _out_buffer << setw(12) << _side[i];
            }
            _out_buffer << "\n";
        } else if (property == "R_CUT") {
            input >> _r_cut;
            _out_buffer << "R_CUT= " << _r_cut << "\n";
        } else if (property == "DELTA") {
            input >> delta;
            _out_buffer << "DELTA= " << delta << "\n";
            _delta = delta;
        } else if (property == "NBLOCKS") {
            input >> _nblocks;
            _out_buffer << "NBLOCKS= " << _nblocks << "\n";
        } else if (property == "NSTEPS") {
            input >> _nsteps;
            _out_buffer << "NSTEPS= " << _nsteps << "\n";
        } else if (property == "ENDINPUT") {
            _out_buffer << "Reading input completed!"
                       << "\n";
            break;
        } else
            cerr << "PROBLEM: unknown input"
                 << "\n";
    }
    input.close();
    // cerr << "file close with success\n";
    this->read_configuration();
    this->initialize_velocities();
    _out_buffer << "System initialized!"
               << "\n";
    // coutf.close();
    return;
}

void System ::initialize_velocities() {
    if (_restart and _sim_type == 0) {
        ifstream cinf;
        cinf.open(_filename + "velocities.in");
        if (cinf.is_open()) {
            double vx, vy, vz;
            for (int i = 0; i < _npart; i++) {
                cinf >> vx >> vy >> vz;
                _particle(i).setvelocity(0, vx);
                _particle(i).setvelocity(1, vy);
                _particle(i).setvelocity(2, vz);
            }
        } else
            cerr << "PROBLEM: Unable to open INPUT file velocities.in"
                 << "\n";
        cinf.close();
    } else {
        vec vx(_npart), vy(_npart), vz(_npart);
        vec sumv(_ndim);
        sumv.zeros();
        for (int i = 0; i < _npart; i++) {
            vx(i) = _rnd.Gauss(0., sqrt(_temp));
            vy(i) = _rnd.Gauss(0., sqrt(_temp));
            vz(i) = _rnd.Gauss(0., sqrt(_temp));
            sumv(0) += vx(i);
            sumv(1) += vy(i);
            sumv(2) += vz(i);
        }
        for (int idim = 0; idim < _ndim; idim++)
            sumv(idim) = sumv(idim) / double(_npart);
        double sumv2 = 0.0, scalef;
        for (int i = 0; i < _npart; i++) {
            vx(i) = vx(i) - sumv(0);
            vy(i) = vy(i) - sumv(1);
            vz(i) = vz(i) - sumv(2);
            sumv2 += vx(i) * vx(i) + vy(i) * vy(i) + vz(i) * vz(i);
        }
        sumv2 /= double(_npart);
        scalef = sqrt(3.0 * _temp / sumv2); // velocity scale factor
        for (int i = 0; i < _npart; i++) {
            _particle(i).setvelocity(0, vx(i) * scalef);
            _particle(i).setvelocity(1, vy(i) * scalef);
            _particle(i).setvelocity(2, vz(i) * scalef);
        }
    }
    if (_sim_type == 0) {
        double xold, yold, zold;
        for (int i = 0; i < _npart; i++) {
            xold = this->pbc(_particle(i).getposition(0, true) - _particle(i).getvelocity(0) * _delta, 0);
            yold = this->pbc(_particle(i).getposition(1, true) - _particle(i).getvelocity(1) * _delta, 1);
            zold = this->pbc(_particle(i).getposition(2, true) - _particle(i).getvelocity(2) * _delta, 2);
            _particle(i).setpositold(0, xold);
            _particle(i).setpositold(1, yold);
            _particle(i).setpositold(2, zold);
        }
    }
    if (_sim_type > 1) {
        for (int i = 0; i < _npart; i++) {
            _particle(i).setvelocity(0, 0.0);
            _particle(i).setvelocity(1, 0.0);
            _particle(i).setvelocity(2, 0.0);
        }
    }
    return;
}

void System ::initialize_buffer(string buf, int which) {

    switch(which){

        case 0:
            _out_buffer << buf;
            break;
        case 1:
            _kinetic_buffer << buf;
            break;
        case 2:
            _potential_buffer << buf;
            break;
        case 3:
            _total_buffer << buf;
            break;
        case 4:
            _temperature_buffer << buf;
            break;
        case 5:
            _pressure_buffer << buf;
            break;
        case 6:
            _magnet_buffer << buf;
            break;
        case 7:
            _cv_buffer << buf;
            break;
        case 8:
            _chi_buffer << buf;
            break;
        case 9:
            _a_buffer << buf;
            break;

        default:
            break;
    }

    return;
}

string System ::print_buffer(int which) const {

    string buf;
    switch(which){

        case 0:
            buf = _out_buffer.str();
            break;
        case 1:
            buf = _kinetic_buffer.str();
            break;
        case 2:
            buf = _potential_buffer.str();
            break;
        case 3:
            buf = _total_buffer.str();
            break;
        case 4:
            buf = _temperature_buffer.str();
            break;
        case 5:
            buf = _pressure_buffer.str();
            break;
        case 6:
            buf = _magnet_buffer.str();
            break;
        case 7:
            buf = _cv_buffer.str();
            break;
        case 8:
            buf = _chi_buffer.str();
            break;
        case 9:
            buf = _a_buffer.str();
            break;

        default:
            break;
    }

    return buf;
}


void System ::initialize_properties(bool verbose) { // Initialize data members used for measurement of properties

    string property;
    int index_property = 0;
    _nprop = 0;

    _measure_penergy = false; // Defining which properties will be measured
    _measure_kenergy = false;
    _measure_tenergy = false;
    _measure_pressure = false;
    _measure_gofr = false;
    _measure_magnet = false;
    _measure_cv = false;
    _measure_chi = false;

    ifstream input(_filename + "properties.dat");
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "POTENTIAL_ENERGY") {
                // ofstream coutp("../OUTPUT/potential_energy.dat");
                // stringstream coutp;
                if (verbose) _potential_buffer << "BLOCK:\tACTUAL_PE:\tPE_AVE:\tERROR:\n";
                // _potential_buffer += coutp.str();
                // coutp.close();
                _nprop++;
                _index_penergy = index_property;
                _measure_penergy = true;
                index_property++;
                _vtail = 8 * M_PI * _rho * ( 1 / (3 * pow(_r_cut, 6)) - 1) / (3 * pow(_r_cut, 3));
            } else if (property == "KINETIC_ENERGY") {
                // ofstream coutk("../OUTPUT/kinetic_energy.dat");
                // stringstream coutk;
                if (verbose) _kinetic_buffer << "BLOCK:\tACTUAL_KE:\tKE_AVE:\tERROR:\n";
                // _kinetic_buffer += coutk.str();
                // coutk.close();
                _nprop++;
                _measure_kenergy = true;
                _index_kenergy = index_property;
                index_property++;
            } else if (property == "TOTAL_ENERGY") {
                // ofstream coutt("../OUTPUT/total_energy.dat");
                // stringstream coutt;
                if (verbose) _total_buffer << "BLOCK:\tACTUAL_TE:\tTE_AVE:\tERROR:\n";
                // _total_buffer += coutt.str();
                // coutt.close();
                _nprop++;
                _measure_tenergy = true;
                _index_tenergy = index_property;
                index_property++;
            } else if (property == "TEMPERATURE") {
                // ofstream coutte("../OUTPUT/temperature.dat");
                // stringstream coutte;
                if (verbose) _temperature_buffer << "BLOCK:\tACTUAL_T:\tT_AVE:\tERROR:\n";
                // _temperature_buffer += coutte.str();
                // coutte.close();
                _nprop++;
                _measure_temp = true;
                _index_temp = index_property;
                index_property++;
            } else if (property == "PRESSURE") {
                // ofstream coutpr("../OUTPUT/pressure.dat");
                // stringstream coutpr;
                if (verbose) _pressure_buffer << "BLOCK:\tACTUAL_P:\tP_AVE:\tERROR:\n";
                // _pressure_buffer += coutpr.str();
                // coutpr.close();
                _nprop++;
                _measure_pressure = true;
                _index_pressure = index_property;
                index_property++;
                _ptail = 32.0 * M_PI * _rho * ( 1/(9 * pow(_r_cut, 9)) - 1/(6 * pow(_r_cut, 3)) );
            } else if (property == "GOFR") {
                _gofr_buffer << "BLOCK:\tDISTANCE:\tACTUAL_G:\tAVE_GOFR:\tERROR:\n";
                //ofstream coutgr("../OUTPUT/gofr.dat");
                //coutgr << "# DISTANCE:     AVE_GOFR:        ERROR:" << "\n";
                //coutgr.close();
                input >> _n_bins;
                _nprop += _n_bins;
                _bin_size = (_halfside.min()) / (double)_n_bins;
                _measure_gofr = true;
                _index_gofr = index_property;
                index_property += _n_bins;
            } else if (property == "MAGNETIZATION") {
                //ofstream coutpr("../OUTPUT/magnetization.dat");
                if (verbose) _magnet_buffer << "BLOCK:\tACTUAL_M:\tM_AVE:\tERROR:\n";
                //coutpr.close();
                _nprop++;
                _measure_magnet = true;
                _index_magnet = index_property;
                index_property++;
            } else if (property == "SPECIFIC_HEAT") {
                //ofstream coutpr("../OUTPUT/specific_heat.dat");
                if (verbose) _cv_buffer << "BLOCK:\tACTUAL_CV:\tCV_AVE:\tERROR:\n";
                //coutpr.close();
                _nprop++;
                _measure_cv = true;
                _index_cv = index_property;
                index_property++;
            } else if (property == "SUSCEPTIBILITY") {
                //ofstream coutpr("../OUTPUT/susceptibility.dat");
                if (verbose) _chi_buffer << "BLOCK:\tACTUAL_X:\tX_AVE:\tERROR:\n";
                //coutpr.close();
                _nprop++;
                _measure_chi = true;
                _index_chi = index_property;
                index_property++;
            } else if (property == "ENDPROPERTIES") {
                stringstream coutf;
                // coutf.open("../OUTPUT/output.dat",ios::app);
                _out_buffer << "Reading properties completed!"
                           << "\n";
                // _out_buffer += coutf.str();
                // coutf.close();
                break;
            } else
                cerr << "PROBLEM: unknown property"
                     << "\n";
        }
        input.close();
    } else
        cerr << "PROBLEM: Unable to open properties.dat"
             << "\n";

    // according to the number of properties, resize the vectors _measurement,_average,_block_av,_global_av,_global_av2
    _measurement.resize(_nprop);
    _average.resize(_nprop);
    _block_av.resize(_nprop);
    _global_av.resize(_nprop);
    _global_av2.resize(_nprop);
    _average.zeros();
    _global_av.zeros();
    _global_av2.zeros();
    _nattempts = 0;
    _naccepted = 0;
    return;
}

void System ::finalize() {
    this->write_configuration();
    _rnd.SaveSeed();
    string outname = _filename + "out/";
    ofstream coutf;

    coutf.open(outname + "output.dat");
    coutf << _out_buffer.str();
    coutf << "Simulation completed!"
          << "\n";
    coutf.close();

    if (_sim_type != 0){  //only to save a few MB of space in the filesystem
        coutf.open(outname + "acceptance.dat");
        coutf << _a_buffer.str();
        coutf.close();
    }

    if (_measure_kenergy){
        coutf.open(outname + "kinetic_energy.dat");
        coutf << _kinetic_buffer.str();
        coutf.close();
    }

    if (_measure_penergy){
        coutf.open(outname + "potential_energy.dat");
        coutf << _potential_buffer.str();
        coutf.close();
    }

    if (_measure_tenergy){
        coutf.open(outname + "total_energy.dat");
        coutf << _total_buffer.str();
        coutf.close();
    }

    if (_measure_temp){
        coutf.open(outname + "temperature.dat");
        coutf << _temperature_buffer.str();
        coutf.close();
    }

    if (_measure_pressure){
        coutf.open(outname + "pressure.dat");
        coutf << _pressure_buffer.str();
        coutf.close();
    }

    if (_measure_magnet){
        coutf.open(outname + "magnetization.dat");
        coutf << _magnet_buffer.str();
        coutf.close();
    }

    if (_measure_cv){
        coutf.open(outname + "specific_heat.dat");
        coutf << _cv_buffer.str();
        coutf.close();
    }

    if (_measure_chi){
        coutf.open(outname + "susceptibility.dat");
        coutf << _chi_buffer.str();
        coutf.close();
    }

    if (_measure_gofr){
        coutf.open(outname + "gofr.dat");
        coutf << _gofr_buffer.str();
        coutf.close();
    }

    return;
}

// Write current configuration as a .xyz file
void System ::write_configuration() {
    ofstream coutf;
    if (_sim_type < 2) {
        coutf.open(_filename + "config.out.xyz");
        if (coutf.is_open()) {
            coutf << _npart << "\n";
            coutf << "#Comment!"
                  << "\n";
            for (int i = 0; i < _npart; i++) {
                coutf << "LJ"
                      << "  "
                      << setw(16) << _particle(i).getposition(0, true) / _side(0)          // x
                      << setw(16) << _particle(i).getposition(1, true) / _side(1)          // y
                      << setw(16) << _particle(i).getposition(2, true) / _side(2) << "\n"; // z
            }
        } else
            cerr << "PROBLEM: Unable to open config.out.xyz"
                 << "\n";
        coutf.close();
        this->write_velocities();
    } else {
        coutf.open(_filename + "config.out.spin");
        for (int i = 0; i < _npart; i++)
            coutf << _particle(i).getspin() << " ";
        coutf.close();
    }
    return;
}

// Write configuration nconf as a .xyz file in directory ../OUTPUT/CONFIG/
void System ::write_XYZ(int nconf) {
    ofstream coutf;
    coutf.open(_filename + "config_" + to_string(nconf) + ".out.xyz");
    if (coutf.is_open()) {
        coutf << _npart << "\n";
        coutf << "#Comment!"
              << "\n";
        for (int i = 0; i < _npart; i++) {
            coutf << "LJ"
                  << "  "
                  << setw(16) << _particle(i).getposition(0, true)          // x
                  << setw(16) << _particle(i).getposition(1, true)          // y
                  << setw(16) << _particle(i).getposition(2, true) << "\n"; // z
        }
    } else
        cerr << "PROBLEM: Unable to open config.out.xyz"
             << "\n";
    coutf.close();
    return;
}

void System ::write_velocities() {
    ofstream coutf;
    coutf.open(_filename + "velocities.out");
    if (coutf.is_open()) {
        for (int i = 0; i < _npart; i++) {
            coutf << setw(16) << _particle(i).getvelocity(0)          // vx
                  << setw(16) << _particle(i).getvelocity(1)          // vy
                  << setw(16) << _particle(i).getvelocity(2) << "\n"; // vz
        }
    } else
        cerr << "PROBLEM: Unable to open velocities.out"
             << "\n";
    coutf.close();
    return;
}

// Read configuration from a .xyz file
void System ::read_configuration() {
    ifstream cinf;
    cinf.open(_filename + "config.in.xyz");
    if (cinf.is_open()) {
        string comment;
        string particle;
        double x, y, z;
        int ncoord;
        cinf >> ncoord;
        if (ncoord != _npart) {
            cerr << "PROBLEM: conflicting number of coordinates in input file & config file not match!"
                 << "\n";
            exit(EXIT_FAILURE);
        }
        cinf >> comment;
        for (int i = 0; i < _npart; i++) {
            cinf >> particle >> x >> y >> z; // units of coordinates in conf.xyz is _side
            _particle(i).setposition(0, this->pbc(_side(0) * x, 0));
            _particle(i).setposition(1, this->pbc(_side(1) * y, 1));
            _particle(i).setposition(2, this->pbc(_side(2) * z, 2));
            _particle(i).acceptmove(); // _x_old = _x_new
        }
    } else
        cerr << "PROBLEM: Unable to open INPUT config file"
             << "\n";
    cinf.close();

    if (_restart and _sim_type > 1) {
        int spin;
        cinf.open(_filename + "config.spin");
        for (int i = 0; i < _npart; i++) {
            cinf >> spin;
            _particle(i).setspin(spin);
        }
        cinf.close();
    }
    return;
}

void System ::block_reset(int blk) { // Reset block accumulators to zero
    // stringstream coutf;
    if (blk > 0) {
        //_out_buffer << "Block completed: " << blk << "\n"; //useless ...
    }
    _block_av.zeros();
    return;
}

void System ::measure() { // Measure properties
    _measurement.zeros();
    // POTENTIAL ENERGY, VIRIAL, GOFR ///////////////////////////////////////////
    int bin;
    vec distance;
    distance.resize(_ndim);
    double penergy_temp = 0.0, dr; // temporary accumulator for potential energy
    double kenergy_temp = 0.0;     // temporary accumulator for kinetic energy
    double tenergy_temp = 0.0, tenergy_temp_2 = 0.0;
    double magnetization = 0.0;
    double pressure_temp{};
    double virial = 0.0;
    if (_measure_penergy or _measure_pressure or _measure_gofr) {
        for (int i = 0; i < _npart - 1; i++) {
            for (int j = i + 1; j < _npart; j++) {
                distance(0) = this->pbc(_particle(i).getposition(0, true) - _particle(j).getposition(0, true), 0);
                distance(1) = this->pbc(_particle(i).getposition(1, true) - _particle(j).getposition(1, true), 1);
                distance(2) = this->pbc(_particle(i).getposition(2, true) - _particle(j).getposition(2, true), 2);
                dr = sqrt(dot(distance, distance));
                bin = static_cast<int>(floor(dr/_bin_size)); 
                if (_measure_gofr && (bin < _n_bins))
                    _measurement(_index_gofr + bin) += 2;
                if (dr < _r_cut) {
                    if (_measure_penergy)
                        penergy_temp += 1.0 / pow(dr, 12) - 1.0 / pow(dr, 6); // POTENTIAL ENERGY
                    if (_measure_pressure)
                        pressure_temp += 1.0 / pow(dr, 12) - 0.5 / pow(dr, 6); // PRESSURE
                }
            }
        }
    }
    // GOFR //////////////////////////////////////////////////////////////////////
    if (_measure_gofr) {
        for (int bin = 0; bin < _n_bins; ++bin){
            double shell = 4 * M_PI * (pow(bin * _bin_size + 1, 3) - pow(bin * _bin_size, 3)) / 3;
            _measurement(_index_gofr + bin) /= _rho * _npart * shell;
        }
    }
    // POTENTIAL ENERGY //////////////////////////////////////////////////////////
    if (_measure_penergy) {
        penergy_temp = _vtail + 4.0 * penergy_temp / double(_npart);
        _measurement(_index_penergy) = penergy_temp;
    }
    // KINETIC ENERGY ////////////////////////////////////////////////////////////
    if (_measure_kenergy) {
        for (int i = 0; i < _npart; i++)
            kenergy_temp += 0.5 * dot(_particle(i).getvelocity(), _particle(i).getvelocity());
        kenergy_temp /= double(_npart);
        _measurement(_index_kenergy) = kenergy_temp;
    }
    // TOTAL ENERGY (kinetic+potential) //////////////////////////////////////////
    if (_measure_tenergy) {
        if (_sim_type < 2)
            _measurement(_index_tenergy) = kenergy_temp + penergy_temp;
        else {
            double s_i, s_j;
            for (int i = 0; i < _npart; i++) {
                s_i = double(_particle(i).getspin());
                s_j = double(_particle(this->pbc(i + 1)).getspin());
                double temp = -_J * s_i * s_j - 0.5 * _H * (s_i + s_j);
                tenergy_temp += temp;
                tenergy_temp_2 += pow(temp, 2);
            }
            tenergy_temp /= double(_npart);
            tenergy_temp_2 /= double(_npart);
            _measurement(_index_tenergy) = tenergy_temp;
        }
    }
    // TEMPERATURE ///////////////////////////////////////////////////////////////
    if (_measure_temp) {
        if (_measure_kenergy)
            _measurement(_index_temp) = (2.0 / 3.0) * kenergy_temp;
        else {
            for (int i = 0; i < _npart; i++)
                kenergy_temp += 0.5 * dot(_particle(i).getvelocity(), _particle(i).getvelocity());
            kenergy_temp /= double(_npart);
            _measurement(_index_temp) = (2.0 / 3.0) * kenergy_temp;
        }
    }
    // PRESSURE //////////////////////////////////////////////////////////////////
    if (_measure_temp and _measure_kenergy and _measure_pressure) {
        pressure_temp /= double(_npart);
        _measurement(_index_pressure) = _rho * (2.0 / 3.0) * kenergy_temp + (48.0 / 3.0) * pressure_temp / _volume;
    }
    // MAGNETIZATION /////////////////////////////////////////////////////////////
    if (_measure_magnet) {
        double s{};
        for (int i = 0; i < _npart; ++i) {
            s += double(_particle(i).getspin());
        }
        _measurement(_index_magnet) = s / double(_npart);
    }
    // SPECIFIC HEAT /////////////////////////////////////////////////////////////
    if (_measure_tenergy and _measure_cv) {
        _measurement(_index_cv) = pow(_beta, 2) * (tenergy_temp_2 - pow(tenergy_temp, 2));
    }
    // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
    if (_measure_chi) {
        double s{};
        for (int i = 0; i < _npart; ++i) {
            s += double(_particle(i).getspin());
        }
        _measurement(_index_chi) = _beta * pow(s, 2) / _npart;
    }

    _block_av += _measurement; // Update block accumulators

    return;
}

void System ::averages(int blk) {

    // stringstream coutf;
    double average, sum_average, sum_ave2;

    _average = _block_av / double(_nsteps);
    _global_av += _average;
    _global_av2 += _average % _average; // % -> element-wise multiplication

    // POTENTIAL ENERGY //////////////////////////////////////////////////////////
    if (_measure_penergy) {
        // coutf.open("../OUTPUT/potential_energy.dat",ios::app);
        average = _average(_index_penergy);
        sum_average = _global_av(_index_penergy);
        sum_ave2 = _global_av2(_index_penergy);
        _potential_buffer << blk
                         << "\t" << average
                         << "\t" << sum_average / double(blk)
                         << "\t" << this->error(sum_average, sum_ave2, blk) << "\n";
        // _potential_buffer += coutf.str();
        // coutf.clear();
        // coutf.close();
    }
    // KINETIC ENERGY ////////////////////////////////////////////////////////////
    if (_measure_kenergy) {
        // coutf.open("../OUTPUT/kinetic_energy.dat",ios::app);
        average = _average(_index_kenergy);
        sum_average = _global_av(_index_kenergy);
        sum_ave2 = _global_av2(_index_kenergy);
        _kinetic_buffer << blk
                       << "\t" << average
                       << "\t" << sum_average / double(blk)
                       << "\t" << this->error(sum_average, sum_ave2, blk) << "\n";
        // _kinetic_buffer += coutf.str();
        // coutf.clear();
        // coutf.close();
    }
    // TOTAL ENERGY //////////////////////////////////////////////////////////////
    if (_measure_tenergy) {
        // coutf.open("../OUTPUT/total_energy.dat",ios::app);
        average = _average(_index_tenergy);
        sum_average = _global_av(_index_tenergy);
        sum_ave2 = _global_av2(_index_tenergy);
        _total_buffer << blk
                     << "\t" << average
                     << "\t" << sum_average / double(blk)
                     << "\t" << this->error(sum_average, sum_ave2, blk) << "\n";
        // _total_buffer += coutf.str();
        // coutf.clear();
        // coutf.close();
    }
    // TEMPERATURE ///////////////////////////////////////////////////////////////
    if (_measure_temp) {
        average = _average(_index_temp);
        sum_average = _global_av(_index_temp);
        sum_ave2 = _global_av2(_index_temp);
        _temperature_buffer << blk
                           << "\t" << average
                           << "\t" << sum_average / double(blk)
                           << "\t" << this->error(sum_average, sum_ave2, blk) << "\n";
    }
    // PRESSURE //////////////////////////////////////////////////////////////////
    if (_measure_pressure) {
        average = _average(_index_pressure);
        sum_average = _global_av(_index_pressure);
        sum_ave2 = _global_av2(_index_pressure);
        _pressure_buffer << blk
                        << "\t" << average
                        << "\t" << sum_average / double(blk)
                        << "\t" << this->error(sum_average, sum_ave2, blk) << "\n";
    }
    // GOFR //////////////////////////////////////////////////////////////////////
    if (_measure_gofr) {
        for (int i = 0; i < _n_bins; ++i){
            average = _average(_index_gofr + i);
            sum_average = _global_av(_index_gofr + i);
            sum_ave2 = _global_av2(_index_gofr + i);
            _gofr_buffer << blk << "\t" << (i+0.5)*_bin_size
                        << "\t" << average
                        << "\t" << sum_average / double(blk)
                        << "\t" << this->error(sum_average, sum_ave2, blk) << "\n";
        }
    }
    // MAGNETIZATION /////////////////////////////////////////////////////////////
    if (_measure_magnet) {
        average = _average(_index_magnet);
        sum_average = _global_av(_index_magnet);
        sum_ave2 = _global_av2(_index_magnet);
        _magnet_buffer << blk
                      << "\t" << average
                      << "\t" << sum_average / double(blk)
                      << "\t" << this->error(sum_average, sum_ave2, blk) << "\n";
    }
    // SPECIFIC HEAT /////////////////////////////////////////////////////////////
    if (_measure_cv) {
        average = _average(_index_cv);
        sum_average = _global_av(_index_cv);
        sum_ave2 = _global_av2(_index_cv);
        _cv_buffer << blk
                  << "\t" << average
                  << "\t" << sum_average / double(blk)
                  << "\t" << this->error(sum_average, sum_ave2, blk) << "\n";
    }
    // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
    if (_measure_chi) {
        average = _average(_index_chi);
        sum_average = _global_av(_index_chi);
        sum_ave2 = _global_av2(_index_chi);
        _chi_buffer << blk
                   << "\t" << average
                   << "\t" << sum_average / double(blk)
                   << "\t" << this->error(sum_average, sum_ave2, blk) << "\n";
    }
    // ACCEPTANCE ////////////////////////////////////////////////////////////////
    double fraction;
    // coutf.open("../OUTPUT/acceptance.dat",ios::app);
    if (_nattempts > 0)
        fraction = double(_naccepted) / double(_nattempts);
    else
        fraction = 0.0;
    _a_buffer << setw(12) << blk << setw(12) << fraction << "\n";
    // _a_buffer += coutf.str();
    // while(coutf >> _a_buffer);
    // coutf.close();

    return;
}

double System ::error(double acc, double acc2, int blk) {
    if (blk <= 1)
        return 0.0;
    else
        return sqrt(fabs(acc2 / double(blk) - pow(acc / double(blk), 2)) / double(blk));
}

int System ::get_nbl() {
    return _nblocks;
}

int System ::get_nsteps() {
    return _nsteps;
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
