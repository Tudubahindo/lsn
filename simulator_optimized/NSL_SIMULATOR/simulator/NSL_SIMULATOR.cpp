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
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {

    std::string read_cycles, read_name;

    if (argc < 3){
        read_cycles = "1";
        read_name = "../l4/md_solid/";
    }
    else {
        read_cycles = argv[1];
        read_name = argv[2];
    }
    int cycles = std::stoi(read_cycles);

    vector<string> bufs(10);

    for (int k = 0; k < cycles; ++k) {
        int nconf = 1;
        System SYS;
        SYS.initialize(k, read_name);
        SYS.initialize_properties( (k == 0) );
        SYS.block_reset(0);

        if (k > 0){
            for (int i = 0; i < bufs.size(); ++i){
                SYS.initialize_buffer(bufs.at(i), i); // che modo merdoso di fare le cose
            }
        }
    
        for (int i = 0; i < SYS.get_nbl(); i++) {        // loop over blocks
            for (int j = 0; j < SYS.get_nsteps(); j++) { // loop over steps in a block
                SYS.step();
                SYS.measure();
                if (j % 10 == 0) {
                    //        SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
                    nconf++;
                }
            }
            SYS.averages(i + 1);
            SYS.block_reset(i + 1);
        }

        for (int i = 0; i < bufs.size(); ++i){
            bufs.at(i) = SYS.print_buffer(i); // che modo STRAmerdoso di fare le cose
        }
        if (k == cycles-1){
            SYS.finalize();
        }
    }

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
