/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "stats.h"

double covariance (std::vector<double> x, std::vector<double> y){
	assert (x.size() == y.size() && "to compute covariance the two input vectors need to be the same size\n");
	double avg_xy{}, avg_x{}, avg_y{};
	for (long unsigned int i=0; i<x.size(); ++i){
		avg_x += x.at(i);
		avg_y += y.at(i);
		avg_xy += x.at(i)*y.at(i);
	}
	avg_x /= (double)x.size();
	avg_y /= (double)x.size();
	avg_xy /= (double)x.size();

	return avg_xy - (avg_x * avg_y);
}

double autocorrelation (std::vector<double> x, long unsigned int lag){
	assert (lag < x.size() && "lag needs to be smaller than vector size\n");
	std::vector<double> front = x;
	std::vector<double> back = x;
	for(long unsigned int i=0; i<lag; ++i){
		front.pop_back();
		//std::reverse(back.begin(), back.end()); // first becomes last, reverses the vector
		//back.pop_back(); // pop last
		//std::reverse(back.begin(), back.end());
		back.erase(back.begin());
	}
	return covariance(front, back);
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
