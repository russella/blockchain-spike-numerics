/*
Copyright [2020] Alexander Russell

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/


#include <iostream>
#include <stdexcept>
#include <cmath>
#include "disttools.h"

using namespace std;

int main(int argc, char **argv)
{
  double hon_stake;
  double adv_stake;
  double hon_prob;
  double adv_prob;
  int delta;
  double f;
  int w, step;
  Distribution* distributions[2];
  
  if (argc < 5) {
    cout << "Usage: " << argv[0] << " <hon_stake * 100> <f * 100> <delta> <walk_length>" << endl;
    return 0;
  }

  hon_stake = atoi(argv[1])/100.0;
  f = atoi(argv[2])/100.0;
  delta = atoi(argv[3]);
  w = atoi(argv[4]);
  
  cout << "hon_stake = " << hon_stake << endl;
  cout << "f         = " << f << endl;
  cout << "delta     = " << delta << endl;
  cout << "w         = " << w << endl;
  
  //cout << "Enter honest stake ratio: (between 0 and 1): ";
  //cin  >> hon_stake;
  //cout << "Enter active slot coefficient (between 0 and 1): ";
  //cin  >> f;
  //cout << "Enter networking delay (Delta, no more than " << maxdelta <<  "): ";
  //cin  >> delta;
  
  hon_prob = 1-pow((1-f),hon_stake);
  cout << "Probability of honest success: " <<  hon_prob << "\n"; 
  cout << "Effective rate of honest advancement: " << 1/(delta-1+1/hon_prob) << "\n";
  
  adv_stake = 1 - hon_stake;
  adv_prob = 1 - pow((1-f),adv_stake);
  cout << "Adversarial stake ratio: " << adv_stake << "\n";
  cout << "Adversarial success probability: " <<  adv_prob << "\n";
  cout << "...equal to expected rate of adversarial advancement: " << adv_prob << "\n";
  

  //cout << "Enter number of steps of evolution (no more than " << maxsteps << "): ";
  //cin  >> w;
  
  distributions[0] = new Distribution(delta,identity);
  cout << "Evolution beginning...\n";
  for (step = 1; step <= w; step++) {
    distributions[step % 2] = evolve(distributions[(step - 1) % 2],
				     adv_prob, hon_prob);
    delete(distributions[(step - 1) % 2]);
    if (step % 10 == 0)
      double new_density = distributions[step % 2]->pdensity();
      //    double new_density_t = distributions[step % 2]->tdensity();
      cout << "(" 
        << "adv. stake: " << adv_stake << ", " 
        << "f: " << f << ", " 
        << "delta: " << delta << ", " 
        << "step: " << step << ", " 
        << "density: " << new_density << ")\n" << std::flush;
  }
  delete(distributions[w % 2]);
  cout <<  "========================================================================================================================" << endl;
  return 0;
}


