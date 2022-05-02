/*
Copyright [2019] XXX

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
#include "barriertools.h"

using namespace std;

int main()
{
  double p;
  int k, w, step;
  BarrierDistribution* stationary;
  BarrierDistribution* spikeshift;
  BarrierDistribution* distributions[maxsteps+1];
  
  cout << "Enter binomial distribution parameter: ";
  cin >> p;
  
  cout << "Enter spike power: ";
  cin >> k;
  
  cout << "Enter walk length (no more than " << maxsteps << "): ";
  cin >> w;
  
  if ((k >= 0) && (w >= 0) && (w <= maxsteps)) {
    stationary = new BarrierDistribution(p,stable);
    spikeshift = new BarrierDistribution(p,spike,k);
    distributions[0] = convolve(stationary,spikeshift);
    delete(stationary);
    delete(spikeshift);
    cout << "Results, of form (length, uncaptured probability)." << "\n";
    for (step = 1; step <= w; step++) {
      distributions[step] = new BarrierDistribution(*distributions[step-1],absorb);
      delete(distributions[step-1]);
      cout << "(" << step << "," << distributions[step]->pdensity() << ")";
      cout << "\n";
    }
    if (w>0) delete(distributions[w]);
  }
  return 0;
}


