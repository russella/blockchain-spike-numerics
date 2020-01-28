/*
Copyright [2019] Alexander Russell

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
  double error_threshold, current_error;
  int k_lower, k_upper, w, step;
  BarrierDistribution* stationary;
  BarrierDistribution* spikeshift;
  BarrierDistribution* distributions[maxsteps+1];
  
  cout << "Enter binomial distribution parameter: ";
  cin >> p;
  cout << "Enter error threshold: ";
  cin >> error_threshold;
  
  cout << "Enter lower spike quota: ";
  cin >> k_lower;
  cout << "Enter upper spike quota: ";
  cin >> k_upper;
  
  cout << "Enter walk length, a conjectured upper bound on how many steps will be necessary to achieve this error (no more than " << maxsteps << "): ";
  cin >> w;
  
  if ((k_lower >= 0) && (w >= 0) && (w <= maxsteps)) {
    stationary = new BarrierDistribution(p,stable);
    for (int k=k_lower; k <= k_upper; k++) { 
      spikeshift = new BarrierDistribution(p,spike,k);
      distributions[0] = convolve(stationary,spikeshift);
      delete(spikeshift);
      current_error = 1;
      step = 1;
      while ((current_error > error_threshold) && (step <= w)) {
	distributions[step] = new BarrierDistribution(*distributions[step-1],absorb);
	delete(distributions[step-1]);
	current_error = distributions[step]->pdensity();
	step++;
      }
      if (current_error <= error_threshold)
	cout << "(" << k << "," << step << ")\n";
      else
	cout << "Underflow: try larger walk.\n";
    }
    delete(stationary); }
  else cout << "Bad parameters./n";
  return 0;
}
  
