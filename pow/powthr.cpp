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
#include "barriertools.h"

using namespace std;

int main()
{
  double hon_param;
  double adv_param;
  int spike_begin, spike_end;
  double approx_error,error_threshold;
  int delta;
  int step;
  BarrierDistribution* distributions[2];
  BarrierDistribution* stationary;
  
  cout << "Enter Poisson parameter for honest distribution: ";
  cin  >> hon_param;
  
  cout << "Enter networking delay (Delta, no more than " << maxdelta <<  "): ";
  cin  >> delta;
  
  cout << "Effective honest unique, isolated probability: " << hon_param * exp(-hon_param * (2 * delta + 1)) << "\n";
  cout << "Enter Poisson parameter of adversarial success: ";
  cin  >> adv_param;
  
  cout << "Enter step-to-step approximation error : ";
  cin  >> approx_error;

  distributions[0] = new BarrierDistribution(delta,identity);
  cout << "Estimating stationary distribution...\n";
  double error = 1;
  for  (step = 1; error > approx_error; step++) {
    distributions[step % 2] = evolve(distributions[(step - 1) % 2],
				     reflect,
				     adv_param, hon_param);
    error = stat_distance(distributions[0],distributions[1]);
    delete(distributions[(step - 1) % 2]);
    cout << "[" << step << ":" << error << "]  \r" << std::flush; }
  cout << "\n";
  stationary = distributions[(step-1) % 2];
  cout << "Stationary approximation complete.\n";

  cout << "Enter desired stabilization error threshold: ";
  cin  >> error_threshold;

  cout << "Enter initial spike (an integer): ";
  cin  >> spike_begin;
  
  cout << "Enter final spike (an integer): ";
  cin  >> spike_end;
  
  for (int spike = spike_begin; spike <= spike_end; spike++) {
    distributions[0] = convolve_spike(stationary,double (spike));
    step = 0; error = 1.0;
    while (error > error_threshold) {
      step++;
      distributions[step % 2] = evolve(distributions[(step - 1) % 2],
				       absorb,
				       adv_param, hon_param);
      delete(distributions[(step - 1) % 2]);
      error = distributions[step % 2]->pdensity(); }
    cout << "(" << spike << ", " << step << ")\n" << std::flush;
    delete(distributions[step % 2]); }
  delete(stationary);
  return 0;
}


