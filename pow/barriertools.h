/*
Copyright [2020] XXX

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

#ifndef __BARRIER_H
#define __BARRIER_H

enum EvolutionType {reflect, absorb};
enum InitializationType {zero, identity};

const int maxsteps = 200;
const int footprint = maxsteps + 1;
const int maxdelta = 30;

class Dist_index {
public:
  const int delta;
  Dist_index(int,   //  initial delta
	     int,   //  initial margin
	     int,   //  initial r_isolation
	     bool); //  initial l_isolated_pending 
  Dist_index* evolve(int,
		     int) const;
  int get_internal() const;
  int get_beta() const;
  void set(int,   //  margin
	   int,   //  r_isolation
	   bool); //  l_isolated_pending
private:
  int  beta;
  int  r_isolation;
  bool l_isolated_pending;
};

class BarrierDistribution {
  
public:
  const int delta;
  BarrierDistribution(int,InitializationType); // delta
  BarrierDistribution(const BarrierDistribution*); // copy constructor
  //
  void show() const;  
  double pdensity() const;
  double tdensity() const;
  friend double stat_distance(const BarrierDistribution*,
			      const BarrierDistribution*);
  friend BarrierDistribution* convolve_spike(const BarrierDistribution*,
					     double);  // spike param
  friend BarrierDistribution* evolve(const BarrierDistribution*,
				     EvolutionType,
				     double,  // adversarial Poisson param
				     double); // honest prob
  
private:
  double sites[footprint][2*maxdelta + 2];
  int    rcheck_beta(int) const;
  int    rcheck_delta(int) const;
  // accessor functions
  double get(Dist_index*) const;    // index object
  void   set(Dist_index*, double);  // index object, value
};

#endif
