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

#ifndef __BARRIER_H
#define __BARRIER_H

enum EvolutionType {reflect, absorb};
enum InitialType   {zero, stable, spike};

const int maxsteps = 2200;
const int footprint = maxsteps + 1;

class BarrierDistribution {
  
public:
  BarrierDistribution(double,InitialType,int shift=0);
  BarrierDistribution(BarrierDistribution const &prev,EvolutionType);
  void show() const;  
  double pdensity() const;
  double tdensity() const;
  const double p;
  friend BarrierDistribution* convolve(const BarrierDistribution*,
				       const BarrierDistribution*);
  friend BarrierDistribution* translate(const BarrierDistribution*,
					int);
  
private:
  double sites[footprint + 1];
  int    rcheck(int) const;
  double stationary(int) const;
  double spikedist(int,int) const;
  double evolve_absorbtion(int) const;
  double evolve_reflection(int) const;
  double get(int) const;
  void   set(int, double);
};

#endif
