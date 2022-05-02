/*
COPYRIGHT 2022 XXX

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


/*
  The BarrierDistribution object can hold a single probability
  distribution supported on {0, ..., maxsteps}.
  
  Given a particular distribution and a parameter 0 < p < 1 (which the
  object holds), one can "evolve" the distribution to yield the next
  one under the random walk with barrier: this means that a particle
  with position given by the distribution moves up with probability p
  and down with probability (1-p); the exception appears at zero,
  where the particle cannot move down, and stays put if "asked to move
  down."  Indeed, the BarrierDistribution object is equipped with a
  constructor function which, given a BarrierDistribution object,
  initializes the new object to be the result of evolving the old one
  for one step.
  
  Note that there are a number of ways one could accelerate this
  code. In particular, distribution objects carry out (hopefully)
  unnecessary range checking, are conservatively initialized to zero,
  and are created anew for each step of the evolution.
*/

//class BarrierDistribution;

double Poisson(double lambda, int k) {
  if (k < 0) throw std::invalid_argument("Poisson index out of range");
  double result = 1.0;
  for (int i=1; i <= k; i++)
    result = result * (lambda / i);
  return( exp(-lambda) * result );
}

void BarrierDistribution::show() const {
  cout << "Distribution contents...\n";
  for (int beta=0; beta < 30; beta++) {
    cout << beta << " : " << sites[beta] << "\n";
  }
}

int BarrierDistribution::rcheck(int beta) const {
  if ((beta < 0) || (beta > footprint)) throw std::invalid_argument("distribution index out of range");
  return beta; }

double BarrierDistribution::get(int beta) const {
  return(sites[rcheck(beta)]);
}

void BarrierDistribution::set(int beta, double value) {
  sites[rcheck(beta)] = value;
}

double BarrierDistribution::evolve_reflection(int beta) const {
  if (beta > 0)
    return(p * get(beta-1) + (1-p) * get(beta+1));
  if (beta == 0)
    return((1-p) * (get(0) + get(1)));
  else
    throw std::invalid_argument("distribution index out of range");
}

double BarrierDistribution::evolve_absorbtion(int beta) const {
  if (beta > 1)
    return(p * get(beta-1) + (1-p) * get(beta+1));
  if (beta == 1)
    return((1-p) * get(2));
  if (beta == 0)
    return(get(0) + (1-p)*get(1));
  else
    throw std::invalid_argument("distribution index out of range");
}

double BarrierDistribution::pdensity() const {
  double result = 0.0;
  int beta;
  for(beta = 1; beta <= footprint; beta++)
    result = result + get(beta);
  if (result < 1)
    return(result);
  else
    return(1);
}

double BarrierDistribution::tdensity() const {
  double result = 0.0;
  int beta;
  for(beta = 0; beta <= footprint; beta++)
    result = result + get(beta);
  return(result);
}

double BarrierDistribution::stationary(int t) const {
  return( pow(p/(1-p),t) * (1-2*p)/ (1-p));
}

double spiketail(int shift, int beta) {
  if (beta <= shift)
    return(1.0);
  else
    return(pow(double(shift)*exp(1.0)/beta,beta) * exp(double(-shift))); }

double BarrierDistribution::spikedist(int shift,int beta) const {
  return(spiketail(shift,beta) - spiketail(shift,beta+1));
}

BarrierDistribution::BarrierDistribution(double parameter,
					 InitialType selection,
					 int shift) : p(parameter) {
  int beta;
  if (selection==stable) {
    for(beta = 0; beta <= footprint; beta++)
      set(beta,stationary(beta)); }
  else if (selection==zero) {
    for(beta = 0; beta <= footprint; beta++)
      set(beta,0.0); }
  else if (selection==spike) {
    for (beta = 0; beta<= footprint; beta++)
      set(beta,spikedist(shift,beta)); }
  else throw std::invalid_argument("initial type unknown");
}

BarrierDistribution::BarrierDistribution(BarrierDistribution const &prev, EvolutionType eselect) : p(prev.p) {
  int beta;
  set(footprint,0);
  if (eselect == reflect) {
    for(beta = 0; beta < footprint; beta++)
      set(beta,prev.evolve_reflection(beta));
    set(footprint,0.0); }
  else if (eselect == absorb) {
    for(beta = 0; beta < footprint; beta++)
      set(beta,prev.evolve_absorbtion(beta));
    set(footprint,0.0); }
  else     
    std::invalid_argument("evolution type unknown");
}

BarrierDistribution* translate(const BarrierDistribution* source,
			       int k) {
  BarrierDistribution* result;
  result = new  BarrierDistribution(source->p,zero);
  for (int beta=k; beta <= footprint; beta++) {
    result->set(beta,source->get(beta-k)); }
  return(result);
}

BarrierDistribution* convolve(const BarrierDistribution* term1,
			      const BarrierDistribution* term2) {
  BarrierDistribution* result;
  result = new BarrierDistribution(term1->p,zero);
  for (int i=0; i <= footprint; i++) {
    for (int j=0; j <= (footprint-i); j++) {
      result->set(i+j,result->get(i+j) + term1->get(i)*term2->get(j) );
    }}
  /* for (int i = 0; i <= footprint; i++)
     if (result->get(i) > 1) result->set(i,1); */
  return(result);
}
