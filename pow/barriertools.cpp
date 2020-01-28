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

//class Dist_index

/* The internal variables for this class are:
   private:
     int  beta;               // current margin
     int  r_isolation;        // right isolation since last honest slot
     bool l_isolated_pending; // uncounted, left-isolated honest slot?
*/

Dist_index::Dist_index(int  init_delta,
		       int  init_margin,
		       int  init_r_isolation,
		       bool init_l_isolated_pending) : delta(init_delta) {
  if (init_delta > maxdelta) throw std::invalid_argument("Delta index out of range");
  set(init_margin,init_r_isolation,init_l_isolated_pending);
}

void Dist_index::set(int   set_margin,
		     int   set_r_isolation,
		     bool  set_l_pending) {
  beta               = set_margin;
  l_isolated_pending = set_l_pending;
  r_isolation        = set_r_isolation;
}

int Dist_index::get_beta() const {
  return(beta);
}

int Dist_index::get_internal() const {
  if (l_isolated_pending)
    return(r_isolation);
  else
    return(r_isolation+delta+1);
  // Assumes array is [0...beta] X [0...2 Delta + 2]
  // so that internal is in [0 ... 2 Delta + 1]
}

Dist_index* Dist_index::evolve(int adversarial,
			       int honest) const {
  int  beta_change_honest = 0;
  int  result_beta;
  int  result_r_isolation;
  bool result_l_isolated_pending;
  
  beta_change_honest = 0;
  if (honest >= 1) {
    result_r_isolation = 0;
    if ((r_isolation >= delta) && (honest == 1)) 
      result_l_isolated_pending = true; 
    else 
      result_l_isolated_pending = false; }
  else // no honest leader in this slot
    if (r_isolation < delta-1) {
      result_r_isolation = r_isolation+1;
      result_l_isolated_pending = l_isolated_pending; }
    else {
      result_r_isolation = delta;
      result_l_isolated_pending = false;
      if (l_isolated_pending) beta_change_honest = -1; }
  // Above computes resulting r_iso, l_iso_pending, and beta_change_honest.
  if (beta > 0) result_beta = beta + beta_change_honest + adversarial;
  else result_beta = adversarial;
  return(new Dist_index(delta,
			result_beta,
			result_r_isolation,
			result_l_isolated_pending));
}

// End of distribution index object.

//class BarrierDistribution;

double Poisson(double lambda, int k) {
  if (k < 0) throw std::invalid_argument("Poisson index out of range");
  double result = 1.0;
  for (int i=1; i <= k; i++)
    result = result * (lambda / i);
  return( exp(-lambda) * result );
}

double gaugetail(double shift, int beta) {
  if (beta <= shift)
    return(1.0);
  else
    return(pow(shift*exp(1.0)/beta,beta) * exp(-shift)); }

/* This is the spike distribution if the adversary has full ability to distribute C units of expectation among the slots. */
//  double spikedist(double shift, int beta) {
//  return(gaugetail(shift,beta) - gaugetail(shift,beta+1));
//  }


double spikedist(double shift, int beta) {
  return(Poisson(shift,beta));
}

void BarrierDistribution::show() const {
  cout << "Distribution contents...\n";
  for (int beta=0; beta < 30; beta++) {
    cout << beta << " : ";
    for (int internal=0; internal <= 2*delta+1; internal++)
      cout << sites[beta][internal] << " ";
    cout << "\n";
  }
}

int BarrierDistribution::rcheck_beta(int beta) const {
  if ((beta < 0) || (beta > maxsteps)) throw std::invalid_argument("distribution index out of range");
  return beta; }

int BarrierDistribution::rcheck_delta(int internal) const {
  if ((internal < 0) || (internal > 2 * delta + 1)) throw std::invalid_argument("distribution index out of range");
  return(internal); }

double BarrierDistribution::get(Dist_index* ind) const {
  return(sites[rcheck_beta(ind->get_beta())]
	 [rcheck_delta(ind->get_internal())]);
}

void BarrierDistribution::set(Dist_index* ind, double value) {
  sites[rcheck_beta(ind->get_beta())]
    [rcheck_delta(ind->get_internal())] = value;
}

double BarrierDistribution::pdensity() const {
  double result = 0.0;
  Dist_index* index;
  index = new Dist_index(delta,0,0,false);
  for(int beta = 1; beta <= maxsteps; beta++)
    for (int r_iso = 0; r_iso <= delta; r_iso++)
      for (bool pend : {false, true}) {
	index->set(beta,r_iso,pend);
	result = result + get(index);
      }
  delete(index);
  return(result); 
}

double BarrierDistribution::tdensity() const {
  double result = 0.0;
  Dist_index* index;
  index = new Dist_index(delta,0,0,true);
  for(int beta = 0; beta <= maxsteps; beta++)
    for (int r_iso = 0; r_iso <= delta; r_iso++)
      for (bool pend : {false, true}) {
	index->set(beta,r_iso,pend);
	result = result + get(index); }
  delete(index);
  return(result); 
}

// Initial constructor, "structure" variable determines if zero or distribution at 0.
BarrierDistribution::BarrierDistribution(int init_delta,InitializationType structure) : delta(init_delta) {
  Dist_index* index;
  if (init_delta > maxdelta) throw std::invalid_argument("Delta index out of range");
  index = new Dist_index(delta,0,0,true);
  for(int beta = 0; beta <= maxsteps; beta++)
    for(int r_iso = 0; r_iso <= delta; r_iso++)
      for(bool pend : {false, true}) {
	index->set(beta,r_iso,pend);
	set(index,0.0);
      }
  if (structure==identity) {
    index->set(0, 0, false);
    set(index,1.0); }
  delete(index);
}

// Copy constructor
BarrierDistribution::BarrierDistribution(const BarrierDistribution* orig) : delta(orig->delta) {
  Dist_index* index;
  index = new Dist_index(delta,0,0,true);
  for(int beta = 0; beta <= maxsteps; beta++)
    for(int r_iso = 0; r_iso <= delta; r_iso++)
      for(bool pend : {false, true}) {
	index->set(beta,r_iso,pend);
	set(index,orig->get(index)); }
  delete(index);
}

double stat_distance(const BarrierDistribution* dista,
		     const BarrierDistribution* distb) {
  double result = 0;
  Dist_index* index;
  index = new Dist_index(dista->delta,0,0,true);
  for (int beta=0; beta <= maxsteps; beta++) 
    for (int r_iso=0; r_iso <= dista->delta; r_iso++)
      for (bool pend : {false, true}) {
	index->set(beta,r_iso,pend);
	result = result + abs(dista->get(index)
			      -distb->get(index)); }
  delete(index);
  return(result/2);
}


BarrierDistribution* convolve_spike(const BarrierDistribution* base,
				    double spike_param) {
  BarrierDistribution* result;
  Dist_index* source_index;
  Dist_index* target_index;
  source_index = new Dist_index(base->delta,0,0,true);
  target_index = new Dist_index(base->delta,0,0,true);
  result = new BarrierDistribution(base->delta,zero);
  for (int beta_a=0; beta_a <= maxsteps; beta_a++) 
    for (int beta_b=0; beta_b <= (maxsteps-beta_a); beta_b++)
      for (int r_iso=0; r_iso <= base->delta; r_iso++)
	for (bool pend : {false, true}) {
	  source_index->set(beta_a,r_iso,pend);
	  target_index->set(beta_a+beta_b,r_iso,pend);
	  result->set(target_index,result->get(target_index)
		      + base->get(source_index)*spikedist(spike_param,beta_b));
	}
  delete(source_index);
  delete(target_index);
  return(result);
}

BarrierDistribution* evolve(const BarrierDistribution* source,
			    EvolutionType convention,
			    double adv_param,
			    double hon_param) {
  double h_transition_pr;
  BarrierDistribution* result;
  Dist_index* source_index;
  Dist_index* target_index;
  int initial_beta;
  if (convention==reflect)
    initial_beta = 0;
  else
    initial_beta = 1;
  source_index = new Dist_index(source->delta,0,0,true);
  result = new BarrierDistribution(source->delta,zero);
  for (int hon : {0, 1, 2}) {            // Honest move
    if (hon == 0) h_transition_pr = exp(-hon_param);
    else if (hon == 1) h_transition_pr = hon_param * exp(-hon_param);
    else h_transition_pr = 1 - exp(-hon_param) * (1 + hon_param);
    for (int adv=0; adv <= maxsteps; adv++)  // Adversarial move
      for (int beta=initial_beta; beta <= (maxsteps-adv); beta++)
	for (int r_iso=0; r_iso <= source->delta; r_iso++)
	  for (bool pend : {false, true}) {
	    source_index->set(beta,r_iso,pend);
	    target_index = source_index->evolve(adv,hon);
	    result->set(target_index,result->get(target_index)
			+ source->get(source_index) * Poisson(adv_param,adv) * h_transition_pr);
	    delete(target_index);
	  }};
  delete(source_index);
  return(result);
}

