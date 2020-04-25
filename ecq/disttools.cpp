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

// The Distribution and Dist_index classes are implemented below.
//
// The Distribution class is intended to maintain a distribution over
// pairs of the form (M,T), where M is a "margin" and T is a
// "transition point." The intention is that these are quantities
// implicitly determined by a string read in thus far: M is equal to
// (A - D(H)), the difference between the number of adversarial slots
// A and the Delta height of the honest slots D(H) observed in the
// string; T is the "transition point" of the string observed so
// far--this is the length(the suffix of the string that has maximum
// height)-1. Thus, if the most recent character increased the Delta
// height, T=0.
//
// This second parameter is sufficient to determine if a new honest
// slot added to the end of the string will increase the Delta height.

//class Dist_index

Dist_index::Dist_index(int  init_delta,
		       int  init_margin,
		       int  init_transition) : delta(init_delta) {
  if (init_delta > maxdelta) throw std::invalid_argument("INDEX CONSTRUCTOR: Delta index out of range");
  if (abs(init_margin) > maxsteps) throw std::invalid_argument("INDEX CONSTRUCTOR: margin index out of range");
  if (init_transition > delta) throw std::invalid_argument("INDEX CONSTRUCTOR: transition index out of range");
  set(init_margin,init_transition);
}

void Dist_index::set(int   set_margin,
		     int   set_transition) {
  beta         = set_margin;
  h_transition = set_transition;
}

int Dist_index::get_internal() const {
  return(beta + maxsteps);
}

int Dist_index::get_transition() const {
  return(h_transition);
}

Dist_index* Dist_index::evolve(int adversarial,
			       int honest) const {
  int  result_beta;
  int  result_h_transition;
  
  result_beta = beta + adversarial;
  if (h_transition < delta-1) {
    result_h_transition = h_transition+1; }
  else
    if (honest >= 1) {
      result_h_transition = 0;
      result_beta = result_beta - 1; }
    else
      result_h_transition = h_transition;
  return(new Dist_index(delta,
			result_beta,
			result_h_transition));
}

// End of distribution index object.

//class Distribution;


void Distribution::show() const {
  cout << "Distribution contents...\n";
  for (int beta=-15; beta < 15; beta++) {
    cout << beta << " : ";
    for (int offset=0; offset <= delta; offset++)
      cout << sites[beta + maxsteps][offset] << " ";
    cout << "\n";
  }
}

int Distribution::rcheck_internal(int internal) const {
  if ((internal < 0) || (internal > 2*maxsteps)) throw std::invalid_argument("ARRAY ACCESS: internal distribution index out of range");
  return internal; }

int Distribution::rcheck_transition(int transition) const {
  if ((transition < 0) || (transition > delta)) throw std::invalid_argument("ARRAY ACCESS: transition distribution index out of range");
  return(transition); }

double Distribution::get(Dist_index* ind) const {
  return(sites[rcheck_internal(ind->get_internal())][rcheck_transition(ind->get_transition())]);
}

void Distribution::set(Dist_index* ind, double value) {
  sites[rcheck_internal(ind->get_internal())][rcheck_transition(ind->get_transition())] = value;
}

double Distribution::pdensity() const {
  double result = 0.0;
  Dist_index* index;
  index = new Dist_index(delta,0,0);
  for(int beta = 0; beta <= maxsteps; beta++)
    for (int transition = 0; transition <= delta; transition++) {
      index->set(beta,transition);
      result = result + get(index);
    }
  delete(index);
  return(result); 
}

double Distribution::tdensity() const {
  double result = 0.0;
  Dist_index* index;
  index = new Dist_index(delta,0,0);
  for(int beta = -maxsteps; beta <= maxsteps; beta++)
    for (int transition = 0; transition <= delta; transition++) {
      index->set(beta,transition);
      result = result + get(index);
    }
  delete(index);
  return(result); 
}

// Initial constructor, "structure" variable determines if zero or distribution at 0.
Distribution::Distribution(int init_delta,InitializationType structure) : delta(init_delta) {
  Dist_index* index;
  if (init_delta > maxdelta) throw std::invalid_argument("DIST CONSTRUCTOR: Delta index out of range");
  index = new Dist_index(delta,0,0);
  for(int beta = -maxsteps; beta <= maxsteps; beta++)
    for(int transition = 0; transition <= delta; transition++) {
      index->set(beta,transition);
      set(index,0.0);
    }
  if (structure==identity) {
    index->set(0, 0);
    set(index,1.0); }
  delete(index);
}

// Copy constructor
Distribution::Distribution(const Distribution* orig) : delta(orig->delta) {
  Dist_index* index;
  index = new Dist_index(delta,0,0);
  for(int beta = -maxsteps; beta <= maxsteps; beta++)
    for(int transition = 0; transition <= delta; transition++) {
      index->set(beta,transition);
      set(index,orig->get(index)); }
  delete(index);
}

Distribution* evolve(const Distribution* source,
		     double adv_prob,
		     double hon_prob) {
  double h_transitions[2];
  double a_transitions[2];
  Distribution* result;
  Dist_index* source_index;
  Dist_index* target_index;
  
  h_transitions[0] = 1- hon_prob;
  h_transitions[1] = hon_prob;
  a_transitions[0] = 1- adv_prob;
  a_transitions[1] = adv_prob;
  source_index = new Dist_index(source->delta,0,0);
  result = new Distribution(source->delta,zero);
  for (int hon : {0, 1}) // Honest move
    for (int adv : {0, 1})  // Adversarial move
      for (int beta=-(maxsteps-hon); beta <= (maxsteps-adv); beta++)
	for (int transition=0; transition <= source->delta; transition++) {
	  source_index->set(beta,transition);
	  target_index = source_index->evolve(adv,hon);
	  result->set(target_index,result->get(target_index)
		      + (source->get(source_index)
			 * a_transitions[adv]
			 * h_transitions[hon])); 
	  delete(target_index);	};
  delete(source_index);
  return(result);
}

