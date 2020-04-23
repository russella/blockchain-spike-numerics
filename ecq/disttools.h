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

#ifndef __DISTCLASS_H
#define __DISTCLASS_H

enum InitializationType {zero, identity};

const int maxsteps = 10000;
const int footprint = 2 * maxsteps + 1;
const int maxdelta = 30;

class Dist_index {
public:
  const int delta;
  Dist_index(int,   //  initial delta
	     int,   //  initial margin
	     int);  //  initial h_transition
  Dist_index* evolve(int,         // adversarial success
		     int) const;  // honest successes
  int get_internal() const;
  int get_transition() const;
  void set(int,   //  margin
	   int);  //  h_transition
private:
  int  beta;          // in range [-maxsteps ... maxsteps]
  int  h_transition;  // in range [0 ... delta]
};

class Distribution {
  
public:
  const int delta;
  Distribution(int,InitializationType); // delta
  Distribution(const Distribution*); // copy constructor
  //
  void show() const;
  double pdensity() const;
  double tdensity() const;
  friend Distribution* evolve(const Distribution*,
			      double,  // adversarial Poisson param
			      double); // honest prob
  
private:
  double sites[footprint][maxdelta+1];
  int    rcheck_internal(int) const;
  int    rcheck_transition(int) const;
  // accessor functions
  double get(Dist_index*) const;    // index object
  void   set(Dist_index*, double);  // index object, value
};

#endif
