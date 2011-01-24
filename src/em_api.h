#ifndef __EM_API_H__
#define __EM_API_H__

#include <iostream>
#include <vector>
#include <map>
#include <cstdio>

extern "C" {
#include "qsopt.h"
};

// API for EarthMover distance calculator
// maintains EarthMover distance from "base" to "current" weight vectors

// Label is the type of the node labels.  It needs to be copyable,
// default constructable, and comparable (or hashable) (see Map below)

// Map is a container class that can map Label to int.
// std::map<Label,int> and __gnu_cxx::hash_map<Label,int> are obvious
// candidates.  Of course, Label must be valid as a key for the map
// (that is, be comparable for std::map, and have a hash function for
// __gnu_cxx::hash_map.

// class EarthMover assumes that there is a function
//   double dist(Label b, Label c)
// that computes the distance between base point "b" and current point
// "c".

template <class Label,
	  class Map = std::map<Label,int> >
class EarthMover {

 public:
  // main constructor.  baseIn gives the labels of the base distribution
  // currentIn gives the labels of the current distribution
  // weightsIn gives the weights for each label.
  //   The same weight is assigned to both distributions.
  // baseIn, currentIn, and weightsIn are required to be the same length.
  // for these constructors, vector could be replaced by list, or
  // an iterator range.

  EarthMover (const std::vector<Label>& baseIn,
	      const std::vector<Label>& currentIn,
	      const std::vector<int>& weightsIn);

  // alternate constructor.
  // This constructor gives different base and current, as well as a
  // set of "good" edges.
  //
  // the size of baseIn does not need to equal the size of currentIn, and
  // the sum of weights in baseIn does not need to equal the sum of weights
  // in currentIn
  EarthMover (const std::vector<std::pair<Label,int> >& baseIn,
	      const std::vector<std::pair<Label,int> >& currentIn,
	      const std::vector<std::pair<Label,Label> >& goodEdgesIn);

  // third constructor.
  // This constructor gives different base and current, but doesn't
  // give a set of "good" edges.
  //
  // the size of baseIn does not need to equal the size of currentIn, and
  // the sum of weights in baseIn does not need to equal the sum of weights
  // in currentIn.
  //
  // this constructor is not efficient - it will use the complete graph as
  // the "good" edges, and is only provided as a preliminary, trial version
  EarthMover (const std::vector<std::pair<Label,int> >& baseIn,
	      const std::vector<std::pair<Label,int> >& currentIn);

  // destructor
  ~EarthMover ();

  // shift "weight" from "from" to "to".  "from" and "to" are also assumed
  // to possibly change location.
  void shiftWeight (Label from, Label to, int weight = 1);

  // splits "old" into two nodes "old" and "new", moving newweight
  // from "old" to "new".  This will call the distance function for "new",
  // so "new" should already have any necessary information for distance.
  void split (Label oldl, Label newl, int newweight);

  // merges "old" into "new", moving any weight on "old" to "new"
  // (so "old" ceases to exist)
  void merge (Label oldl, Label newl);

  // returns the current distance from base to current.
  double distance ();
  
  // returns true if the current distance is <= thresh
  // this can be much more efficient than (distance() <= thresh)
  // because it can use a bound instead of the actual distance.
  bool distance_le (double thresh);
  
  // returns true if the current distance is < thresh
  // I wasn't sure which would be more appropriate.
  bool distance_lt (double thresh);

 private:
 struct node {
   int weight;
   Label lbl;
   int lpnum;
 };

 QSprob lp;
 Map label2base;
 Map label2current;
 std::vector<node> base;
 std::vector<node> current;
 int totbase;
 int totcurrent;

 void build_lp_complete ();
 
};

template <class Label,
	  class Map>
EarthMover<Label,Map>::EarthMover (const std::vector<Label>& baseIn,
				   const std::vector<Label>& currentIn,
				   const std::vector<int>& weightsIn)
{
  base.resize (baseIn.size());
  totbase = 0;
  for (size_t i = 0; i < baseIn.size(); ++i) {
    label2base[baseIn[i]] = i;
    base[i].lbl = baseIn[i];
    base[i].weight = weightsIn[i];
    totbase += base[i].weight;
  }
  current.resize (currentIn.size());
  totcurrent = 0;
  for (size_t i = 0; i < currentIn.size(); ++i) {
    label2current[currentIn[i]] = i;
    current[i].lbl = currentIn[i];
    current[i].weight = weightsIn[i];
    totcurrent += current[i].weight;
  }

  //  build_lp_matching ();
  build_lp_complete ();
}

template <class Label,
	  class Map>
EarthMover<Label,Map>::EarthMover (const std::vector<std::pair<Label,int> >& baseIn,
			const std::vector<std::pair<Label,int> >& currentIn,
			const std::vector<std::pair<Label,Label> >& goodEdgesIn)
{
  lp = QScreate_prob ("earthmover", QS_MIN);
  base.resize (baseIn.size());
  totbase = 0;
  for (size_t i = 0; i < baseIn.size(); ++i) {
    label2base[baseIn[i].first] = i;
    base[i].lbl = baseIn[i].first;
    base[i].weight = baseIn[i].second;
    totbase += base[i].weight;
  }
  current.resize (currentIn.size());
  totcurrent = 0;
  for (size_t i = 0; i < currentIn.size(); ++i) {
    label2current[currentIn[i].first] = i;
    current[i].lbl = currentIn[i].first;
    current[i].weight = currentIn[i].second;
    totcurrent += current[i].weight;
  }

  //  build_lp_edges (goodEdgesIn);
  build_lp_complete ();
}

template <class Label,
	  class Map>
EarthMover<Label,Map>::EarthMover (const std::vector<std::pair<Label,int> >& baseIn,
			const std::vector<std::pair<Label,int> >& currentIn)
{
  lp = QScreate_prob ("earthmover", QS_MIN);
  base.resize (baseIn.size());
  totbase = 0;
  for (size_t i = 0; i < baseIn.size(); ++i) {
    label2base[baseIn[i].first] = i;
    base[i].lbl = baseIn[i].first;
    base[i].weight = baseIn[i].second;
    totbase += base[i].weight;
  }
  current.resize (currentIn.size());
  totcurrent = 0;
  for (size_t i = 0; i < currentIn.size(); ++i) {
    label2current[currentIn[i].first] = i;
    current[i].lbl = currentIn[i].first;
    current[i].weight = currentIn[i].second;
    totcurrent += current[i].weight;
  }

  build_lp_complete ();
}

template <class Label,
	  class Map>
EarthMover<Label,Map>::~EarthMover ()
{
  QSfree_prob (lp);
}

template <class Label,
	  class Map>
void EarthMover<Label,Map>::build_lp_complete ()
{
  extern double dist (Label a, Label b);

  int rval;
  lp = QScreate_prob ("earthmover", QS_MIN);
  for (size_t i = 0; i<base.size(); ++i) {
    char nm[32];
    sprintf (nm, "base_%ld", (long int) i);
    rval = QSnew_row (lp, double(base[i].weight), 'E', nm);
    if (rval) std::cerr << "QSnew_row " << nm << " failed: " << rval << "\n";
  }
  double cscale = -double(totbase)/double(totcurrent);
  for (size_t i = 0; i<current.size()-1; ++i) {
    char nm[32];
    sprintf (nm, "current_%ld", (long int) i);
    rval = QSnew_row (lp, current[i].weight * cscale, 'E', nm);
    if (rval) std::cerr << "QSnew_row " << nm << " failed: " << rval << "\n";
  }

  for (size_t i = 0; i<base.size(); ++i) {
    for (size_t j = 0; j < current.size(); ++j) {
      int cnt;
      int matind[2];
      double matval[2];
      char nm[45];
      double obj = dist (base[i].lbl, current[j].lbl) / totbase;
      matind[0] = i;
      matval[0] = 1;
      if (j == current.size() - 1) {
	cnt = 1;
      } else {
	cnt = 2;
	matind[1] = j + base.size();
	matval[1] = -1;
      }
      sprintf (nm, "e_%ld_%ld", (long int) i, (long int) j);
      rval = QSadd_col (lp, cnt, matind, matval, obj, 0.0, QS_MAXDOUBLE, nm);
      if (rval) std::cerr << "QSadd_col " << nm << " failed: " << rval << "\n";
    }
  }
}

template <class Label,
	  class Map>
double EarthMover<Label,Map>::distance ()
{
  int status = 0;
  int rval = QSopt_primal (lp, &status);
  if (rval) std::cerr << "QSopt_primal failed: " << rval << "\n";
  if (status != QS_LP_OPTIMAL) std::cerr << "QSopt_primal status: "
					 << status << "\n";
  double objval;
  rval = QSget_objval (lp, &objval);
  if (rval) std::cerr << "QSget_objval failed: " << rval << "\n";
  return objval;
}

#endif /* __EM_API_H__ */
