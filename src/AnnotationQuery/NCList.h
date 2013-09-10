#ifndef _NC_LIST_H_
#define _NC_LIST_H_

#include <string>
#include <list>
#include "base/SVector.h"

//======================================================
/** 
  * Each sublist contains a set of contiguous intervals.
  */
template<class IntervalType>
class Sublist
{
public:
  Sublist(): intervals() {}
  /** ctor1 used for initial memory allocation */
  Sublist(int initRes): intervals() { reserve(initRes); }

  ~Sublist() {}

  /** Add a new interval to the list */
  void addInterval(IntervalType* i) {
    intervals.push_back(i);
  }

  /** Reserve memory for the interval vector */
  void reserve(int numOfItems) { intervals.reserve(numOfItems); }

  /** Checks to see if the latest interval contains the given interval */
  bool contains(const IntervalType& other) const {
    if(intervals.size()>0) { return intervals.back()->contains(other); }
    else { return false; }
  }
  
 /** Use to check if the latest interval has an associated sublist */
 bool hasSublist() const {
    if(intervals.size()>0) { return intervals.back()->hasSublist(); }
    else { return false; }
 }
   
 /** Use to get the associated sublist of the latest interval if it has one */
 int getSublist() const {
    if(intervals.size()>0) { return intervals.back()->getSublist(); }
    else { return -1; }
 }
 
 /** Use to set the given sublist index for the latest interval */
 void setSublist(int sIndex) {
    if(intervals.size()>0) { intervals.back()->setSublist(sIndex); }
 }

 /** Get all the intervals in the list that have any overlap with the given subject */
 int getAnyOverlaps(IntervalType* subject, svec<IntervalType*>& results) const; 
   
private:
  svec<IntervalType*> intervals;
};



//======================================================
/** NCList - Nested Containment Lists
 * Based on http://bioinformatics.oxfordjournals.org/content/23/11/1386
 */
template<class IntervalType>
class NCList {
public:
  NCList():sublists() {}

  void clear() { sublists.clear(); }

  /** Given a vector of IntervalType items, this function constructs the relevant sublists 
   *  Note that you should make sure that the sublists in the input IntervalType objects are all unset. 
 TODO decouple sublists from underlying input ds by adding them into the NCList class instead.
   */
  void constructSublists(const svec<IntervalType*>& input); 
  int getAnyOverlapping(IntervalType* subject, svec<IntervalType*>& results) const; 

private:
  /** Recursion for finding all overlaps in the nested containment lists */
  void getOverlapsFromSublist(IntervalType* subject, int sublistIndex, svec<IntervalType*>& results) const;

  /** 
   * Adds a new sublist to the records (i.e sublists vector) containing the given interval.
   * initRes tells the function how much initial memory to reserve in terms of number of items
   * Returns the index to the added sublist item in the sublists vector
   * so that it can be used for links to transcripts
   */
  int addNewSublist(IntervalType* interval, int initRes);

  svec< Sublist<IntervalType> > sublists; /// Container for holding all the nested sublists
};


//======================================================
template<class IntervalType>
int Sublist<IntervalType>::getAnyOverlaps(IntervalType* subject, svec<IntervalType*>& results) const {
  typename svec<IntervalType*>::const_iterator lBound = lower_bound(intervals.begin(), intervals.end(),
                                                                subject, IntervalType());
  // Scan forward from the lBound to find the overlapping items 
  for(typename svec<IntervalType*>::const_iterator it = lBound; it != intervals.end() ; ++it) {
    if(subject->getCoords().hasOverlap((*it)->getCoords())) { results.push_back(*it); }
    else { break; } // Gone beyond any possible overlap
  }
  // Scan backward from the lBound to find the overlapping items
  typename svec<IntervalType*>::const_reverse_iterator rit(lBound);
  for(; rit != intervals.rend() ; ++rit) {
    if(subject->getCoords().hasOverlap((*rit)->getCoords())) { results.push_back(*rit); }
    else { break; } // Gone beyond any possible overlap
  }
  return results.isize();   
}
   

//======================================================
template<class IntervalType>
void NCList<IntervalType>::constructSublists(const svec<IntervalType*>& input) 
{ 
  if(input.empty()) { return; }
  sublists.clear();                  // Make sure sublists don't exist from a previous creation 
  sublists.reserve(input.isize()/2); // Rough estimate to avoid fragmented memory allocation
  std::list<int> stack;              // Depth first stack for setting up the sublists
  typename svec<IntervalType*>::const_iterator it = input.begin();
  int slIndex = addNewSublist(*it, input.isize()/2);
  stack.push_front(slIndex);
  for(it = it + 1; it != input.end() ; ++it) {
    while(!stack.empty()) {
      if(sublists[stack.front()].contains(**it)) {
        if(sublists[stack.front()].hasSublist()) {
          sublists[sublists[stack.front()].getSublist()].addInterval(*it);
          stack.push_front(sublists[stack.front()].getSublist()); 
        } else { 
          int slIndex = addNewSublist(*it, 100); 
          sublists[stack.front()].setSublist(slIndex);
          stack.push_front(slIndex);
        } 
        break;
      } else {
        stack.pop_front();
        if(stack.empty()) {// the item belongs to the top level sublist
          sublists[0].addInterval(*it);
          stack.push_front(0);
          break;
        } 
      }
    }
  }
}

template<class IntervalType>
int NCList<IntervalType>::getAnyOverlapping(IntervalType* subject, svec<IntervalType*>& results) const {
  int sublistIndex = 0;
  getOverlapsFromSublist(subject, sublistIndex, results);
  return results.isize();
}

template<class IntervalType>
void NCList<IntervalType>::getOverlapsFromSublist(IntervalType* subject, int sublistIndex, svec<IntervalType*>& results) const {
  Sublist<IntervalType> slist = sublists[sublistIndex];
  svec<IntervalType*> out;
  slist.getAnyOverlaps(subject, out);
  for(typename svec<IntervalType*>::iterator it = out.begin(); it != out.end() ; ++it) {
    results.push_back(*it);
    if((*it)->hasSublist()) { getOverlapsFromSublist(subject, (*it)->getSublist(), results); } // Recurse with the sublist
  }
  return;
}

template<class IntervalType>
int NCList<IntervalType>::addNewSublist(IntervalType* interval, int initRes) { 
  Sublist<IntervalType> sl(initRes);
  sl.addInterval(interval);
  sublists.push_back(sl);
  int slIndex = (int) sublists.size() -1; //index of added sublist
  return slIndex;
}

#endif //_NC_LIST_H_
