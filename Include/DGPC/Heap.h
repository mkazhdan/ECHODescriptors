/************************************************************
 * This file is part of the DGPC library. The library computes
 * Discrete Geodesic Polar Coordinates on a polygonal mesh.
 *
 * More info:
 *   http://folk.uio.no/eivindlm/dgpc/
 *
 * Authors: Eivind Lyche Melvær and Martin Reimers
 * Centre of Mathematics and Department of Informatics
 * University of Oslo, Norway, 2012
 ************************************************************/
#ifndef DGPC_HEAP_H
#define DGPC_HEAP_H

#include <queue>
#include <functional>

namespace DGPC {

template<class real>
class HeapNode {

public:
  int   idx_;
  real  key_;

  HeapNode( int idx, real key) { idx_ = idx; key_ = key;}
  ~HeapNode(){}

  bool operator >  ( const HeapNode<real>& x) const { return (this->key_ >  x.key_);}
  bool operator >= ( const HeapNode<real>& x) const { return (this->key_ >= x.key_);}
  bool operator <  ( const HeapNode<real>& x) const { return (this->key_ <  x.key_);}
  bool operator <= ( const HeapNode<real>& x) const { return (this->key_ <= x.key_);}
};

template<class real>
class Heap {

  typedef std::priority_queue< HeapNode<real>, std::vector< HeapNode<real> >, std::greater<HeapNode<real> > > Heap_t;

  Heap_t heap_;
  std::vector<bool>   flags_;
  std::vector<real>*  keys_;

public:
  Heap( ) { keys_ = 0; }
  Heap( std::vector<real>*    keys) { initialize(keys);}
  ~Heap() {}

  void initialize( std::vector<real>*    keys) {	keys_ = keys; resize (keys_->size()); std::fill(flags_.begin(),flags_.end(),false);}
  void resize( int size) { flags_.resize(size); }
  bool isInHeap( int idx) const { return flags_[idx];}
  void removeFromHeap( int idx) { flags_[idx]=false;}
  void push( int idx) { heap_.push(HeapNode<real>(idx,(*keys_)[idx])); flags_[idx]=true; }
  bool empty () { return (top()==-1);}
  int  size () const { return heap_.size();}

  int getCandidate() {
	if (!heap_.empty()) {
		int ret = top(); 
		pop();
		return ret;
	} else 
		return -1;
  }

  void pop() { 
    while (!heap_.empty() && !isInHeap(heap_.top().idx_)) 
      heap_.pop();
    if (!heap_.empty()) { 
      flags_[heap_.top().idx_]=false; 
      heap_.pop();
    }
  } 

  int  top() { 
    while (!heap_.empty() && !isInHeap(heap_.top().idx_)) 
      heap_.pop();
    if (!heap_.empty())
      return heap_.top().idx_;
    else
      return -1;
  } 

};

}; //End namespace DGPC

#endif // DGPC_HEAP_H
