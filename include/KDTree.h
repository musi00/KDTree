/********************************************************************
 * File: KDTree.h
 * Author:  Mustafa Mohamad
 *
 * An interface representing a kd-tree in some number of dimensions.
 * The tree can be constructed from a set of value and then queried
 * for membership and nearest neighbors.
 */

#ifndef KDTree_Included
#define KDTree_Included

#include "Point.h"
#include "BoundedPQueue.h"
#include <stdexcept>
#include <cmath>
#include <memory>
#include <unordered_map>
using namespace std;

template <size_t N, typename ElemType> class KDTree {
public:
	/* Constructor: KDTree();
	 * Usage: KDTree myTree;
	 * ----------------------------------------------------
	 * Constructs an empty KDTree.
	 */
	KDTree();

	/* Destructor: ~KDTree()
	 * Usage: (implicit)
	 * ----------------------------------------------------
	 * Cleans up all resources used by the KDTree.
	 */
	/* ~KDTree(); */

	/* KDTree(const KDTree& rhs);
	 * KDTree& operator= (const KDTree& rhs);
	 * Usage: KDTree one = two;
	 * Usage: one = two;
	 * -----------------------------------------------------
	 * Deep-copies the contents of another KDTree into this
	 * one.
	 */
	KDTree(const KDTree& rhs);
	KDTree& operator= (const KDTree& rhs);

	/* size_t dimension() const;
	 * Usage: size_t dim = kd.dimension();
	 * ----------------------------------------------------
	 * Returns the dimension of the points stored in this
	 * KDTree.
	 */
	size_t dimension() const;

	/* size_t size()  const;
	 * bool   empty() const;
	 * Usage: if (kd.empty())
	 * ----------------------------------------------------
	 * Returns the number of elements in the kd-tree and
	 * whether the tree is empty.
	 */
	size_t size() const;
	bool   empty() const;

	/* bool contains(const Point<N>& pt) const;
	 * Usage: if (kd.contains(pt)) { ... }
	 * ----------------------------------------------------
	 * Returns whether the specified point is contained in 
	 * the KDTree.
	 */
	bool contains(const Point<N>& pt) const;

	/* void insert(const Point<N>& pt, const ElemType& value);
	 * Usage: kd.insert(v, "This value is associated with v.");
	 * ----------------------------------------------------
	 * Inserts the point pt into the KDTree, associating it
	 * with the specified value.  If the element already existed
	 * in the tree, the new value will overwrite the existing
	 * one.
	 */
	void insert(const Point<N>& pt, const ElemType& value);

	/* ElemType& operator[] (const Point<N>& pt);
	 * Usage: kd[v] = "Some Value";
	 * ----------------------------------------------------
	 * Returns a reference to the value associated with point
	 * pt in the KDTree.  If the point does not exist, then
	 * it is added to the KDTree using the default value of
	 * ElemType as its key.
	 */
	ElemType& operator[] (const Point<N>& pt);

	/* ElemType& at(const Point<N>& pt);
	 * const ElemType& at(const Point<N>& pt) const;
	 * Usage: cout << kd.at(v) << endl;
	 * ----------------------------------------------------
	 * Returns a reference to the key associated with the point
	 * pt.  If the point is not in the tree, this function throws
	 * an out_of_range exception.
	 */
	ElemType& at(const Point<N>& pt);
	const ElemType& at(const Point<N>& pt) const;

	/* ElemType kNNValue(const Point<N>& key, size_t k) const
	 * Usage: cout << kd.kNNValue(v, 3) << endl;
	 * ----------------------------------------------------
	 * Given a point v and an integer k, finds the k points
	 * in the KDTree nearest to v and returns the most common
	 * value associated with those points.  In the event of
	 * a tie, one of the most frequent value will be chosen.
	 */
	ElemType kNNValue(const Point<N>& key, size_t k) const;

private:
  /* TODO: Add implementation details here. */
  
  enum ChildIndicator {UNSET, LEFT, RIGHT};

  typedef struct KDTreeNode {
  public:
    KDTreeNode() : value() {}

    /* Copy constructor */
    KDTreeNode (const KDTreeNode& other) {
      copyOther(other);
    }

    /* assignment operator */
    KDTreeNode& operator=(const KDTreeNode& other) {
      /* destroy old children */
      left.reset(nullptr);
      right.reset(nullptr);

      if (this != &other) {
        copyOther(other);
      }

      return *this;
    }

    void copyOther (const KDTreeNode& other) {
      point = other.point;
      value = other.value;

      /* Copy children if they exist */
      if (other.left) {
        left.reset(new KDTreeNode);
        *left = *(other.left);
      }
      if (other.right) {
        right.reset(new KDTreeNode);
        *right = *(other.right);
      }
    }

    /* Members */
    // TODO: Add underscore to all memebers and everywhere they are used
    Point<N> point;
    ElemType value;
    unique_ptr<KDTreeNode> left;
    unique_ptr<KDTreeNode> right;
  } KDTreeNode;

  void copyOther (const KDTree& rhs);
  KDTreeNode* findNode(const Point<N>& pt, KDTreeNode** parent_ptr = nullptr, ChildIndicator* which_child = nullptr) const;
  
  KDTreeNode* addNewNode (const Point<N>& pt, const ElemType& value, KDTreeNode* parent, const ChildIndicator& which_child);
  
  void kNNValueHelper(const Point<N>& key, KDTreeNode* curr, BoundedPQueue<ElemType>* bounded_priority_q, int cmp_index) const;
 
  ElemType mostFrequentElement(BoundedPQueue<ElemType>* bounded_priority_q) const;

  /* Members */
 
  unique_ptr<KDTreeNode> root_;
 
  /* Dimension of points used as keys for the tree nodes */
  size_t dim_;

  /* Number of nodes in the tree */
  size_t size_;

};

/* * * * * Implementation Below This Point * * * * */

template <size_t N, typename ElemType> KDTree<N, ElemType>::KDTree() {
  dim_ = N;
  size_ = 0;
}

/*template <size_t N, typename ElemType> KDTree<N, ElemType>::~KDTree() {
  // TODO: Fill this in.
}*/

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree& rhs) {
 copyOther(rhs); 
}

template <size_t N, typename ElemType>
typename KDTree<N, ElemType>::KDTree& KDTree<N, ElemType>::operator= (const KDTree& rhs) {
  root_.reset(nullptr);
  if (this != &rhs) {
    copyOther(rhs);
  }
  return *this;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::copyOther (const KDTree& rhs) {
  dim_ = rhs.dim_;
  size_ = rhs.size_;
  if (rhs.root_) {
    root_.reset(new KDTreeNode);
    /* use the KDTreeNode asssignment operator */
    *root_ = *rhs.root_;
  }
}


template <size_t N, typename ElemType> size_t KDTree<N, ElemType>::dimension() const {
  return dim_;
}


template <size_t N, typename ElemType> 
size_t KDTree<N, ElemType>::size() const {
  return size_;
}

template <size_t N, typename ElemType> bool KDTree<N, ElemType>::empty() const {
  return (size_ == 0);
}

template <size_t N, typename ElemType> 
bool KDTree<N, ElemType>::contains(const Point<N>& pt) const {
  return (findNode(pt) != nullptr);
}


template <size_t N, typename ElemType> 
void KDTree<N, ElemType>::insert(const Point<N>& pt, const ElemType& value) {
  KDTreeNode* parent;
  KDTreeNode** parent_ptr = &parent;
  ChildIndicator which_child = UNSET;
  KDTreeNode* node = findNode(pt, parent_ptr, &which_child);
  if (!node) {
    addNewNode(pt, value, *parent_ptr, which_child); 
  }
  else {
    node->value = value;
  }
}



template <size_t N, typename ElemType> 
ElemType& KDTree<N, ElemType>::operator[] (const Point<N>& pt) {
  KDTreeNode* parent;
  KDTreeNode** parent_ptr = &parent;
  ChildIndicator which_child = UNSET;
  KDTreeNode* node = findNode(pt, parent_ptr, &which_child);
  if (!node) {
    ElemType value = ElemType();
    node = addNewNode(pt, value, *parent_ptr, which_child); 
  }
  return node->value;
}

template <size_t N, typename ElemType> 
ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) {
  KDTreeNode* node = findNode(pt);
  if (!node)
    throw out_of_range("Point not in k-d tree");
  else
    return node->value;
}

template <size_t N, typename ElemType> 
const ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) const {
  KDTreeNode* node = findNode(pt);
  if (!node)
    throw out_of_range("Point not in k-d tree");
  else
    return node->value;
}

template <size_t N, typename ElemType> 
ElemType KDTree<N, ElemType>::kNNValue(const Point<N>& key, size_t k) const {
  if (empty()) {
    throw out_of_range("Can't perform k-nearest neighbour search on empty tree");
  }

  KDTreeNode* curr = root_.get();
  BoundedPQueue<ElemType> bounded_priority_q(k);
  int starting_cmp_index = 0;
  kNNValueHelper(key, curr, &bounded_priority_q, starting_cmp_index);

  /* Find the most frequent element in the queue */
  return mostFrequentElement(&bounded_priority_q);
}


template <size_t N, typename ElemType> 
typename KDTree<N, ElemType>::KDTreeNode* KDTree<N, ElemType>::findNode (const Point<N>& pt,
                                                             KDTreeNode** parent_ptr,
                                                             ChildIndicator* which_child) const {
  KDTreeNode* curr = root_.get();
  if (parent_ptr)
    *parent_ptr = nullptr;
  ChildIndicator tmp = UNSET;
  if (!which_child)
    which_child = &tmp;  
  int cmp_index = 0;

  while (curr != nullptr && curr->point != pt) {
   if (parent_ptr)
    *parent_ptr = curr;

   if ( pt[cmp_index] <= curr->point[cmp_index]) {
      curr = curr->left.get();
      *which_child = LEFT;
    }
    else {
      curr = curr->right.get();
      *which_child = RIGHT;
    }
    cmp_index = (cmp_index + 1) % N;
  }

  return curr;
}


template <size_t N, typename ElemType> 
typename KDTree<N, ElemType>::KDTreeNode* KDTree<N, ElemType>::addNewNode (const Point<N>& pt, const ElemType& value, KDTreeNode* parent, const ChildIndicator& which_child) {
  unique_ptr<KDTreeNode> new_node (new KDTreeNode);
  new_node->point = pt;
  new_node->value = value;
  KDTreeNode* new_node_raw_ptr = new_node.get();
  if(!parent) {
    root_ = move(new_node);
  }
  else {
    if (which_child == LEFT)
      parent->left = move(new_node);
    else
      parent->right = move(new_node);
  }
  ++size_;
  return new_node_raw_ptr;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::kNNValueHelper(const Point<N>& key, KDTreeNode* curr, BoundedPQueue<ElemType>* bounded_priority_q, int cmp_index) const {
  if (!curr)
    return;
  
  /* distance from node's point to the key point is the priority in the queue*/
  double priority = Distance(curr->point, key);
  bounded_priority_q->enqueue(curr->value, priority);

  /* Recursively search the tree for the test point (key) */
  if (key[cmp_index] <= curr->point[cmp_index])
    kNNValueHelper(key, curr->left.get(), bounded_priority_q, (cmp_index + 1) % N);
  else
    kNNValueHelper(key, curr->right.get(), bounded_priority_q, (cmp_index + 1) % N);

  /* if the queue isn't full or the current direction of the candidate hypersphere 
   * is less than the distance of the farthest k-neighrest neighbour then examine
   * the other subtree */
  if (!bounded_priority_q->full() || fabs(curr->point[cmp_index] - key[cmp_index]) < bounded_priority_q->worst()) {
    /* recursivly search the other subtree */
    if (key[cmp_index] <= curr->point[cmp_index])
      kNNValueHelper(key, curr->right.get(), bounded_priority_q, (cmp_index + 1) % N);
    else
      kNNValueHelper(key, curr->left.get(), bounded_priority_q, (cmp_index + 1) % N);
  }

}


template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::mostFrequentElement(BoundedPQueue<ElemType>* bounded_priority_q) const {
  unordered_map<ElemType,int> elem_count;

  /* Create a table of element frequency using a hash table.  When the same element
   * is encountered again it hashes to the same cell causing the count value to go
   * up by one */
  while(!bounded_priority_q->empty()) {
    ++elem_count[bounded_priority_q->top().second];
    bounded_priority_q->pop();
  }

  auto cmp = [] (const pair<ElemType,int>& a, const pair<ElemType,int>& b)
    {
      return a.second < b.second;
    };

  return max_element(elem_count.begin(), elem_count.end(), cmp)->first;
}


 




#endif
