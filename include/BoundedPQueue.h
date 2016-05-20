/**************************************************************
 * File: BoundedPQueue.h
 * Author: Mustafa Mohamad, based on the implementation of 
 * Keith Schwarz (htiek@cs.stanford.edu)
 *
 * An implementation of the bounded priority queue abstraction.
 * A bounded priority queue is in many ways like a regular priority
 * queue.  It stores a collection of elements tagged with a real-
 * valued priority, and allows for access to the element whose
 * priority is the smallest.  However, unlike a regular priority
 * queue, the number of elements in a bounded priority queue has
 * a hard limit that is specified in the constructor.  Whenever an
 * element is added to the bounded priority queue such that the
 * size exceeds the maximum, the element with the highest priority
 * value will be ejected from the bounded priority queue.  In this
 * sense, a bounded priority queue is like a high score table for
 * a video game that stores a fixed number of elements and deletes
 * the least-important entry whenever a new value is inserted.
 *
 * When creating a bounded priority queue, you must specify the
 * maximum number of elements to store in the queue as an argument
 * to the constructor.  For example:
 *
 * BoundedPQueue<int> bpq(15); // Holds up to fifteen values.
 *
 * The maximum size of the bounded priority queue can be obtained
 * using the maxSize() function, as in
 *
 * size_t k = bpq.maxSize();
 *
 * Beyond these restrictions, the bounded priority queue behaves
 * similarly to other containers.  You can query its size using
 * size() and check whether it is empty using empty().  You
 * can push an element into the bounded priority queue by
 * writing
 *
 * bpq.push(elem, priority);
 *
 * Note that after enqueuing the element, there is no guarantee 
 * that the value will actually be in the queue.  If the queue
 * is full and the new element's priority exceeds the largest
 * priority in the container, it will not be added.
 *
 * You can dequeue elements from a bounded priority queue using
 * the dequeueMin() function, as in
 *
 * int val = bpq.dequeueMin();
 *
 * The bounded priority queue also allows you to query the min
 * and max priorities of the values in the queue.  These values
 * can be queried using the best() and worst() functions, which
 * return the smallest and largest priorities in the queue,
 * respectively.
 */

#ifndef BoundedPQueue_Included
#define BoundedPQueue_Included

#include <queue>
#include <algorithm>
#include <limits>
using namespace std;

template <typename T> 
class BoundedPQueue {
public:

	/* Typedef and compare functor of the Queue elements */
  typedef pair<double,T> QueueType;
  class QueueTypeComparator {
  public:
    bool operator() (const QueueType& a, const QueueType& b) {
      return a.first < b.first;
    }
  };

  /* Constructor: BoundedPQueue(size_t maxSize);
	 * Usage: BoundedPQueue<int> bpq(15);
	 * --------------------------------------------------
	 * Constructs a new, empty BoundedPQueue with
	 * maximum size equal to the constructor argument.
	 */
	explicit BoundedPQueue(size_t maxSize);

	/* void push(const T& value, double priority);
	 * Usage: bpq.push("Hi!", 2.71828);
	 * --------------------------------------------------
	 * pushs a new element into the BoundedPQueue with
	 * the specified priority.  If this overflows the maximum
	 * size of the queue, the element with the highest
	 * priority will be deleted from the queue.  Note that
	 * this might be the element that was just added.
	 */
	void push(const T& value, double priority);

  void pop(); 
	/* size_t size() const;
	 * bool  empty() const;
   * bool full() const;
	 * Usage: while (!bpq.empty()) { ... }
	 * --------------------------------------------------
	 * Returns the number of elements in the queue and whether
	 * the queue is emptyi or full, respectively.
	 */
	size_t size() const;
	bool empty() const;
  bool full() const;
	
	/* size_t maxSize() const;
	 * Usage: size_t queueSize = bpq.maxSize();
	 * --------------------------------------------------
	 * Returns the maximum number of elements that can be
	 * stored in the queue.
	 */
	size_t maxSize() const;
	
	/* double worst() const;
	 * Usage: double highestPriority = bpq.worst();
	 * --------------------------------------------------
	 * worst() returns the largest priority of an element
	 * stored in the container.  If an element is pushd
	 * with a priority above this value, it will automatically
	 * be deleted from the queue.  The function returns
	 * numeric_limits<double>::infinity() if the queue is
	 * empty.
	 */
	double worst() const;

  /* const QueueType& top() const;
   * Usage QueueType& top_of_queue = bpq.top() 
   * --------------------------------------------------
   * Returns a reference to top of the queue (i.e the max
   * priority element)
   */
  const QueueType& top() const;
	

private:

  /* This class is layered on top of an stl priority queue */
	priority_queue<QueueType, vector<QueueType>, QueueTypeComparator> elems;

  /* Maximum number of elements the queue can hold */
	size_t maximumSize;
};

/* * * * * Implementation Below This Point * * * * */

template <typename T> BoundedPQueue<T>::BoundedPQueue(size_t maxSize) {
	maximumSize = maxSize;
}

template <typename T> void 
BoundedPQueue<T>::push(const T& value, double priority) {
	/* Add the element to the collection. */
	elems.emplace(priority, value);

	/* If there are too many elements in the queue, remove the highest priority element */
	if (size() > maxSize())
	{
          elems.pop();
	}
}

template <typename T> 
void BoundedPQueue<T>::pop() {
  if (size() > 0) {
    elems.pop();
  }
}

template <typename T> const typename BoundedPQueue<T>::QueueType& 
BoundedPQueue<T>::top() const {
  return elems.top();
}

template <typename T> size_t 
BoundedPQueue<T>::size() const {
	return elems.size();
}

template <typename T> bool 
BoundedPQueue<T>::empty() const {
	return elems.empty();
}

template <typename T> bool 
BoundedPQueue<T>::full() const {
	return size() == maxSize();
}

template <typename T> size_t 
BoundedPQueue<T>::maxSize() const {
	return maximumSize;
}

template <typename T> double 
BoundedPQueue<T>::worst() const
{
	return empty()? numeric_limits<double>::infinity() : elems.top().first;
}


#endif
