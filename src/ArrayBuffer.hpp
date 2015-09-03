#ifndef ARRAYQUEUE_H    // include guard.
#define ARRAYQUEUE_H 1

#include <cstddef> // For NULL pointer.


struct QueueEntry {
  int size;
  double *Array;
  QueueEntry *next;
};

class ArrayQueue {
private:
  // Information storage:
  int QueueSize;
  QueueEntry *first;
  QueueEntry *last;
public:
  ArrayQueue();
  void   put(int size);
  int    take();
  void   writeLast(int i);
  double readFirst(int i);
};


#endif
