ArrayQueue::ArrayQueue() {
  QueueSize = 0;
  first     = NULL;
  last      = NULL:
}
 

// Add another entry (which is an double array) to the queue:
void ArrayQueue::put(int size) {
  
  // If it is the first, set its address as first and last:
  if (QueueSize==0) {
    first = (QueueEntry *) new QueueEntry;
    last  = first;
  }
  // If not, put it behind the last and define it as the new last one:
  else {
    last->next  = (QueueEntry *) new QueueEntry;
    last        = last->next;
  }
  
  // Set properties of the new entry (which is the last):
  last->next  = NULL;
  last->size  = size;
  last->Array = vector<double>(0, size-1);
  // Increase queue size:
  QueueSize++;
}

// Remove the first entry from the queue: 
int ArrayQueue::take() {
  
  

}

