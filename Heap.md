# Heap Data Structure in Java

## Overview
A heap is a complete binary tree that satisfies the heap property:
- **Min Heap**: Parent nodes have values less than or equal to their children (ascending order)
- **Max Heap**: Parent nodes have values greater than or equal to their children (descending order)

Heaps are commonly implemented as priority queues in Java.

## Implementation

### Using PriorityQueue
Java provides `PriorityQueue` which implements a min-heap by default:

```java
// Min Heap (default)
Queue<Integer> minHeap = new PriorityQueue<>();

// Max Heap (using custom comparator)
Queue<Integer> maxHeap = new PriorityQueue<>((a, b) -> b - a);
```

### Common Operations
```java
// Insert element
minHeap.offer(5);  // O(log n)

// Remove and return top element
int top = minHeap.poll();  // O(log n)

// Peek top element without removing
int top = minHeap.peek();  // O(1)

// Get size
int size = minHeap.size();  // O(1)

// Check if empty
boolean isEmpty = minHeap.isEmpty();  // O(1)
```

### Array-Based Implementation
For a custom heap implementation using arrays:

```java
class MaxHeap {
    private List<Integer> heap;
    
    public MaxHeap() {
        this.heap = new ArrayList<>();
    }
    
    // Helper methods
    private int left(int i) { return 2 * i + 1; }
    private int right(int i) { return 2 * i + 2; }
    private int parent(int i) { return (i - 1) / 2; }
    
    // Peek top element
    public int peek() {
        if (heap.isEmpty()) throw new NoSuchElementException();
        return heap.get(0);
    }
    
    // Insert element
    public void push(int val) {
        heap.add(val);
        siftUp(heap.size() - 1);
    }
    
    private void siftUp(int i) {
        while (i > 0) {
            int p = parent(i);
            if (heap.get(i) <= heap.get(p)) break;
            Collections.swap(heap, i, p);
            i = p;
        }
    }
    
    // Remove top element
    public int pop() {
        if (isEmpty()) throw new NoSuchElementException();
        Collections.swap(heap, 0, heap.size() - 1);
        int val = heap.remove(heap.size() - 1);
        siftDown(0);
        return val;
    }
    
    private void siftDown(int i) {
        while (true) {
            int l = left(i);
            int r = right(i);
            int largest = i;
            
            if (l < heap.size() && heap.get(l) > heap.get(largest))
                largest = l;
            if (r < heap.size() && heap.get(r) > heap.get(largest))
                largest = r;
            if (largest == i) break;
                
            Collections.swap(heap, i, largest);
            i = largest;
        }
    }
}
```

## Time Complexities
| Operation       | Time Complexity |
|----------------|-----------------|
| Insert (push)  | O(log n)        |
| Remove (pop)   | O(log n)        |
| Peek           | O(1)            |
| Heapify        | O(n)            |

## Applications
1. **Priority Queues**: Natural implementation for priority queues
2. **Heap Sort**: O(n log n) sorting algorithm
3. **Top-K Elements**: Efficiently find K largest/smallest elements
4. **Graph Algorithms**: Dijkstra's, Prim's algorithms

## Key Points
- Complete binary tree property makes array representation efficient
- Heapify operations maintain heap property after insertions/removals
- Java's `PriorityQueue` is a good default choice
- Custom implementations allow more control when needed

## Example Usage
```java
// Find top 3 elements
Queue<Integer> minHeap = new PriorityQueue<>();
for (int num : nums) {
    minHeap.offer(num);
    if (minHeap.size() > 3) minHeap.poll();
}
// minHeap now contains top 3 elements
```
