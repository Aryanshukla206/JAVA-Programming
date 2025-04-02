Queue
---

### **1. Core Queue Operations (FIFO)**
| Method | Description | Throws Exception | Returns Special Value |
|--------|------------|------------------|-----------------------|
| `add(E e)` | Inserts element | `IllegalStateException` if full | `offer(E e)` |
| `remove()` | Removes and returns head | `NoSuchElementException` if empty | `poll()` |
| `element()` | Returns head (no removal) | `NoSuchElementException` if empty | `peek()` |

---

### **2. Implementation-Specific Methods**
#### **LinkedList/ArrayDeque**
| Method | Description |
|--------|------------|
| `offerFirst(E e)` | Adds to front (deque operation) |
| `offerLast(E e)` | Adds to end (standard queue) |
| `pollFirst()` | Removes from front |
| `pollLast()` | Removes from end |
| `peekFirst()` | Views front element |
| `peekLast()` | Views end element |

#### **PriorityQueue**
| Method | Description |
|--------|------------|
| `comparator()` | Returns ordering comparator |
| `spliterator()` | Creates late-binding spliterator |

---

### **3. Bulk Operations**
| Method | Description |
|--------|------------|
| `addAll(Collection<? extends E> c)` | Adds all elements |
| `clear()` | Removes all elements |
| `contains(Object o)` | Checks for element |

---

### **4. Java 8+ Features**
| Method | Description |
|--------|------------|
| `removeIf(Predicate<? super E> filter)` | Conditionally removes |
| `forEach(Consumer<? super E> action)` | Iterates elements |
| `stream()` | Creates sequential stream |
| `parallelStream()` | Creates parallel stream |

---

### **5. Queue Implementations Comparison**
| Feature | `LinkedList` | `ArrayDeque` | `PriorityQueue` |
|---------|-------------|-------------|----------------|
| **Ordering** | FIFO | FIFO | Priority-based |
| **Null Allowed** | Yes | No | No |
| **Thread-Safe** | No | No | No |
| **Best For** | General purpose | High-performance | Sorted processing |
| **Time Complexity** | O(1) ops | O(1) ops | O(log n) insert/remove |

---

### **Example Usage**
```java
Queue<String> queue = new LinkedList<>();

// Core operations
queue.add("First");       // Throws exception if full
queue.offer("Second");    // Returns false if full
String head = queue.poll(); // "First"

// PriorityQueue example
Queue<Integer> pq = new PriorityQueue<>();
pq.add(3); pq.add(1); pq.add(2);
while (!pq.isEmpty()) {
    System.out.println(pq.poll()); // Prints 1, 2, 3
}

// ArrayDeque as double-ended queue
Deque<String> deque = new ArrayDeque<>();
deque.offerFirst("Front");
deque.offerLast("End");
```

---

### **Key Characteristics**
1. **FIFO Principle**: First-In-First-Out (except `PriorityQueue`)
2. **Blocking Variants**: `BlockingQueue` implementations (`ArrayBlockingQueue`, `LinkedBlockingQueue`) for thread-safe operations
3. **Preferred Implementations**:
   - `ArrayDeque` for single-threaded performance
   - `ConcurrentLinkedQueue` for thread-safe needs
   - `PriorityQueue` for sorted processing

---

### **Important Notes**
- Always check `isEmpty()` before `remove()`/`element()` to avoid exceptions
- `PriorityQueue` orders elements either:
  - Natural ordering (via `Comparable`)
  - Custom `Comparator`
- For stack-like behavior (LIFO), use `Deque` interface methods:
  ```java
  Deque<String> stack = new ArrayDeque<>();
  stack.push("A"); // LIFO
  stack.pop();
  ```

For complete details, refer to the [Java 17 Queue Javadoc](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/util/Queue.html).
