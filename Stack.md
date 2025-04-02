Stack
---

### **1. Core Stack Operations (LIFO)**
| Method | Description | Time Complexity |
|--------|------------|-----------------|
| `push(E item)` | Pushes element onto stack | O(1) |
| `pop()` | Removes and returns top element | O(1) |
| `peek()` | Returns top element without removal | O(1) |
| `empty()` | Checks if stack is empty | O(1) |

---

### **2. Vector Inheritance (Legacy Methods)**
| Method | Description |
|--------|------------|
| `search(Object o)` | Returns 1-based position from top | O(n) |
| `addElement(E item)` | Legacy alias for `push()` |
| `elementAt(int index)` | Gets element at position (0=bottom) |
| `removeElementAt(int index)` | Removes element at position |

---

### **3. List Operations (Inherited from Vector)**
| Method | Description |
|--------|------------|
| `add(E e)` | Appends element (to top) |
| `remove(int index)` | Removes element at position |
| `get(int index)` | Gets element at position |
| `size()` | Returns element count |

---

### **4. Iteration Methods**
| Method | Description |
|--------|------------|
| `iterator()` | Returns iterator (top-to-bottom) |
| `listIterator()` | Returns bidirectional iterator |
| `forEach(Consumer<? super E> action)` | Java 8 iteration |

---

### **5. Utility Methods**
| Method | Description |
|--------|------------|
| `clear()` | Removes all elements |
| `contains(Object o)` | Checks if element exists |
| `copyInto(Object[] anArray)` | Copies elements to array |

---

### **Example Usage**
```java
Stack<String> stack = new Stack<>();

// Core stack operations
stack.push("First");
stack.push("Second");
String top = stack.peek(); // "Second"
String popped = stack.pop(); // "Second"

// Legacy methods
int pos = stack.search("First"); // 1 (top=1, bottom=N)

// Modern iteration
stack.forEach(System.out::println);
```

---

### **Key Characteristics**
1. **Thread-Safe**: Inherits synchronization from `Vector` (but consider `ArrayDeque` for better performance)
2. **LIFO Order**: Last-In-First-Out behavior
3. **Legacy Class**: Consider `Deque` implementations (e.g., `ArrayDeque`) for new code
4. **Null Elements**: Allows multiple `null` values

---

### **Important Notes**
- The `Stack` class extends `Vector`, which gives it many legacy collection methods
- For new code, Oracle recommends using `ArrayDeque` instead (more performant)
- Java 17 maintains backward compatibility but hasn't added new methods to `Stack`

For complete details, refer to the [Java 17 Stack Javadoc](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/util/Stack.html).
