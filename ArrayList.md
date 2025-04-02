ArrayList
---

### **1. Core List Operations**
| Method | Description | Time Complexity |
|--------|------------|-----------------|
| `add(E e)` | Appends element | O(1) amortized |
| `add(int index, E e)` | Inserts at position | O(n) |
| `get(int index)` | Returns element | O(1) |
| `set(int index, E e)` | Replaces element | O(1) |
| `remove(int index)` | Removes by position | O(n) |
| `remove(Object o)` | Removes first occurrence | O(n) |
| `clear()` | Removes all elements | O(n) |

---

### **2. Bulk Operations**
| Method | Description |
|--------|------------|
| `addAll(Collection<? extends E> c)` | Appends collection |
| `addAll(int index, Collection<? extends E> c)` | Inserts collection at position |
| `removeAll(Collection<?> c)` | Removes all elements in collection |
| `retainAll(Collection<?> c)` | Keeps only elements in collection |
| `removeIf(Predicate<? super E> filter)` | Removes matching elements (Java 8+) |

---

### **3. Size/Capacity Management**
| Method | Description |
|--------|------------|
| `size()` | Returns element count |
| `isEmpty()` | Checks if empty |
| `ensureCapacity(int minCapacity)` | Pre-allocates space |
| `trimToSize()` | Trims capacity to current size |

---

### **4. Search Operations**
| Method | Description |
|--------|------------|
| `contains(Object o)` | Checks for element |
| `indexOf(Object o)` | Returns first index |
| `lastIndexOf(Object o)` | Returns last index |

---

### **5. Iteration/Conversion**
| Method | Description |
|--------|------------|
| `iterator()` | Returns basic iterator |
| `listIterator()` | Returns bidirectional iterator |
| `listIterator(int index)` | Starts iterator at position |
| `toArray()` | Converts to Object[] |
| `toArray(T[] a)` | Converts to typed array |
| `spliterator()` | Parallel iteration (Java 8+) |

---

### **6. Java 8+ Features**
| Method | Description |
|--------|------------|
| `forEach(Consumer<? super E> action)` | Performs action on elements |
| `replaceAll(UnaryOperator<E> operator)` | Replaces all elements |
| `sort(Comparator<? super E> c)` | Sorts with comparator |
| `stream()` | Returns sequential stream |
| `parallelStream()` | Returns parallel stream |

---

### **7. Sub-List Operations**
| Method | Description |
|--------|------------|
| `subList(int fromIndex, int toIndex)` | Returns view of portion |

---

### **Example Usage**
```java
ArrayList<String> colors = new ArrayList<>();

// Core operations
colors.add("Red");
colors.add(0, "Blue"); // Insert at front
colors.set(1, "Green"); // Replace

// Bulk operations
colors.addAll(List.of("Yellow", "Purple"));

// Java 8+ features
colors.removeIf(c -> c.startsWith("B"));
colors.forEach(System.out::println);

// Conversion
String[] colorArray = colors.toArray(new String[0]);
```

---

### **Key Characteristics**
1. **Dynamic Resizing**: Grows by ~50% when full (default capacity: 10)
2. **Random Access**: Implements `RandomAccess` marker interface
3. **Null Elements**: Allows multiple `null` values
4. **Not Thread-Safe**: Use `Collections.synchronizedList()` for thread safety

---

### **Performance Notes**
- **Indexed Access**: O(1) time for `get/set`
- **Insertion/Removal**: O(n) time for middle operations
- **Iteration**: Faster than `LinkedList` for sequential access

---

### **Alternatives**
- **`LinkedList`**: Better for frequent insertions/deletions
- **`CopyOnWriteArrayList`**: Thread-safe for read-heavy workloads
- **`Vector`**: Thread-safe legacy version (slower due to synchronization)

For complete details, refer to the [Java 17 ArrayList Javadoc](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/util/ArrayList.html).
