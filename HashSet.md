HashSet
---

### **1. Core Operations**
| Method | Description | Time Complexity |
|--------|------------|-----------------|
| `add(E e)` | Adds element if not present | O(1) avg |
| `remove(Object o)` | Removes element | O(1) avg |
| `contains(Object o)` | Checks if element exists | O(1) avg |
| `clear()` | Removes all elements | O(n) |

---

### **2. Bulk Operations**
| Method | Description |
|--------|------------|
| `addAll(Collection<? extends E> c)` | Adds all elements from collection |
| `removeAll(Collection<?> c)` | Removes all elements in collection |
| `retainAll(Collection<?> c)` | Keeps only elements in collection |
| `removeIf(Predicate<? super E> filter)` | Removes elements matching condition |

---

### **3. Set Characteristics**
| Method | Description |
|--------|------------|
| `size()` | Returns number of elements |
| `isEmpty()` | Checks if set is empty |
| `iterator()` | Returns iterator over elements |

---

### **4. Conversion Methods**
| Method | Returns |
|--------|---------|
| `toArray()` | Object array of elements |
| `toArray(T[] a)` | Typed array of elements |
| `spliterator()` | Parallel iteration support (Java 8+) |

---

### **5. Java 8+ Features**
| Method | Description |
|--------|------------|
| `forEach(Consumer<? super E> action)` | Performs action on each element |
| `stream()` | Returns sequential Stream |
| `parallelStream()` | Returns parallel Stream |

---

### **6. Special Cases**
| Method | Description |
|--------|------------|
| `clone()` | Shallow copy (avoid in production) |

---

### **Example Usage**
```java
HashSet<String> fruits = new HashSet<>();

// Core operations
fruits.add("Apple");
fruits.add("Banana");
boolean hasApple = fruits.contains("Apple"); // true

// Bulk operations
fruits.addAll(List.of("Cherry", "Date"));

// Java 8+ features
fruits.removeIf(f -> f.startsWith("A")); // Removes "Apple"

// Iteration
fruits.forEach(System.out::println);
```

---

### **Key Characteristics**
1. **Null Elements**: Allows exactly one `null` element
2. **Unordered**: No guarantee of iteration order (use `LinkedHashSet` for insertion-order)
3. **Thread-Unsafe**: Use `Collections.synchronizedSet()` or `ConcurrentHashMap.newKeySet()` for thread safety
4. **Backed by HashMap**: Implemented using a `HashMap<E, Object>` internally

---

### **Important Notes**
- All `O(1)` operations assume good hash distribution (worst case is `O(n)` for collisions)
- Methods like `iterator()` return fail-fast iterators (throw `ConcurrentModificationException` if modified during iteration)
- Java 17 added no new HashSet methods (last major update was Java 8)

For complete details, refer to the [Java 17 HashSet Javadoc](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/util/HashSet.html).
