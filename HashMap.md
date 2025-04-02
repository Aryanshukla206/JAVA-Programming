HashMap
---

### **1. Core Operations**
| Method | Description | Time Complexity |
|--------|------------|-----------------|
| `put(K key, V value)` | Inserts key-value pair | O(1) avg |
| `get(Object key)` | Returns value for key | O(1) avg |
| `remove(Object key)` | Removes key-value pair | O(1) avg |
| `containsKey(Object key)` | Checks if key exists | O(1) avg |
| `containsValue(Object value)` | Checks if value exists | O(n) |
| `clear()` | Removes all mappings | O(n) |

---

### **2. Bulk Operations**
| Method | Description |
|--------|------------|
| `putAll(Map<? extends K, ? extends V> m)` | Copies all from another map |
| `replaceAll(BiFunction<? super K, ? super V, ? extends V> function)` | Replaces all values |

---

### **3. View Collections**
| Method | Returns |
|--------|---------|
| `keySet()` | Set of keys |
| `values()` | Collection of values |
| `entrySet()` | Set of key-value pairs |

---

### **4. Conditional Updates (Java 8+)**
| Method | Description |
|--------|------------|
| `putIfAbsent(K key, V value)` | Adds if key absent |
| `compute(K key, BiFunction<? super K, ? super V, ? extends V> remappingFunction)` | Computes new value |
| `computeIfAbsent(K key, Function<? super K, ? extends V> mappingFunction)` | Computes if key absent |
| `computeIfPresent(K key, BiFunction<? super K, ? super V, ? extends V> remappingFunction)` | Computes if key present |
| `merge(K key, V value, BiFunction<? super V, ? super V, ? extends V> remappingFunction)` | Merges values |

---

### **5. Size/Status Checks**
| Method | Description |
|--------|------------|
| `size()` | Returns number of mappings |
| `isEmpty()` | Checks if map is empty |

---

### **6. Default Methods (Java 8+)**
| Method | Description |
|--------|------------|
| `getOrDefault(Object key, V defaultValue)` | Returns value or default |
| `forEach(BiConsumer<? super K, ? super V> action)` | Iterates all entries |

---

### **7. Special Cases**
| Method | Description |
|--------|------------|
| `clone()` | Shallow copy (avoid in production) |
| `replace(K key, V oldValue, V newValue)` | Conditional replace |
| `replace(K key, V value)` | Unconditional replace |

---

### **Example Usage**
```java
HashMap<String, Integer> map = new HashMap<>();

// Core operations
map.put("Apple", 1);
map.put("Banana", 2);
int val = map.get("Apple"); // 1

// Bulk operations
map.putAll(Map.of("Cherry", 3, "Date", 4));

// Java 8+ features
map.computeIfAbsent("Elderberry", k -> k.length()); // Adds "Elderberry" → 10

// Iteration
map.forEach((k, v) -> System.out.println(k + ": " + v));
```

---

### **Key Characteristics**
1. **Null Keys/Values**: Allows one null key and multiple null values
2. **Unordered**: No guarantee of iteration order (use `LinkedHashMap` for insertion-order)
3. **Thread-Unsafe**: Use `ConcurrentHashMap` or `Collections.synchronizedMap()` for thread safety
4. **Load Factor**: Default 0.75 (resizes when 75% full)

---

### **Important Notes**
- All `O(1)` operations assume good hash distribution (worst case is `O(n)` for collisions)
- Methods like `entrySet()` return views that reflect map changes
- Java 17 added no new HashMap methods (last major update was Java 8)

For complete details, refer to the [Java 17 HashMap Javadoc](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/util/HashMap.html).
