Here's a comprehensive list of all Java data types with their ranges, explanations, and important edge cases:

### **1. Primitive Data Types**

| Data Type | Size (bits) | Range | Default Value | Notes & Edge Cases |
|-----------|------------|-------|--------------|--------------------|
| **byte**  | 8          | -128 to 127 | 0 | Overflow wraps around (e.g., `(byte)128 → -128`) |
| **short** | 16         | -32,768 to 32,767 | 0 | Rarely used; prefer `int` unless memory-critical |
| **int**   | 32         | -2³¹ to 2³¹-1 (-2.1B to 2.1B) | 0 | Common choice; `Integer.MAX_VALUE + 1 → Integer.MIN_VALUE` |
| **long**  | 64         | -2⁶³ to 2⁶³-1 | 0L | Suffix with `L` (e.g., `10000000000L`); time calculations often use this |
| **float** | 32         | ±1.4E-45 to ±3.4E+38 | 0.0f | 7 decimal digits precision; `0.1f + 0.2f != 0.3f` (floating-point precision) |
| **double**| 64         | ±4.9E-324 to ±1.8E+308 | 0.0d | 15 decimal digits precision; default for decimals |
| **char**  | 16         | '\u0000' to '\uffff' (0 to 65,535) | '\u0000' | Represents Unicode; `char c = 65;` → `'A'` |
| **boolean**| 1 (JVM-dependent) | `true` or `false` | false | Not convertible to/from integers |

---

### **2. Reference Data Types**

| Data Type | Storage | Default Value | Notes & Edge Cases |
|-----------|--------|--------------|--------------------|
| **String** | Heap (immutable) | `null` | `"a" + "b"` creates new object; `==` compares references |
| **Array**  | Heap   | `null`       | `int[] arr = new int[5];` (fixed size) |
| **Class Objects** | Heap | `null` | Requires `new` keyword; e.g., `Object obj = new Object();` |
| **Interfaces** | Heap | `null` | Cannot be instantiated directly |

---

### **3. Special Notes & Edge Cases**

#### **Integer Overflow**
```java
int max = Integer.MAX_VALUE;  // 2,147,483,647
max += 1;                     // Wraps to -2,147,483,648
```

#### **Floating-Point Precision**
```java
System.out.println(0.1 + 0.2);  // Prints 0.30000000000000004
```

#### **Char Arithmetic**
```java
char c = 'A';
c += 32;  // Now 'a' (ASCII manipulation)
```

#### **Boolean Size**
- JVM-dependent (typically 1 bit, but may use 1 byte)

#### **Null Pitfalls**
```java
String s = null;
System.out.println(s.length());  // NullPointerException
```

#### **Autoboxing**
```java
Integer a = 100, b = 100;
System.out.println(a == b);   // true (cached values)
Integer c = 200, d = 200;
System.out.println(c == d);   // false (new objects)
```

---

### **4. Memory Comparison**

| Type       | Example              | Memory (approx) |
|------------|----------------------|----------------|
| `byte`     | `byte b = 100;`      | 1 byte         |
| `int`      | `int i = 100_000;`   | 4 bytes        |
| `double`   | `double d = 3.14;`   | 8 bytes        |
| `String`   | `String s = "Hi";`   | ~40 bytes (object overhead) |

---

### **When to Use Which?**
- **For precise calculations**: `BigDecimal` (not primitive)
- **For memory efficiency**: `byte`/`short` in large arrays
- **General-purpose**: `int`, `double`, `boolean`
- **Text processing**: `char`/`String`

This covers all 8 primitives + reference types in Java with their key edge cases. Let me know if you'd like examples for any specific scenario!
