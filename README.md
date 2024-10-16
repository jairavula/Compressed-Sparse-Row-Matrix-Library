# Compressed Sparse Row Matrix Library (CSR) - C++

This library implements a **Compressed Sparse Row (CSR) matrix** in C++, optimizing storage and operations on sparse matrices. It supports addition, subtraction, multiplication, elementary row operations, and row reduction.

## Features

- **Efficient Storage**: Stores only non-zero elements.
- **Matrix Operations**: Add, subtract, and multiply matrices.
- **Elementary Row Operations**: Swap, scale, and add rows.
- **Row Reduction**: Perform Gaussian elimination for solving linear systems.

## How It Works

The CSR format uses:
1. **`values[]`**: Stores non-zero matrix values.
2. **`col_ind[]`**: Column indices corresponding to `values[]`.
3. **`row_ptr[]`**: Pointers to the start of each row.

### Matrix Operations

- **Addition/Subtraction**: Efficient operations across matching non-zero elements.
- **Multiplication**: Optimized for sparse data.
- **Row Operations**: Easily manipulate rows for row reduction.

## Benefits

- **Memory Efficient**: Ideal for sparse matrices.
- **Fast Access**: Efficient row-based modifications.
- **Performance**: Minimizes operations, improving performance.

## Example Usage

```cpp
#include "CSRMatrix.h"

CSRMatrix matrixA(...);
CSRMatrix result = matrixA.add(matrixB);
matrixA.rowReduce();
