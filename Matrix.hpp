#ifndef _MATRIX_HPP_
#define _MATRIX_HPP_

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

/**
 * For mathematical classes in general we pay particular attention to the underlying data type, as we rely on knowing the type to optimize and perform error correction/detection
 */
typedef std::int_fast32_t Elem; //a matrix element, an integer at least 32 bits in size
typedef std::int_fast32_t Scalar; //a scalar value by which an element may be multiplied

/**
 * A basic C++ class for two-dimensional matrix arithmetic that uses the compressed sparse row (CSR) format for storing the matrix.
 *
 * NOTES: 
 *   1. Your matrix must be stored in the CSR format and you must use the private data members provided without modification
 *   2. You may *not* add additional data members, either public or private, as doing so will result in incomplete copies/assignments
 *   3. You may *not* add/modify public methods but you may define your own private methods
 */

class Matrix
{
public:
  /**
   * Default constructor. It should create an empty 2-by-2 matrix.
   */ 
  Matrix();

  /**
   * Parameterized constructor.  It should create an empty matrix of user-supplied dimensions.
   * @param m - number of row for the new matrix.
   * @param n - number of columns for the new matrix.
   */ 
  Matrix(std::size_t m, std::size_t n);
  
  /**
   * Parameterized constructor.  Use the parameters to set the matrix elements; if parameters are inconsistent then create a 0-by-0 (empty) matrix.
   * @param A - values for matrix elements, specified using a dense, row-major order vector
   * @param n - number of columns for the new matrix.
   */ 
  Matrix(const std::vector<Elem> &A, std::size_t n);

  /**
   * Another parameterized constructor.  Use the parameters to set the matrix element; if parameters are inconsistent then create a 0-by-0 (empty) matrix.
   * @param ptr_A - pointer to a dense, row-major order array of matrix elements
   * @param m - number of rows for the new matrix.
   * @param n - number of columns for the new matrix.
   */ 
  Matrix(const Elem *ptr_A, std::size_t m, std::size_t n);

  /**
   * Copy constructor.  Necessary if we have const data members; do not modify.
   */
  Matrix(const Matrix &A);

  /**
   * Copy assignment operator.  Necessary if we have const data members; do not modify.
   */
  Matrix& operator=(Matrix A);
  
  /**
   * Returns the element at specified row, column index.
   * @param i - row index of object.
   * @param j - column index of object.
   * @return element at specified row, column index or smallest possible value for Elem if index is invalid.
   */ 
  Elem e(std::size_t i, std::size_t j) const;
  
  /**
   * Sets the element at specified row, column index to given value; if either index is invalid matrix should not be modified.
   * @param i - row index of object to set.
   * @param j - column index of object to set.
   * @param aij - value for element at index i, j
   * @return true if set is successful, false otherwise.
   */
  bool e(std::size_t i, std::size_t j, Elem aij);

  /**
   * Returns the size of the matrix along a given dimension (i.e., number of row(s) or column(s))
   * @param dim - 1 for row, 2 for column
   * @return the length of the dimension specified, if dimension is not valid return 0
   */ 
  std::size_t size(std::size_t dim) const;
  
 /**
   * Returns true if the matrices this and rhs are the same, false otherwise.
   * @param rhs - the Matrix object to compare to this object.
   * @return true if all the elements in both objects are the same, false otherwise.
   */ 
  bool equal( const Matrix& rhs ) const;

  /**
   * Creates and returns a new Matrix object representing the matrix addition of two Matrix objects.
   * @return a new Matrix object that contains the appropriate summed elements, a 0-by-0 matrix if matrices can't be added.
   * @param rhs - the Matrix object to add to this object.
   */
  Matrix add( const Matrix &rhs ) const;

  /**
   * Creates and returns a new Matrix object representing the matrix subtraction of two Matrix objects.
   * @return a new Matrix object that contains the appropriate difference elements, a 0-by-0 matrix if matrices can't be subtracted.
   * @param rhs - the Matrix object to subtract from this object.
   */
  Matrix sub( const Matrix &rhs ) const;

  /**
   * Creates and returns a new Matrix object that is the multiplication of this and the given Matrix object.
   * @return a new Matrix object that contains the multiplication of this and the given Matrix object, a 0-by-0 matrix if matrices can't be multiplied.
   * @param rhs - the Matrix object to multiply with this object.
   */
  Matrix mult( const Matrix &rhs ) const;

  /**
   * Creates and returns a new Matrix object that is the multiplication of this and the given scalar.
   * @return a new Matrix object that contains the multiplication of this and the given scalar.
   * @param rhs - the scalar value to multiply with this object.
   */
  Matrix mult( Scalar k ) const;

  /**
   * Creates and returns a new Matrix object that is the power of this.
   * @return a new Matrix object that raises this and to the given power.
   * @param k - the power to which this object should be raised, a 0-by-0 matrix if matrix can't be raised to power.
   */
  Matrix pow( Scalar k ) const;

  /**
   * Creates and returns a new Matrix object that is the transpose of this.
   * @return a new Matrix object that is the transpose of this object.
   */
  Matrix trans() const;

  /**
   * Creates and returns a new Matrix object that is the concatenation of this and the given Matrix object.
   * @return a new Matrix object that is the vertical or horizontal concatenation of two matrices, a 0-by-0 matrix if matrices can't be concatenated.
   * @param dim - 1 for vertical cat (RHS below this), 2 for horizontal cat (RHS to the right of this)
   */
  Matrix cat(const Matrix &rhs, std::size_t dim) const;

  /**
   * Creates and returns a new Matrix object that is the modular multiplicative inverse, using the modulus m, of this
   * @return a new Matrix object A such that this*A = 1 (mod m), a 0-by-0 matrix if this is non-invertible.
   * @param m - modulus for the modular multiplicative inverse
   */
  Matrix invMod(Scalar m) const;
  
  /**
   * Switch (swap) rows within the matrix (in-place operation)
   * Ri <=> Rj
   * @return true if rows i and j were switched
   * @param i - row number
   * @param j - row number
   */
  bool rowSwitch(std::size_t i, std::size_t j);

  /**
   * Multiply a row within the matrix by a scalar (in-place operation)
   * Ri = k*Ri
   * @return true if rows i was multiplied by scalar k
   * @param i - row number
   * @param k - scalar to multiply row i by
   */
  bool rowMult(std::size_t i, Scalar k);

  /**
   * Replace a row with the sum of itself and the scalar multiple of another row (in-place operation)
   * Ri = Ri + k*Rj
   * @return true if row i Ri
   * @param i - row number
   * @param k - scalar to multiply row i by
   */
  bool rowAdd(std::size_t i, std::size_t j, Scalar k);
  
  /**
   * Outputs this Matrix object on the given ostream (for debugging).
   * @param os - the ostream object to use to output.
   * @param A - the Matrix we want to `display' or `print'
   * @return - Matrix elements in two-dimensional form
   */ 
  friend std::ostream& operator<<(std::ostream& os, const Matrix &A);

private:
  std::vector<Elem> elements; //The elements of our matrix; size is equal to the number of non-zero elements
  const std::size_t m; //number of rows in the matrix (const as matrices shouldn't change size after being defined)
  std::vector<std::size_t> col_ind; //column for the i-th element of elements (i.e., col_ind[i] is the column for elements[i])
  std::vector<std::size_t> row_ind; //number of non-zero elements above row i (i.e., row_ind[i] = sum(non-zero elements in
                                    //rows zero to i-1)), remember this should have an extra entry at the end of the vector where
                                    //we store the total number of non-zero entries in elements
  //Feel free to add your own private methods below
};
#endif
