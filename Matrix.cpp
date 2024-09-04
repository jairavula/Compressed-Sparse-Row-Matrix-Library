#include <algorithm>
#include <limits> // std::numeric_limits

#include "Matrix.hpp"
//your code here...

// Default constructor
Matrix::Matrix() : elements(), col_ind(), row_ind(3), n(2) {}

Matrix::Matrix(std::size_t m, std::size_t n): elements(), col_ind(), row_ind(m + 1, 0), n(n){};


  /**
   * Parameterized constructor.  Use the parameters to set the matrix elements; if parameters are inconsistent then create a 0-by-0 (empty) matrix.
   * @param A - values for matrix elements, specified using a dense, row-major order vector
   * @param n - number of columns for the new matrix.
   */ 
  Matrix::Matrix(const std::vector<Elem> &A, std::size_t n) : n(n), row_ind(A.size()/n,0) {

      //populate the elements array
      for(int i = 0; i < A.size(); i++){
        if (A[i] != 0){
          elements.push_back(A[i]);
        }
      }

      // populate the column indices array
      for(int i = 0; i < A.size(); i++) {
        if (A[i] != 0) {
          col_ind.push_back(i % n);
        }

        //populate the row pointer array
        if (i % n == 0 && i != 0) {
          row_ind.push_back(i);
        }
      }
      // last element in row pointer array holds the number of elements + 1
      row_ind.push_back(elements.size() + 1);

  };

  //  /**
  //  * Another parameterized constructor.  Use the parameters to set the matrix element; if parameters are inconsistent then create a 0-by-0 (empty) matrix.
  //  * @param ptr_A - pointer to a dense, row-major order array of matrix elements
  //  * @param m - number of rows for the new matrix.
  //  * @param n - number of columns for the new matrix.
  //  */ 
  // Matrix::Matrix(const Elem *ptr_A, std::size_t m, std::size_t n) : n(n){

  // }


  Elem Matrix::e(std::size_t i, std::size_t j) const{
    // 1. Ensure i is a valid row number
    // 2. Iterate through the column index using the possible entry indices from the row index
    // 3. Look for a column index matching j


    if (i < row_ind.size()-1){
      for (std::size_t index = row_ind[i]; index < row_ind[i+1]; index++) {
      if (col_ind[index] == j){
        return elements[index];
        }
      }
    } 
    // return 0 if value not found
    return 0;
  }
  

  bool Matrix::e(std::size_t i, std::size_t j, Elem aij){

    // TODO: Double check logic here
    if (j > n || i > row_ind.size()-1) {
      return false;
    }
  
    if (i < row_ind.size()-1){
      for (std::size_t index = row_ind[i]; index < row_ind[i+1]; index++) {
      if (col_ind[index] == j){
        elements[index] = aij;
        return true;
        }
      }
    }
    return false;
  
  }

  std::size_t Matrix::size(std::size_t dim) const {
    if (dim == 1){
      return (row_ind.size()-1);
    } else if (dim == 2) {
      return n;
    } else return 0;
  }

  bool Matrix::equal( const Matrix& rhs ) const {
    // Ensure full CSR representations vectors are equal in both matrices
    if ((rhs.elements == elements) && (rhs.col_ind == col_ind) && (rhs.row_ind == row_ind)){
      return true;
    }
    else return false;
  }

  Matrix Matrix::add( const Matrix &rhs ) const {
    return Matrix();
  }
  Matrix Matrix::sub( const Matrix &rhs ) const{
    return Matrix();
  }
  Matrix Matrix::mult( const Matrix &rhs ) const{
    return Matrix();
  }
  Matrix Matrix::mult( Scalar k ) const{
    return Matrix();
  }
  Matrix Matrix::pow( Scalar k ) const{
    return Matrix();
  }
  Matrix Matrix::trans() const{
    return Matrix();
  }
  Matrix Matrix::cat(const Matrix &rhs, std::size_t dim) const{
    return Matrix();
  }

  bool Matrix::rowSwitch(std::size_t i, std::size_t j){
    return true;
  }

  bool Matrix::rowMult(std::size_t i, Scalar k){
    return true;
  }
  

  bool Matrix::rowAdd(std::size_t i, std::size_t j, Scalar k){
    return true;
  }
  











//You'll implement this method in part two of the project
Matrix Matrix::invMod(Scalar m) const
{
  //For part two of the project you will use the elementary row operations to peform
  //Gauss-Jordan Elimination to find the inverse of a matrix...you may start early, if you wish
  return Matrix();
}

//Instructor provided methods and functions below; do not modify
Matrix::Matrix(const Matrix &A): n(A.n) 
{
  elements = A.elements;
  col_ind = A.col_ind;
  row_ind = A.row_ind;
}

Matrix& Matrix::operator=(Matrix A)
{
  std::swap(elements, A.elements);
  std::swap(col_ind, A.col_ind);
  std::swap(row_ind, A.row_ind);

  return *this;
}
  
std::ostream& operator<<(std::ostream& os, const Matrix &A)
{
  //number of columns
  std::size_t n = A.elements.size() / A.n;
  
  for(std::size_t i = 0; i < A.elements.size(); i++)
    {
      //beginning of a row (column zero): insert a line break unless we're in the first row
      if(i != 0 && i % A.n == 0)
	os << std::endl;
      
      os << A.elements[i] << " ";
    }

  return os;
}
