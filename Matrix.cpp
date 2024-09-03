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
      elements = A;
      // populate the column indices array
      int row_num = 0;
      for(int i = 0; i < elements.size(); i++) {
        if (elements[i] != 0) {
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
