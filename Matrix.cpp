#include <algorithm>
#include <limits> // std::numeric_limits

#include "Matrix.hpp"
//your code here...



//You'll implement this method in part two of the project
Matrix Matrix::invMod(Scalar m) const
{
  //For part two of the project you will use the elementary row operations to peform
  //Gauss-Jordan Elimination to find the inverse of a matrix...you may start early, if you wish
  return Matrix();
}

//Instructor provided methods and functions below; do not modify
Matrix::Matrix(const Matrix &A): m(A.m) 
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
  std::size_t n = A.elements.size() / A.m;
  
  for(std::size_t i = 0; i < A.elements.size(); i++)
    {
      //beginning of a row (column zero): insert a line break unless we're in the first row
      if(i != 0 && i % n == 0)
	os << std::endl;
      
      os << A.elements[i] << " ";
    }

  return os;
}
