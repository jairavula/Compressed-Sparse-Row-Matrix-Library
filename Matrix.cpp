#include <algorithm>
#include <limits> // std::numeric_limits

#include "Matrix.hpp"
//your code here...

// Default constructor
Matrix::Matrix() : elements(), col_ind(), row_ind(3), n(2) {}

Matrix::Matrix(std::size_t m, std::size_t n): elements(), col_ind(), row_ind(m + 1, 0), n(n){};

  Matrix::Matrix(const std::vector<Elem> &A, std::size_t n) : n(n), row_ind() {

    std::size_t nonzero_elements = 0;
    row_ind.push_back(0);

      //populate the elements array and corresponding column index. 
      for(int i = 0; i < A.size(); i++){
        if (A[i] != 0){
          elements.push_back(A[i]);
          col_ind.push_back(i % n);
          nonzero_elements++;
        }
        //push back the index of non-zero elements to the row index array when a row finishes
        if ((i + 1) % n == 0) {
            row_ind.push_back(nonzero_elements);
        }
      }
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
    // Ensure full CSR representation vectors are equal in both matrices
    if ((rhs.elements == elements) && (rhs.col_ind == col_ind) && (rhs.row_ind == row_ind)){
      return true;
    }
    else return false;
  }

  Matrix Matrix::add( const Matrix &rhs ) const {
    if( row_ind.size() != rhs.row_ind.size() || n != rhs.n){
      return Matrix(0,0);
    }
    // create empty m x n array and start row_ind array with a 0
    Matrix summed_matrix(row_ind.size()-1, n); 
    summed_matrix.row_ind.push_back(0);

    // iteration through row index array
    for(int i=0; i < row_ind.size()-1; i++){
      // scope on one row at a time
      int matrix1_row_start = row_ind[i];
      int matrix1_row_end = row_ind[i+1];

      int matrix2_row_start = rhs.row_ind[i];
      int matrix2_row_end = rhs.row_ind[i+1];
      
      while (matrix1_row_start != matrix1_row_end || matrix2_row_start != matrix2_row_end) {
        if (matrix1_row_start == matrix1_row_end) {
        summed_matrix.elements.push_back(rhs.elements[matrix2_row_start]);
        summed_matrix.col_ind.push_back(rhs.col_ind[matrix2_row_start]);
        matrix2_row_start++;
        } else if (matrix2_row_start == matrix2_row_end) {
        summed_matrix.elements.push_back(elements[matrix1_row_start]);
        summed_matrix.col_ind.push_back(col_ind[matrix1_row_start]);
        matrix1_row_start++;
        } else if(col_ind[matrix1_row_start] == rhs.col_ind[matrix2_row_start]){
          summed_matrix.elements.push_back(elements[matrix1_row_start] + rhs.elements[matrix2_row_start]);
          summed_matrix.col_ind.push_back(col_ind[matrix1_row_start]);
          matrix1_row_start++;
          matrix2_row_start++;
        } else if(col_ind[matrix1_row_start] < rhs.col_ind[matrix2_row_start]) {
          summed_matrix.elements.push_back(elements[matrix1_row_start]);
          summed_matrix.col_ind.push_back(col_ind[matrix1_row_start]);
          matrix1_row_start++;
        } else if(col_ind[matrix1_row_start] > rhs.col_ind[matrix2_row_start]){
          summed_matrix.elements.push_back(elements[matrix2_row_start]);
          summed_matrix.col_ind.push_back(col_ind[matrix2_row_start]);
          matrix2_row_start++;
        }
      }
      summed_matrix.row_ind.push_back(summed_matrix.elements.size()-1);
    }  
    return summed_matrix;
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
    // Loop through each row
    for (std::size_t i = 0; i < A.row_ind.size() - 1; i++) {
        std::size_t row_start = A.row_ind[i];      // Start of the current row's non-zero elements
        std::size_t row_end = A.row_ind[i + 1];    // End of the current row's non-zero elements
        std::size_t column_pointer = 0;            // Column index tracker
        
      // Execute while in scope of row
        while (column_pointer < A.n) {
            
            if (row_start < row_end && A.col_ind[row_start] == column_pointer) {
                //if there is a non-zero entry in the specified column, print the entry
                os << A.elements[row_start] << " ";
                row_start++;  
            } else {
                //if not print a 0
                os << 0 << " ";
            }
            // go to next column
            column_pointer++; 
        }
        // new line after row is completed
        os << std::endl; 
    }
    return os;
}

