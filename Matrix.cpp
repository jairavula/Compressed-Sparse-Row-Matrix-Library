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

Matrix Matrix::add(const Matrix &rhs) const {
    if (row_ind.size() != rhs.row_ind.size() || n != rhs.n) {
        return Matrix(0, 0);  // Inconsistent dimensions
    }

    // Create the result matrix and initialize with empty row_ind
    Matrix summed_matrix(0, n);

    // Iterate through each row of the matrix
    for (std::size_t i = 0; i < row_ind.size() - 1; ++i) {
        std::size_t matrix1_row_start = row_ind[i];
        std::size_t matrix1_row_end = row_ind[i + 1];

        std::size_t matrix2_row_start = rhs.row_ind[i];
        std::size_t matrix2_row_end = rhs.row_ind[i + 1];

        // Add non-zero elements from both matrices
        while (matrix1_row_start != matrix1_row_end || matrix2_row_start != matrix2_row_end) {
            if (matrix1_row_start == matrix1_row_end) {
                // If matrix 1 has no more elements, take from matrix 2
                summed_matrix.elements.push_back(rhs.elements[matrix2_row_start]);
                summed_matrix.col_ind.push_back(rhs.col_ind[matrix2_row_start]);
                matrix2_row_start++;
            } else if (matrix2_row_start == matrix2_row_end) {
                // If matrix 2 has no more elements, take from matrix 1
                summed_matrix.elements.push_back(elements[matrix1_row_start]);
                summed_matrix.col_ind.push_back(col_ind[matrix1_row_start]);
                matrix1_row_start++;
            } else if (col_ind[matrix1_row_start] == rhs.col_ind[matrix2_row_start]) {
                // If both matrices have elements at the same column, add them
                summed_matrix.elements.push_back(elements[matrix1_row_start] + rhs.elements[matrix2_row_start]);
                summed_matrix.col_ind.push_back(col_ind[matrix1_row_start]);
                matrix1_row_start++;
                matrix2_row_start++;
            } else if (col_ind[matrix1_row_start] < rhs.col_ind[matrix2_row_start]) {
                // Take element from matrix 1 if its column is smaller
                summed_matrix.elements.push_back(elements[matrix1_row_start]);
                summed_matrix.col_ind.push_back(col_ind[matrix1_row_start]);
                matrix1_row_start++;
            } else {
                // Take element from matrix 2 if its column is smaller
                summed_matrix.elements.push_back(rhs.elements[matrix2_row_start]);
                summed_matrix.col_ind.push_back(rhs.col_ind[matrix2_row_start]);
                matrix2_row_start++;
            }
        }

        // After finishing this row, update `row_ind` with the current number of elements
        summed_matrix.row_ind.push_back(summed_matrix.elements.size());
    }
    return summed_matrix;
}


Matrix Matrix::sub(const Matrix &rhs) const {
    if (row_ind.size() != rhs.row_ind.size() || n != rhs.n) {
        return Matrix(0, 0);  // Inconsistent dimensions
    }

    // Create the result matrix and initialize with empty row_ind
    Matrix summed_matrix(0, n);

    // Iterate through each row of the matrix
    for (std::size_t i = 0; i < row_ind.size() - 1; ++i) {
        std::size_t matrix1_row_start = row_ind[i];
        std::size_t matrix1_row_end = row_ind[i + 1];

        std::size_t matrix2_row_start = rhs.row_ind[i];
        std::size_t matrix2_row_end = rhs.row_ind[i + 1];

        // Track the number of elements at the start of the current row
        std::size_t summed_matrix_start = summed_matrix.elements.size();

        // Subtract non-zero elements from both matrices
        while (matrix1_row_start != matrix1_row_end || matrix2_row_start != matrix2_row_end) {
            if (matrix1_row_start == matrix1_row_end) {
                // If matrix 1 has no more elements, subtract from 0
                if (rhs.elements[matrix2_row_start] != 0) {
                    summed_matrix.elements.push_back(-rhs.elements[matrix2_row_start]);
                    summed_matrix.col_ind.push_back(rhs.col_ind[matrix2_row_start]);
                }
                matrix2_row_start++;
            } else if (matrix2_row_start == matrix2_row_end) {
                // If matrix 2 has no more elements, take from matrix 1
                if (elements[matrix1_row_start] != 0) {
                    summed_matrix.elements.push_back(elements[matrix1_row_start]);
                    summed_matrix.col_ind.push_back(col_ind[matrix1_row_start]);
                }
                matrix1_row_start++;
            } else if (col_ind[matrix1_row_start] == rhs.col_ind[matrix2_row_start]) {
                // If both matrices have elements at the same column, subtract them
                Elem result = elements[matrix1_row_start] - rhs.elements[matrix2_row_start];
                if (result != 0) {
                    summed_matrix.elements.push_back(result);
                    summed_matrix.col_ind.push_back(col_ind[matrix1_row_start]);
                }
                matrix1_row_start++;
                matrix2_row_start++;
            } else if (col_ind[matrix1_row_start] < rhs.col_ind[matrix2_row_start]) {
                // Take element from matrix 1 if its column is smaller
                if (elements[matrix1_row_start] != 0) {
                    summed_matrix.elements.push_back(elements[matrix1_row_start]);
                    summed_matrix.col_ind.push_back(col_ind[matrix1_row_start]);
                }
                matrix1_row_start++;
            } else {
                // Take element from matrix 2 if its column is smaller, subtract from 0
                if (rhs.elements[matrix2_row_start] != 0) {
                    summed_matrix.elements.push_back(-rhs.elements[matrix2_row_start]);
                    summed_matrix.col_ind.push_back(rhs.col_ind[matrix2_row_start]);
                }
                matrix2_row_start++;
            }
        }

        // After finishing this row, update `row_ind` if non-zero elements were added
        if (summed_matrix.elements.size() > summed_matrix_start) {
            summed_matrix.row_ind.push_back(summed_matrix.elements.size());
        } else {
            // Ensure the row_ind array is consistent even for empty rows
            summed_matrix.row_ind.push_back(summed_matrix_start);
        }
    }

    return summed_matrix;
}

  Matrix Matrix::mult( const Matrix &rhs ) const{
    // TODO: add logic for invalid multiplication

    //resultant matrix in Row-Major order
    std::vector<Elem> mult_result((row_ind.size()-1) * rhs.n, 0);
    std::size_t result_pointer = 0;
    for (std::size_t i = 0; i < row_ind.size()-1; i++){

      int row_start = row_ind[i];
      int row_end = row_ind[i+1];

      // buffers for multiplication operation
      std::vector<int> matrix1_row(n, 0);

      int column_index = 0;

      //update corresponding values in row buffer based on CSR
      while(row_start != row_end){
        if (col_ind[row_start] == column_index){
           matrix1_row[column_index] = elements[row_start];
           row_start++;
        }
        column_index++;
      }


      for(std::size_t col = 0; col < rhs.n; col++){

        std::vector<int> matrix2_column(rhs.row_ind.size()-1,0);
        int matrix2_column_pos = 0;

        // build vector of columns
        for(std::size_t j = 0; j < rhs.row_ind.size()-1; j++){
          int rhs_row_start = rhs.row_ind[j];
          int rhs_row_end = rhs.row_ind[j+1];
          while(rhs_row_start != rhs_row_end){
            //if the non-zero element is present at the column index, add it to the column vector  
            if(rhs.col_ind[rhs_row_start] == col){
            matrix2_column[matrix2_column_pos] = rhs.elements[rhs_row_start];
            }
            rhs_row_start++;
          }
          // increment matrix2_column index to next item
          matrix2_column_pos++;
        }
        std::size_t total = 0;
        for(std::size_t item = 0; item < n; item++){ 
          total += matrix1_row[item] * matrix2_column[item];
        }
        mult_result[result_pointer] = total;
        result_pointer++;
      }
    }
    
    return Matrix(mult_result, rhs.n);
  }

Matrix Matrix::mult(Scalar k) const {
    std::vector<Elem> new_elements(elements.size());
    //copy the existing matrix and scale non-zero elements.
    for (std::size_t entry = 0; entry < elements.size(); entry++) {
        new_elements[entry] = elements[entry] * k;
    }

    Matrix A(new_elements, n); 

    A.col_ind = col_ind;
    A.row_ind = row_ind;

    return A;
}


  Matrix Matrix::pow( Scalar k ) const{

    Matrix A(*this);

    A.col_ind = col_ind;
    A.row_ind = row_ind;

    for (int i = 0; i < k-1; i++){
      A = A.mult(*this);
    }

    return A;
  }

  Matrix Matrix::trans() const{

    //return empty matrix with same dim as original if matrix is all zero
    if(elements.size() == 0){
      return *this;
    }

    // element buffer used for matrix reconstruction
    std::vector<Elem> transpose_elems(elements.size()); 
    int transpose_elems_index = 0;

    //grab each column out of the matrix
    for(std::size_t col = 0; col < n; col++){

        //column buffer
        std::vector<int> matrix2_column(row_ind.size()-1,0);
        int matrix2_column_pos = 0;


        // build vector of columns
        for(std::size_t j = 0; j < row_ind.size()-1; j++){
          int row_start = row_ind[j];
          int row_end = row_ind[j+1];
          while(row_start != row_end){
            //if the non-zero element is present at the column index, add it to the column vector  
            if(col_ind[row_start] == col){
            matrix2_column[matrix2_column_pos] = elements[row_start];
            }
            row_start++;
          }
          // increment matrix2_column index to next item
          matrix2_column_pos++;
        }
        // "push back" columns as rows in new elements buffer
        for(int k = 0; k < matrix2_column.size(); k++){
          transpose_elems[transpose_elems_index] = matrix2_column[k];
          transpose_elems_index++;
        }
      }
      //reconstruct transpose matrix
      Matrix A = Matrix(transpose_elems, row_ind.size()-1);

    return A;
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

