#include <algorithm>
#include <limits> // std::numeric_limits

#include "Matrix.hpp"
// your code here...

// Default constructor
Matrix::Matrix() : elements(), col_ind(), row_ind(3), n(2) {}

Matrix::Matrix(std::size_t m, std::size_t n) : elements(), col_ind(),
 row_ind((m == 0 || n == 0) ? std::vector<std::size_t>(1, 0) : std::vector<std::size_t>(m + 1, 0)), n((m == 0) ? 0 : n) {};

Matrix::Matrix(const std::vector<Elem> &A, std::size_t n) : n((A.size() == 0 || n == 0 || A.size() % n != 0) ? 0 : n), row_ind()
{
if (n == 0 || A.size() % n != 0) {
        // Invalid dimensions, set matrix to 0x0
        elements.clear();
        col_ind.clear();
        row_ind.clear();
        row_ind.push_back(0);  // Ensures row_ind has at least one value
        return;
    }

  std::size_t nonzero_elements = 0;
  row_ind.push_back(0);

  // Populate the elements array and corresponding column index.
  for (int i = 0; i < A.size(); i++)
  {
    if (A[i] != 0)
    {
      elements.push_back(A[i]);
      col_ind.push_back(i % n);
      nonzero_elements++;
    }
    // Push back the index of non-zero elements to the row index array when a row finishes
    if ((i + 1) % n == 0)
    {
      row_ind.push_back(nonzero_elements);
    }
  }
};

Matrix::Matrix(const Elem *ptr_A, std::size_t m, std::size_t n) : elements(), col_ind(), n((m == 0) ? 0 : n), row_ind(m + 1, 0) {
  // clear out any entries for invalid dimensions
    if (m == 0 || n == 0) {
        elements.clear();
        col_ind.clear();
        row_ind.clear();
        row_ind.push_back(0);  // Ensures row_ind has at least one value
        return;
    }
    // fill elements and corresponding column index with array from pointer
    for (std::size_t row = 0; row < m; row++) {
        for (std::size_t col = 0; col < n; col++) {
            Elem value = ptr_A[row * n + col];  
            if (value != 0) {
                elements.push_back(value);  
                col_ind.push_back(col);     
            }
        }
        // update row index
        row_ind[row + 1] = elements.size(); 
    }
}

Elem Matrix::e(std::size_t i, std::size_t j) const
{
  // 1. Ensure i is a valid row number
  // 2. Iterate through the column index using the possible entry indices from the row index
  // 3. Look for a column index matching j

  // return minimum value if invalid index
  if (i >= row_ind.size() - 1 || j >= n) {
        return std::numeric_limits<int>::min();
    }

  if (i < row_ind.size() - 1 && i >= 0)
  {
    for (std::size_t index = row_ind[i]; index < row_ind[i + 1]; index++)
    {
      if (col_ind[index] == j)
      {
        return elements[index];
      }
    }
  }
  return 0;
}

bool Matrix::e(std::size_t i, std::size_t j, Elem aij)
{
  // return false for invalid dimensions
    if (j >= n || i >= row_ind.size() - 1) {
        return false; 
    }

    // search for the element in the current row
    for (std::size_t index = row_ind[i]; index < row_ind[i + 1]; index++) {
        if (col_ind[index] == j) {
            if (aij == 0) {
                //If setting to zero, remove the element from elements and corresponding columsn index
                elements.erase(elements.begin() + index);
                col_ind.erase(col_ind.begin() + index);

                //update the row index as well
                for (std::size_t row = i + 1; row < row_ind.size(); row++) {
                    row_ind[row]--;
                }
                return true; 
            } else {
                // set the element at the index
                elements[index] = aij;
                return true;  
            }
        }
    }

    if (aij == 0) {
        return true;  // dont add a 0 element
    }

    // insert new non-zero element
    elements.insert(elements.begin() + row_ind[i + 1], aij);   
    col_ind.insert(col_ind.begin() + row_ind[i + 1], j);     

    for (std::size_t row = i + 1; row < row_ind.size(); row++) {
        row_ind[row]++;  // increment row pointers for subsequent rows
    }

    return true;  
}


std::size_t Matrix::size(std::size_t dim) const
{ // 0x0 matrix
  if (row_ind.size() == 0 && n == 0){
    return 0;
  }
  if (dim == 1)
  { // return rows
    return (row_ind.size() - 1);
  }
  else if (dim == 2)
  { //return columns
    return n;
  }
  else
    return 0;
}

bool Matrix::equal(const Matrix &rhs) const
{
  // Ensure full CSR representation vectors are equal in both matrices
  if ((rhs.elements == elements) && (rhs.col_ind == col_ind) && (rhs.row_ind == row_ind))
  {
    return true;
  }
  else
    return false;
}

Matrix Matrix::add(const Matrix &rhs) const
{
    // Check if dimensions match
    if (row_ind.size() != rhs.row_ind.size() || n != rhs.n)
    {
        return Matrix(0, 0); // Inconsistent dimensions
    }

    // Initialize the result matrix 
    Matrix summed_matrix(row_ind.size() - 1, n);
    
    // ensure correct dim for row index
    summed_matrix.row_ind.resize(row_ind.size());

    // Iterate through each row of the matrix
    for (std::size_t i = 0; i < row_ind.size() - 1; ++i)
    {
        std::size_t matrix1_row_start = row_ind[i];
        std::size_t matrix1_row_end = row_ind[i + 1];

        std::size_t matrix2_row_start = rhs.row_ind[i];
        std::size_t matrix2_row_end = rhs.row_ind[i + 1];

        // Add non-zero elements from both matrices
        while (matrix1_row_start < matrix1_row_end || matrix2_row_start < matrix2_row_end)
        {
            if (matrix1_row_start == matrix1_row_end)
            {
                // Matrix 1 is exhausted, take from Matrix 2
                summed_matrix.elements.push_back(rhs.elements[matrix2_row_start]);
                summed_matrix.col_ind.push_back(rhs.col_ind[matrix2_row_start]);
                matrix2_row_start++;
            }
            else if (matrix2_row_start == matrix2_row_end)
            {
                // Matrix 2 is exhausted, take from Matrix 1
                summed_matrix.elements.push_back(elements[matrix1_row_start]);
                summed_matrix.col_ind.push_back(col_ind[matrix1_row_start]);
                matrix1_row_start++;
            }
            else if (col_ind[matrix1_row_start] == rhs.col_ind[matrix2_row_start])
            {
                // Columns match, add the elements
                Elem sum = elements[matrix1_row_start] + rhs.elements[matrix2_row_start];
                if (sum != 0)  // only insert non-zero
                {
                    summed_matrix.elements.push_back(sum);
                    summed_matrix.col_ind.push_back(col_ind[matrix1_row_start]);
                }
                matrix1_row_start++;
                matrix2_row_start++;
            }
            else if (col_ind[matrix1_row_start] < rhs.col_ind[matrix2_row_start])
            {
                // Column of Matrix 1 is smaller, take the element from Matrix 1
                summed_matrix.elements.push_back(elements[matrix1_row_start]);
                summed_matrix.col_ind.push_back(col_ind[matrix1_row_start]);
                matrix1_row_start++;
            }
            else
            {
                // Column of Matrix 2 is smaller, take the element from Matrix 2
                summed_matrix.elements.push_back(rhs.elements[matrix2_row_start]);
                summed_matrix.col_ind.push_back(rhs.col_ind[matrix2_row_start]);
                matrix2_row_start++;
            }
        }

        // Update the row index for the summed matrix
        summed_matrix.row_ind[i + 1] = summed_matrix.elements.size();
    }

    return summed_matrix;
}

Matrix Matrix::sub(const Matrix &rhs) const
{
    // Check if dimensions match
    if (row_ind.size() != rhs.row_ind.size() || n != rhs.n)
    {
        return Matrix(0, 0); // Inconsistent dimensions
    }

    // Initialize the result matrix 
    Matrix diff_matrix(row_ind.size() - 1, n);
    
    // ensure correct dim for row index
    diff_matrix.row_ind.resize(row_ind.size());

    // Iterate through each row of the matrix
    for (std::size_t i = 0; i < row_ind.size() - 1; ++i)
    {
        std::size_t matrix1_row_start = row_ind[i];
        std::size_t matrix1_row_end = row_ind[i + 1];

        std::size_t matrix2_row_start = rhs.row_ind[i];
        std::size_t matrix2_row_end = rhs.row_ind[i + 1];

        // Subtract non-zero elements from both matrices
        while (matrix1_row_start < matrix1_row_end || matrix2_row_start < matrix2_row_end)
        {
            if (matrix1_row_start == matrix1_row_end)
            {
                // Matrix 1 is exhausted, take from Matrix 2
                diff_matrix.elements.push_back(-rhs.elements[matrix2_row_start]);
                diff_matrix.col_ind.push_back(rhs.col_ind[matrix2_row_start]);
                matrix2_row_start++;
            }
            else if (matrix2_row_start == matrix2_row_end)
            {
                // Matrix 2 is exhausted, take from Matrix 1
                diff_matrix.elements.push_back(elements[matrix1_row_start]);
                diff_matrix.col_ind.push_back(col_ind[matrix1_row_start]);
                matrix1_row_start++;
            }
            else if (col_ind[matrix1_row_start] == rhs.col_ind[matrix2_row_start])
            {
                // Columns match, subtract the elements
                Elem diff = elements[matrix1_row_start] - rhs.elements[matrix2_row_start];
                if (diff != 0)  // only non-zero elements
                {
                    diff_matrix.elements.push_back(diff);
                    diff_matrix.col_ind.push_back(col_ind[matrix1_row_start]);
                }
                matrix1_row_start++;
                matrix2_row_start++;
            }
            else if (col_ind[matrix1_row_start] < rhs.col_ind[matrix2_row_start])
            {
                // Column of Matrix 1 is smaller, take the element from Matrix 1
                diff_matrix.elements.push_back(elements[matrix1_row_start]);
                diff_matrix.col_ind.push_back(col_ind[matrix1_row_start]);
                matrix1_row_start++;
            }
            else
            {
                // Column of Matrix 2 is smaller, take the negated element from Matrix 2
                diff_matrix.elements.push_back(-rhs.elements[matrix2_row_start]);
                diff_matrix.col_ind.push_back(rhs.col_ind[matrix2_row_start]);
                matrix2_row_start++;
            }
        }

        // Update the row index for the diff matrix
        diff_matrix.row_ind[i + 1] = diff_matrix.elements.size();
    }

    return diff_matrix;
}

Matrix Matrix::mult(const Matrix &rhs) const
{
  if (n != rhs.row_ind.size()-1){
    Matrix A = Matrix(0,0);
    return A;
  }

  // resultant matrix in Row-Major order
  std::vector<Elem> mult_result((row_ind.size() - 1) * rhs.n, 0);
  std::size_t result_pointer = 0;
  for (std::size_t i = 0; i < row_ind.size() - 1; i++)
  {

    int row_start = row_ind[i];
    int row_end = row_ind[i + 1];

    // buffers for multiplication operation
    std::vector<int> matrix1_row(n, 0);

    int column_index = 0;

    // update corresponding values in row buffer based on CSR
    while (row_start != row_end)
    {
      if (col_ind[row_start] == column_index)
      {
        matrix1_row[column_index] = elements[row_start];
        row_start++;
      }
      column_index++;
    }
    // construct a vector from each column in right hand side
    for (std::size_t col = 0; col < rhs.n; col++)
    {

      std::vector<int> matrix2_column(rhs.row_ind.size() - 1, 0);
      int matrix2_column_pos = 0;

      // build vector of columns
      for (std::size_t j = 0; j < rhs.row_ind.size() - 1; j++)
      {
        int rhs_row_start = rhs.row_ind[j];
        int rhs_row_end = rhs.row_ind[j + 1];
        while (rhs_row_start != rhs_row_end)
        {
          // if the non-zero element is present at the column index, add it to the column vector
          if (rhs.col_ind[rhs_row_start] == col)
          {
            matrix2_column[matrix2_column_pos] = rhs.elements[rhs_row_start];
          }
          rhs_row_start++;
        }
        // increment matrix2_column index to next item
        matrix2_column_pos++;
      }
      
      // perform the vector multiplication on each index and add them
      std::size_t total = 0;
      for (std::size_t item = 0; item < n; item++)
      {
        total += matrix1_row[item] * matrix2_column[item];
      } 
      //set the correct index of the resultant vector with the total (row-major order)
      mult_result[result_pointer] = total;
      result_pointer++;
    }
  }

  return Matrix(mult_result, rhs.n);
}

Matrix Matrix::mult(Scalar k) const
{
  // row major order vector buffer
  std::vector<Elem> new_elements(elements.size());
  // copy the existing matrix and scale non-zero elements.
  for (std::size_t entry = 0; entry < elements.size(); entry++)
  {
    new_elements[entry] = elements[entry] * k;
  }
  // construct the new matrix
  Matrix A(*this);
  // copy over col ind and row ind
  A.elements = new_elements;
  A.col_ind = col_ind;
  A.row_ind = row_ind;

  return A;
}

Matrix Matrix::pow(Scalar k) const
{
  // return a 0x0 matrix for invalid exponent
  if (row_ind.size()-1 != n || k < 0){
    return Matrix(0,0);
  }
  // copy the matrix
  Matrix A(*this);
  // copy fields
  A.col_ind = col_ind;
  A.row_ind = row_ind;

    //create the identity matrix if k =0 0
    if (k == 0)
    {   
        Matrix identity(row_ind.size() - 1, n);  
        for (std::size_t i = 0; i < row_ind.size() - 1; ++i)
        {
            identity.elements.push_back(1);     
            identity.col_ind.push_back(i);    
            identity.row_ind[i + 1] = i + 1;     
        }
        return identity;
    }

  // else run the multiplication operation for as many times as k
  for (int i = 0; i < k - 1; i++)
  {
    A = A.mult(*this);
  }

  return A;
}

Matrix Matrix::trans() const
{

  // return empty matrix with same dim as original if matrix is all zero
  if (elements.size() == 0)
  {
    return *this;
  }

  // element buffer used for matrix reconstruction
  std::vector<Elem> transpose_elems(elements.size());
  int transpose_elems_index = 0;

  // grab each column out of the matrix
  for (std::size_t col = 0; col < n; col++)
  {

    // column buffer
    std::vector<int> matrix2_column(row_ind.size() - 1, 0);
    int matrix2_column_pos = 0;

    // build vector of columns
    for (std::size_t j = 0; j < row_ind.size() - 1; j++)
    {
      int row_start = row_ind[j];
      int row_end = row_ind[j + 1];
      while (row_start != row_end)
      {
        // if the non-zero element is present at the column index, add it to the column vector
        if (col_ind[row_start] == col)
        {
          matrix2_column[matrix2_column_pos] = elements[row_start];
        }
        row_start++;
      }
      // increment matrix2_column index to next item
      matrix2_column_pos++;
    }
    // "push back" columns as rows in new elements buffer
    for (int k = 0; k < matrix2_column.size(); k++)
    {
      transpose_elems[transpose_elems_index] = matrix2_column[k];
      transpose_elems_index++;
    }
  }
  // reconstruct transpose matrix
  Matrix A = Matrix(transpose_elems, row_ind.size() - 1);

  return A;
}



bool Matrix::rowSwitch(std::size_t i, std::size_t j) {

    // Ensure valid dims and valid row swap references
    if (row_ind.size() <= 2 || elements.size() <= 1 || i >= row_ind.size() - 1 || j >= row_ind.size() - 1) {
        return false;  // Invalid input or no rows to swap
    }

    // Start and end of each row
    int row1_start = row_ind[i];   
    int row1_end = row_ind[i + 1]; 
    int row2_start = row_ind[j];   
    int row2_end = row_ind[j + 1]; 

    if (row1_start == row1_end && row2_start == row2_end) {
        return false;  // Nothing to swap, both rows are empty
    }

    // Swap elements and column indices if rows have the same size
    std::swap_ranges(elements.begin() + row1_start, elements.begin() + row1_end, elements.begin() + row2_start);
    std::swap_ranges(col_ind.begin() + row1_start, col_ind.begin() + row1_end, col_ind.begin() + row2_start);

    // Adjust row sizes if both are different
    int row1_size = row1_end - row1_start;
    int row2_size = row2_end - row2_start;

    if (row1_size != row2_size) {
        row_ind[i + 1] = row_ind[i] + row2_size; 
        row_ind[j + 1] = row_ind[j] + row1_size; // Update row_ind for j
    }

    return true;  
}

bool Matrix::rowMult(std::size_t i, Scalar k)
{
    if (i >= row_ind.size() - 1)  // Ensure i is within bounds
    {
        return false;
    }

    // Start and end of row i
    int row_start = row_ind[i];
    int row_end = row_ind[i + 1];

    // no rows to multiply
    if (elements.size() == 0)
    {
        return false;
    }

    if (k == 0)
    {
        // Remove all elements in this row since k == 0
        elements.erase(elements.begin() + row_start, elements.begin() + row_end);
        col_ind.erase(col_ind.begin() + row_start, col_ind.begin() + row_end);

        int num_removed_elements = row_end - row_start;

        // Adjust row_ind for subsequent rows
        for (std::size_t j = i + 1; j < row_ind.size(); j++)
        {
            row_ind[j] -= num_removed_elements;
        }

        // Remove the current row from row_ind
        row_ind.erase(row_ind.begin() + i);
        return true;
    }

    // Multiply all elements in the row by k
    for (int index = row_start; index < row_end; index++)
    {
        elements[index] *= k;
    }

    return true;
}

bool Matrix::rowAdd(std::size_t i, std::size_t j, Scalar k) {

    // return false for invalid dimensions or indices
    if (i < 0 || i > row_ind.size()-1 || j < 0 || j > row_ind.size()-1){
      return false;
    }
    if(row_ind.size() <= 2 || elements.size() <= 1) {
        return false;
    }

    // start and end indices for both rows
    int dest_row_start = row_ind[i];  
    int dest_row_end = row_ind[i + 1];       

    int row_start = row_ind[j];      
    int row_end = row_ind[j + 1];            

    std::vector<Elem> dest_row(n, 0);  
    std::vector<Elem> src_row(n, 0);   

    // Fill the destination row 
    for (int l = dest_row_start; l < dest_row_end; l++) {
        dest_row[col_ind[l]] = elements[l];
    }

    // Multiply the source row by scalar k and fill
    for (int m = row_start; m < row_end; m++) {
        src_row[col_ind[m]] = elements[m] * k;
    }

    // Add the source row to the destination row
    for (int o = 0; o < n; o++) {
        dest_row[o] += src_row[o];
    }

    // Erase the old destination row from the matrix and insert the updated destination row
    elements.erase(elements.begin() + dest_row_start, elements.begin() + dest_row_end);
    col_ind.erase(col_ind.begin() + dest_row_start, col_ind.begin() + dest_row_end);

    int insert_pos = dest_row_start;
    for (std::size_t col = 0; col < n; col++) {
        if (dest_row[col] != 0) {  // Only add non-zero elements
            elements.insert(elements.begin() + insert_pos, dest_row[col]);
            col_ind.insert(col_ind.begin() + insert_pos, col);
            insert_pos++;
        }
    }

    // Adjust row_ind if the row size changed
    int new_row_size = insert_pos - dest_row_start;
    int old_row_size = dest_row_end - dest_row_start;
    int size_difference = new_row_size - old_row_size;

    for (std::size_t r = i + 1; r < row_ind.size(); r++) {  
        row_ind[r] += size_difference;
    }

    return true; 
}

Matrix Matrix::cat(const Matrix &rhs, std::size_t dim) const
{
    // Vertical concatenation
    if (dim == 1 && n == rhs.n) // Check for same num of columns
    {
        std::vector<Elem> row_major; 

        // Copy elements from the current matrix to row-major format
        for (std::size_t i = 0; i < row_ind.size() - 1; ++i)
        {
            std::vector<Elem> row(n, 0); // row buffer
            for (std::size_t j = row_ind[i]; j < row_ind[i + 1]; ++j)
            {
                row[col_ind[j]] = elements[j];
            }
            row_major.insert(row_major.end(), row.begin(), row.end());  // Append the row to row_major
        }

        // Copy elements from the rhs matrix to row-major format
        for (std::size_t i = 0; i < rhs.row_ind.size() - 1; ++i)
        {
            std::vector<Elem> row(rhs.n, 0);  // row buffer
            for (std::size_t j = rhs.row_ind[i]; j < rhs.row_ind[i + 1]; ++j)
            {
                row[rhs.col_ind[j]] = rhs.elements[j];
            }
            row_major.insert(row_major.end(), row.begin(), row.end());  // Append the row to row_major
        }

        // Construct the final matrix 
        return Matrix(row_major, n);  
    }
    // Horizontal concatenation
    else if (dim == 2 && row_ind.size() == rhs.row_ind.size())  // Check if both matrices have the same number of rows
    {
        std::vector<Elem> row_major;

        // Append each row from both matrices side by side
        for (std::size_t i = 0; i < row_ind.size() - 1; ++i)
        {
            std::vector<Elem> row(n + rhs.n, 0);  // Zero-filled row for both matrices

            // Fill the row with elements from the first matrix
            for (std::size_t j = row_ind[i]; j < row_ind[i + 1]; ++j)
            {
                row[col_ind[j]] = elements[j];
            }

            // Fill the row with elements from the rhs matrix, adjust for column offset
            for (std::size_t j = rhs.row_ind[i]; j < rhs.row_ind[i + 1]; ++j)
            {
                row[rhs.col_ind[j] + n] = rhs.elements[j];  // Shift the column index by n
            }

            row_major.insert(row_major.end(), row.begin(), row.end());  // Append the row to row_major
        }

        // Construct the final matrix with n + rhs.n columns
        return Matrix(row_major, n + rhs.n);  
    }
    else
    {
        // Return a 0x0 matrix for invalid dimensions
        return Matrix(0, 0);
    }
}






// You'll implement this method in part two of the project
Matrix Matrix::invMod(Scalar m) const
{
  Matrix Identity(row_ind.size()-1, row_ind.size()-1);
  for (std::size_t i = 0; i < m; ++i) {
        Identity.elements.push_back(1);   // Elements on the diagnol are 1
        Identity.col_ind.push_back(i);    // Column index is the same as the row index
        Identity.row_ind[i + 1] = i + 1;  //  Increment by 1 after every row
    }

  Matrix AugmentedMatrix = this->cat(Identity, 2);
  std::cout << AugmentedMatrix << std::endl;

  for (std::size_t i = 0; i < AugmentedMatrix.row_ind.size() - 1; ++i) {
    std::size_t pivot_index = -1;
    // Step 1: Find the pivot element in the current row
    for (int j = AugmentedMatrix.row_ind[i]; j < AugmentedMatrix.row_ind[i+1]; j++) {
        if (AugmentedMatrix.elements[j] != 0 && AugmentedMatrix.col_ind[j] == i) {  // Check column index matches
            pivot_index = j;
            std::cout << "Pivot found at index: " << pivot_index << " in column " << i << std::endl;
            break;
        }
    }
    Elem pivot_value = AugmentedMatrix.elements[pivot_index];  // The pivot element
    Elem pivot_inverse = AugmentedMatrix.modularInverse(pivot_value, m); 
    
    if (pivot_value == 0 || pivot_inverse == -1) {
    std::cout << "Pivot is zero or has no modular inverse, matrix is not invertible!" << std::endl;
    return Matrix(0, 0);  // Return an empty matrix to indicate non-invertibility
}

    AugmentedMatrix.rowMult(i,pivot_inverse);

    for (int j = AugmentedMatrix.row_ind[i]; j < AugmentedMatrix.row_ind[i+1]; j++) {
        AugmentedMatrix.elements[j] = AugmentedMatrix.elements[j] % m;
    }

    std::cout << AugmentedMatrix << std::endl;

   // Perform elimination below the pivot row
for (std::size_t row = i + 1; row < AugmentedMatrix.row_ind.size() - 1; ++row) {
    std::size_t row_start = AugmentedMatrix.row_ind[row];
    std::size_t row_end = AugmentedMatrix.row_ind[row + 1];

    Elem element_below_pivot = -1;

    // Find the element in the same column as the pivot (pivot column = i)
    for (std::size_t j = row_start; j < row_end; ++j) {
        if (AugmentedMatrix.col_ind[j] == i) {  // Find element in the pivot column
            element_below_pivot = AugmentedMatrix.elements[j];
            break;
        }
    }

    // If there's a non-zero element below the pivot, perform row elimination
    if (element_below_pivot != -1 && element_below_pivot != 0) {
        Elem factor = element_below_pivot;

        // Subtract factor * pivot row from the current row using rowAdd
        AugmentedMatrix.rowAdd(row, i, (-factor % m + m) % m);  // Eliminate below the pivot

        // Ensure all elements in the current row are reduced modulo m
        std::size_t new_row_start = AugmentedMatrix.row_ind[row];
        std::size_t new_row_end = AugmentedMatrix.row_ind[row + 1];

        for (std::size_t j = new_row_start; j < new_row_end; ++j) {
            AugmentedMatrix.elements[j] = (AugmentedMatrix.elements[j] % m + m) % m;
        }

        // Double check row after elimination and mod operation
        std::cout << "Row " << row << " after row elimination and mod m: " << AugmentedMatrix << std::endl;
    }
}
// Eliminate elements above the pivot
for (std::size_t row = i; row-- > 0;) {  // Work upwards
    std::size_t row_start = row_ind[row];
    std::size_t row_end = row_ind[row + 1];

    Elem element_above_pivot = -1;

    // Find the element in the same column as the pivot (pivot column = i)
    for (std::size_t j = row_start; j < row_end; ++j) {
        if (AugmentedMatrix.col_ind[j] == i) {  // Find element in the pivot column
            element_above_pivot = AugmentedMatrix.elements[j];
            break;
        }
    }

    // If there's a non-zero element above the pivot, perform row elimination
    if (element_above_pivot != -1 && element_above_pivot != 0) {
        Elem factor = element_above_pivot;

        // Subtract factor * pivot row from the current row using rowAdd
        std::cout << "Performing rowAdd on row " << row << " with factor: " << factor << std::endl;
        AugmentedMatrix.rowAdd(row, i, (-factor % m + m) % m);  // Eliminate above the pivot

        // Ensure all elements in the current row are reduced modulo m
        std::size_t new_row_start = AugmentedMatrix.row_ind[row];
        std::size_t new_row_end = AugmentedMatrix.row_ind[row + 1];

        for (std::size_t j = new_row_start; j < new_row_end; ++j) {
            AugmentedMatrix.elements[j] = (AugmentedMatrix.elements[j] % m + m) % m;
        }

        // Output the current state of the row after elimination
        std::cout << "Row " << row << " after row elimination and mod m: " << AugmentedMatrix << std::endl;
    }
}
  }
  Matrix inverse = AugmentedMatrix.extractInverse(n);  // n is the size of the original matrix
std::cout << inverse << std::endl;
  return Matrix();
}

int Matrix::extendedEuclidean(int a, int m, int &x, int &y) {
    if (a == 0) {
        x = 0;
        y = 1;
        return m;
    }

    int x1, y1;
    int gcd = extendedEuclidean(m % a, a, x1, y1);

    x = y1 - (m / a) * x1;
    y = x1;

    return gcd;
}

int Matrix::modularInverse(int a, int m) {
    int x, y;
    int g = extendedEuclidean(a, m, x, y);
    if (g != 1) {
        return -1;  
    } else {
        // Ensure x is positive
        return (x % m + m) % m;
    }
}

Matrix Matrix::extractInverse(std::size_t n) const {
    Matrix inverseMatrix(n, n);  // Create a new matrix to hold the inverse
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t row_start = row_ind[i];
        std::size_t row_end = row_ind[i + 1];

        // Loop through each row in the augmented matrix
        for (std::size_t j = row_start; j < row_end; ++j) {
            // Check if the column is in the range of the right-hand side
            if (col_ind[j] >= n && col_ind[j] < 2 * n) {
                // Adjust the column index to fit into the inverse matrix
                std::size_t inverse_col = col_ind[j] - n;
                inverseMatrix.elements.push_back(elements[j]);
                inverseMatrix.col_ind.push_back(inverse_col);
            }
        }
        // Update the row index for the inverse matrix
        inverseMatrix.row_ind[i + 1] = inverseMatrix.elements.size();
    }
    return inverseMatrix;
}

// Instructor provided methods and functions below; do not modify
Matrix::Matrix(const Matrix &A) : n(A.n)
{
  elements = A.elements;
  col_ind = A.col_ind;
  row_ind = A.row_ind;
}

Matrix &Matrix::operator=(Matrix A)
{
  std::swap(elements, A.elements);
  std::swap(col_ind, A.col_ind);
  std::swap(row_ind, A.row_ind);

  return *this;
}

std::ostream &operator<<(std::ostream &os, const Matrix &A)
{
  // Loop through each row
  for (std::size_t i = 0; i < A.row_ind.size() - 1; i++)
  {
    std::size_t row_start = A.row_ind[i];   // Start of the current row's non-zero elements
    std::size_t row_end = A.row_ind[i + 1]; // End of the current row's non-zero elements
    std::size_t column_pointer = 0;         // Column index tracker

    // Execute while in scope of row
    while (column_pointer < A.n)
    {

      if (row_start < row_end && A.col_ind[row_start] == column_pointer)
      {
        // if there is a non-zero entry in the specified column, print the entry
        os << A.elements[row_start] << " ";
        row_start++;
      }
      else
      {
        // if not print a 0
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
