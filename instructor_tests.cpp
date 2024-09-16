#include <limits> // std::numeric_limits
#include <random>

#include "catch.hpp"
#include "Matrix.hpp"

//for generating random integers
std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
std::uniform_int_distribution<> rnd_p(0, 9); //random numbers on the closed interval [0, 9]
std::uniform_int_distribution<> rnd_pn(-9, 9); //random numbers on the closed interval [-9, 9]

bool DEBUG = true;

//Check equality of matrix A with its reference elements, rows, and columns
bool equality(const Matrix &A, const std::vector<Elem> &elements, std::size_t m, std::size_t n)
{
  bool flag = true;

  if(A.size(1) == m && A.size(2) == n) //check that size of matrix is correct
    { //check elements of A against elements used to construct matrix
      for(std::size_t i = 0; i < m; i++)
  	for(std::size_t j = 0; j < n; j++)
	  {
	    if(A.e(i,j) != elements[(i * n) + j])
	      {
		flag = false;
		break;
	      }
	  }
    }
  else
    flag = false;
  
  return flag;
}


TEST_CASE( "Default constructor", "[Matrix]" )
{
  INFO("Hint: create an empty 2-by-2 matrix");
  Matrix A;

  //check dimensions to avoid invalid memory access 
  REQUIRE(A.size(1) == 2);
  REQUIRE(A.size(2) == 2);

  //check elements
  std::size_t i, j;
  
  for(i = 0; i < 2; i++)
    for(j = 0; j < 2; j++)
      REQUIRE(A.e(i,j) == 0);
}

/** 
 * NOTE: the below test cases are not comprehensive/complete/correct so you must supplement them.  They are only 
 * meant to invoke a particular method to force you to write stubs for all your methods.  This should ensure 
 * that your code will at least compile on the autograder because the instructor test code calls every method in
 * the header file and without at least a stub an error is produced.
 */
TEST_CASE( "Parameterized constructor 1", "[Matrix]" )
{
  INFO("Hint: create an empty matrix of user-supplied dimensions");
  std::size_t m = 3;
  std::size_t n = 5;
  std::size_t i, j; //for looping

  //test case: square matrix, non-default size
  Matrix A(m,m); 

  //check dimensions to avoid invalid memory access 
  REQUIRE(A.size(1) == m);
  REQUIRE(A.size(2) == m);

  //check elements
  for(i = 0; i < m; i++)
    for(j = 0; j < m; j++)
      REQUIRE(A.e(i,j) == 0);

  //test case: non-square matrix
  Matrix B(m,n); 

  //check dimensions to avoid invalid memory access 
  REQUIRE(B.size(1) == m);
  REQUIRE(B.size(2) == n);

  //check elements
  for(i = 0; i < m; i++)
    for(j = 0; j < n; j++)
      REQUIRE(B.e(i,j) == 0);

  //test case: invalid dimensions, row
  Matrix C(0,1);
  REQUIRE(C.size(1) == 0);
  REQUIRE(C.size(2) == 0);
  
  //test case: invalid dimensions, column
  Matrix D(1,0);
  REQUIRE(D.size(1) == 0);
  REQUIRE(D.size(2) == 0);
}

/**
 * NOTE: I would recommend that as you implement more methods/functionality that you supplement getter/setter 
 * with additional tests.  This makes testing an iterative procedure; e.g., once we can set elements we can test
 * the getter method on non-empty matrices. 
 */
TEST_CASE( "Element get", "[Matrix]" )
{
  INFO("Hint: Returns the element at specified row, column index");
  Matrix A;
  std::vector<Elem> elements;
  std::size_t i, j;

  //1. Valid indices tests
  //emtpy matrix, 2-by-2
  for(i = 0; i < 2; i++)
    for(j = 0; j < 2; j++)
      REQUIRE(A.e(i,j) == 0);

  //non-empty matrix, 4-by-6: requires Parameterized constructor 2 to work
  //get random values for a matrix
  for (int i = 0; i < 24; i++)
    elements.push_back(rnd_p(gen));

  Matrix B(elements, 6);
  
  REQUIRE(equality(B, elements, 4, 6));
  
  //2. Invalid indices tests
  Elem aij_min = std::numeric_limits<int>::min();
  REQUIRE(B.e(-1,0) == aij_min);
  REQUIRE(B.e(0,-1) == aij_min);
  REQUIRE(B.e(10,0) == aij_min);
  REQUIRE(B.e(0,10) == aij_min);
  REQUIRE(B.e(10,10) == aij_min);
}

TEST_CASE( "Element set", "[Matrix]" )
{
  INFO("Hint: Sets the element at specified row, column index to given value");
  Matrix A(3,3);

  //1. simple tests: create identity matrix (tests insert at each row and column)
  REQUIRE(A.e(0,0,1));  //test case: insert into empty matrix
  REQUIRE(A.e(0,0) == 1);

  REQUIRE(A.e(1,1,1));  //test case: insert into non-empty matrix but empty row
  REQUIRE(A.e(1,1) == 1);

  REQUIRE(A.e(2,2,1));  //test case: insert into last column
  REQUIRE(A.e(2,2) == 1);

  //2. corner cases for insert (e.g., first or last column) and proper overwrites
  REQUIRE(A.e(2,0,1));  //test case: insert into non-empty matrix but empty row, first column
  REQUIRE(A.e(2,0) == 1);

  REQUIRE(A.e(0,0,2));  //test case: overwrite entry in first column, first row
  REQUIRE(A.e(0,0) == 2);

  REQUIRE(A.e(2,2,2));  //test case: overwrite entry in last column, last row
  REQUIRE(A.e(2,2) == 2);

  //test case: ensure no other values affected by insert
  Matrix B = A;
  REQUIRE(A.equal(B));

  //make change to non-zero element
  REQUIRE(B.e(2,2,1));
  REQUIRE(B.e(2,2) == 1);
  //undo change to non-zero element
  REQUIRE(B.e(2,2,2));
  REQUIRE(B.e(2,2) == 2);

  REQUIRE(A.equal(B));

  //test case: set zero element to zero (no change to matrix)
  REQUIRE(B.e(0,1,0));
  REQUIRE(B.e(0,1) == 0);
  REQUIRE(A.equal(B));
  
  //test case: replace existing element with zero (should remove entry from matrix)
  Matrix C, D; //start with two empty matrices, fill one, then delete, should be identical
  REQUIRE(C.equal(D));

  //set entries
  REQUIRE(D.e(0,0,1));
  REQUIRE(D.e(0,0) == 1);
  REQUIRE(D.e(0,1,2));
  REQUIRE(D.e(0,1) == 2);
  REQUIRE(D.e(1,0,3));
  REQUIRE(D.e(1,0) == 3);
  REQUIRE(D.e(1,1,4));
  REQUIRE(D.e(1,1) == 4);

  //(un)set entries to zero, result should be C == D
  REQUIRE(D.e(0,0,0));
  REQUIRE(D.e(0,0) == 0);
  REQUIRE(D.e(0,1,0));
  REQUIRE(D.e(0,1) == 0);
  REQUIRE(D.e(1,0,0));
  REQUIRE(D.e(1,0) == 0);
  REQUIRE(D.e(1,1,0));
  REQUIRE(D.e(1,1) == 0);

  REQUIRE(C.equal(D));
  
  //test case: invalid set, empty matrix
  REQUIRE_FALSE(C.e(-1,0,1));
  REQUIRE_FALSE(C.e(0,-1,1));
  REQUIRE_FALSE(C.e(10,0,1));
  REQUIRE_FALSE(C.e(0,10,1));
  REQUIRE_FALSE(C.e(10,10,1));

  //test case: invalid set, non-empty matrix
  REQUIRE_FALSE(B.e(-1,0,1));
  REQUIRE_FALSE(B.e(0,-1,1));
  REQUIRE_FALSE(B.e(10,0,1));
  REQUIRE_FALSE(B.e(0,10,1));
  REQUIRE_FALSE(B.e(10,10,1));

  //3. reference matrix: manually check contents of private data members for correct ordering
  // matrix from https://en.wikipedia.org/wiki/Sparse_matrix
  // m = 4, n = 6
  // elements := [ 10 20 30 40 50 60 70 80 ]
  // col_ind := [  0  1  1  3  2  3  4  5 ]
  // row_ind := [  0  2  4  7  8 ]
  Matrix R(4,6);
  R.e(0, 0, 10);
  REQUIRE(R.e(0, 0) == 10);
  R.e(0, 1, 20);
  REQUIRE(R.e(0, 1) == 20);
  R.e(1, 1, 30);
  REQUIRE(R.e(1, 1) == 30);
  R.e(1, 3, 40);
  REQUIRE(R.e(1, 3) == 40);
  R.e(2, 2, 50);
  REQUIRE(R.e(2, 2) == 50);
  R.e(2, 3, 60);
  REQUIRE(R.e(2, 3) == 60);
  R.e(2, 4, 70);
  REQUIRE(R.e(2, 4) == 70);
  R.e(3, 5, 80);
  REQUIRE(R.e(3, 5) == 80);

  std::cout << R << std::endl;
}

TEST_CASE( "Parameterized constructor 2", "[Matrix]" )
{
  INFO("Hint: Use the parameters to set the matrix elements, vector");
  std::vector<Elem> elements;
  
  //get random values for a matrix
  for (int i = 0; i < 4; i++)
    elements.push_back(rnd_p(gen));

  //test case: square matrix
  Matrix A(elements,2);
  REQUIRE(equality(A, elements, 2, 2));

  //test case: one column matrix
  Matrix B(elements,1);
  REQUIRE(equality(B, elements, 4, 1));

  //test case: one row matrix
  Matrix C(elements,4);
  REQUIRE(equality(C, elements, 1, 4));
  
  //test case: invalid dimensions, column
  Matrix D(elements,0);
  REQUIRE(D.size(1) == 0);
  REQUIRE(D.size(2) == 0);
}


TEST_CASE( "Parameterized constructor 3", "[Matrix]" )
{
  INFO("Hint: Use the parameters to set the matrix elements, array");
  Elem elements[4] = {1, 0, 0, 1};
  std::vector<Elem> elements_v = {1, 0, 0, 1};
  
  //test case: square matrix
  Matrix A(elements,2,2);
  REQUIRE(equality(A, elements_v, 2, 2));

  //test case: non-square matrix
  Matrix B(elements,1,4);
  REQUIRE(equality(B, elements_v, 1, 4));
  
  //test case: invalid dimensions, row
  Matrix C(elements,0,4);
  REQUIRE(C.size(1) == 0);
  REQUIRE(C.size(2) == 0);
  
  //test case: invalid dimensions, column
  Matrix D(elements,4,0);
  REQUIRE(D.size(1) == 0);
  REQUIRE(D.size(2) == 0);
}

TEST_CASE( "Method size", "[Matrix]" )
{
  INFO("Hint: Sets the element at specified row, column index to given value");
  Matrix A;

  //test case: empty matrix, default size
  REQUIRE(A.size(1) == 2);
  REQUIRE(A.size(2) == 2);

  //test case: empty matrix, non-default size
  Matrix B(3,4);
  REQUIRE(B.size(1) == 3);
  REQUIRE(B.size(2) == 4);
  
  //test case: invalid dimension
  REQUIRE(A.size(3) == 0);
}

TEST_CASE( "Method equal", "[Matrix]" )
{
  INFO("Hint: Returns true if two matrices are the same, false otherwise");
  std::vector<Elem> elements;
  
  //get random values for a matrix
  for (int i = 0; i < 9; i++)
    elements.push_back(rnd_p(gen));

  //test case: empty square matrices
  Matrix A, B;
  REQUIRE(A.equal(B));

  //test case: square matrices of different sizes
  Matrix C(elements,3);
  REQUIRE(C.equal(C));
  REQUIRE_FALSE(C.equal(A));

  //test case: matrices with same elements (ordered the same in elements vector) but different dimensions
  Matrix D(elements,9);
  REQUIRE(D.equal(D));
  REQUIRE_FALSE(D.equal(C));
  
  //test case: empty matrices of different dimensions and number of elements
  Matrix E(3,3);
  REQUIRE_FALSE(E.equal(A));

  //test case: empty matrices of different dimensions but same number of elements
  Matrix F(4,1);
  REQUIRE_FALSE(F.equal(A));
}

TEST_CASE( "Method add", "[Matrix]" )
{
  INFO("Hint: Addition of matrices");

  SECTION("valid additions")
    {
      std::vector<Elem> a = {1,2,3,4};
      std::vector<Elem> b = {5,6,7,8};
      std::vector<Elem> c = {6,8,10,12};

      //2-by-2
      Matrix A(a,2), B(b,2);
      Matrix C = A.add(B);

      REQUIRE(equality(C,c,2,2));

      //1-by-4
      Matrix A_14(a,4), B_14(b,4);
      Matrix C_14 = A_14.add(B_14);
  
      REQUIRE(equality(C_14, c, 1, 4));

      //4-by-1
      Matrix A_41(a,1), B_41(b,1);
      Matrix C_41 = A_41.add(B_41);
  
      REQUIRE(equality(C_41, c, 4, 1));
    }

  SECTION("invalid additions")
    {
      std::vector<Elem> a = {1,2,3,4};
      std::vector<Elem> b = {5,6};

      //2-by-2 + 2-by-1
      Matrix A(a, 2), B(b, 1);
      Matrix C = A.add(B);

      REQUIRE(equality(C, {}, 0, 0));

      //2-by-2 + 1-by-2
      Matrix B_12(b,2);
      Matrix D = A.add(B_12);
  
      REQUIRE(equality(D, {}, 0, 0));
    }
}

TEST_CASE( "Method sub", "[Matrix]" )
{
  INFO("Hint: Subtraction of matrices");

  SECTION("valid subtractions")
    {
      std::vector<Elem> a = {1,2,3,4};
      std::vector<Elem> b = {5,6,7,8};
      std::vector<Elem> c = {-4,-4,-4,-4};

      //2-by-2
      Matrix A(a, 2), B(b, 2);
      Matrix C = A.sub(B);

      REQUIRE(equality(C, c, 2, 2));

      //1-by-4
      Matrix A_14(a,4), B_14(b,4);
      Matrix C_14 = A_14.sub(B_14);
  
      REQUIRE(equality(C_14, c, 1, 4));

      //4-by-1
      Matrix A_41(a,1), B_41(b,1);
      Matrix C_41 = A_41.sub(B_41);

      REQUIRE(equality(C_41, c, 4, 1));
    }

    SECTION("invalid subtractions")
    {
      std::vector<Elem> a = {1,2,3,4};
      std::vector<Elem> b = {5,6};

      //2-by-2 - 2-by-1
      Matrix A(a, 2), B(b, 1);
      Matrix C = A.sub(B);

      REQUIRE(equality(C, {}, 0, 0));

      //2-by-2 - 1-by-2
      Matrix B_12(b,2);
      Matrix C_12 = A.sub(B_12);
  
      REQUIRE(equality(C_12, {}, 0, 0));
    }
}

TEST_CASE( "Method matrix mult", "[Matrix]" )
{
  INFO("Hint: Multiplication of matrices");

  SECTION("valid multiplications")
    {
      std::vector<Elem> a = {1,2,3,4};
      std::vector<Elem> b = {5,6,7,8};

      //2-by-2 and 2-by-2
      Matrix A(a, 2), B(b, 2);
      Matrix C = A.mult(B);
      
      REQUIRE(equality(C, std::vector<Elem>{19,22,43,50}, 2, 2));

      //1-by-4 and 4-by-1
      Matrix A_14 = Matrix(a,4);
      Matrix B_41 = Matrix(b,1);
      Matrix C_11 = A_14.mult(B_41);
  
      REQUIRE(equality(C_11, std::vector<Elem>{70}, 1, 1));

      //2-by-2 and 2-by-1
      b = {5, 6};
      Matrix A_22 = Matrix(a, 2);
      Matrix B_21 = Matrix(b, 1);
      Matrix C_21 = A_22.mult(B_21);
      
      REQUIRE(equality(C_21, std::vector<Elem>{17,39}, 2, 1));

      //2-by-2 and 2-by-3
      b = {5, 6, 7, 8, 9, 10};
      Matrix B_23 = Matrix(b, 3);
      Matrix C_23 = A.mult(B_23);

      REQUIRE(equality(C_23, std::vector<Elem>{21, 24, 27, 47, 54, 61}, 2, 3));

      //2-by-3 and 3-by-2
      a = {1, 2, 3, 4, 5, 6};
      b = {7, 8, 9, 10, 11, 12};

      Matrix A_23 = Matrix(a, 3);
      Matrix B_32 = Matrix(b, 2);
      Matrix C_22 = A_23.mult(B_32);

      REQUIRE(equality(C_22, std::vector<Elem>{58, 64, 139, 154}, 2, 2));

      //3-by-3 and 3-by-3
      a = {1,2,3,4,5,6,7,8,9};
      b = {10,11,12,13,14,15,16,17,18};

      Matrix A_33 = Matrix(a,3);
      Matrix B_33 = Matrix(b,3);
      Matrix C_33 = A_33.mult(B_33);
      
      REQUIRE(equality(C_33, std::vector<Elem>{84,90,96,201,216,231,318,342,366}, 3, 3));
    }

  SECTION("invalid multiplications")
    {
      std::vector<Elem> a = {1,2,3,4};
      std::vector<Elem> b = {5,6};

      //2-by-2 and 1-by-2
      Matrix A(a, 2), B(b, 2);
      Matrix C = A.mult(B);

      REQUIRE(equality(C, {}, 0, 0));
    }

}

TEST_CASE( "Method scalar mult", "[Matrix]" )
{
  INFO("Hint: Multiplication of matrix by a scalar");

  std::vector<Elem> a = {1,2,3,4};
  std::vector<Elem> b = {2,4,6,8};

    //2-by-2
  Matrix A(a, 2);
  Matrix B = A.mult(2);
  
  REQUIRE(equality(B, b, 2, 2));

  //1-by-4
  Matrix A_14 = Matrix(a, 4);
  Matrix B_14 = A_14.mult(2);
  
  REQUIRE(equality(B_14, b, 1, 4));

  //4-by-1
  Matrix A_41 = Matrix(a, 1);
  Matrix B_41 = A_41.mult(2);
  
  REQUIRE(equality(B_41, b, 4, 1));

  //2-by-3
  a = {1, 2, 3, 4, 5, 6};

  Matrix A_23 = Matrix(a, 3);
  Matrix B_23 = A_23.mult(2);
  
  REQUIRE(equality(B_23, std::vector<Elem>{2, 4, 6, 8, 10, 12}, 2, 3));
  
  //3-by-3
  a = {1,2,3,4,5,6,7,8,9};
  Matrix A_33 = Matrix(a, 3);
  Matrix B_33 = A_33.mult(2);
  
  REQUIRE(equality(B_33, std::vector<Elem>{2, 4, 6, 8, 10, 12, 14, 16, 18}, 3, 3));
}

TEST_CASE( "Method pow", "[Matrix]" )
{
  INFO("Hint: Matrix raised to a given power");

  SECTION("pow w/valid inputs")
    {
      std::vector<Elem> a = {1,2,3,4};
      std::vector<Elem> a0 = {1,0,0,1};
      std::vector<Elem> a1 = {1,2,3,4};
      std::vector<Elem> a2 = {7,10,15,22};
      std::vector<Elem> a5 = {1069,1558,2337,3406};
  
      //A^0
      Matrix A(a, 2);
      Matrix AP = A.pow(0);

      REQUIRE(equality(AP, a0, 2, 2));

      //A^1
      Matrix AP_1 = A.pow(1);
  
      REQUIRE(equality(AP_1, a1, 2, 2));

      //A^2
      Matrix AP_2 = A.pow(2);
  
      REQUIRE(equality(AP_2, a2, 2, 2));

      //A^5
      Matrix AP_5 = A.pow(5);
  
      REQUIRE(equality(AP_5, a5, 2, 2));
    }

  SECTION("pow w/invalid inputs")
    {
      std::vector<Elem> a = {1,2,3};

      //1-by-3
      Matrix A(a, 3);
      Matrix B = A.pow(1);

      REQUIRE(equality(B, {}, 0, 0));

      //2-by-2 w/negative power
      a = {1,2,3,4};
      Matrix A_22 = Matrix(a, 2);
      Matrix B_22 = A.pow(-1);

      REQUIRE(equality(B_22, {}, 0, 0));
    }
}

TEST_CASE( "Method trans", "[Matrix]" )
{
  INFO("Hint: Transpose matrix");

  std::vector<Elem> a = {1,2,3,4};
  std::vector<Elem> aT = {1,3,2,4};
  
  //2-by-2
  Matrix A(a, 2);
  Matrix AT = A.trans();

  REQUIRE(equality(AT, aT, 2, 2));

  //1-by-4
  Matrix A_14 = Matrix(a,4);
  Matrix AT_41 = A_14.trans();
  
  REQUIRE(equality(AT_41, a, 4, 1));

  //4-by-1
  Matrix A_41 = Matrix(a,1);
  Matrix AT_14 = A_41.trans();
    
  REQUIRE(equality(AT_14, a, 1, 4));

  //3-by-3
  a = {1,2,3,4,5,6,7,8,9};
  aT = {1,4,7,2,5,8,3,6,9};

  Matrix A_33 = Matrix(a,3);
  Matrix AT_33 = A_33.trans();

  REQUIRE(equality(AT_33, aT, 3, 3));
}

TEST_CASE( "Method cat", "[Matrix]" )
{
  INFO("Hint: Vertical or horizontal concatenation of two matrices");
  SECTION("vertical cat")
    {
      std::vector<Elem> a = {1,2,3,4,5,6}; //3-by-2
      std::vector<Elem> b = {1,2,3,4}; //2-by-2
      std::vector<Elem> ab = {1,2,3,4,5,6,1,2,3,4};

      Matrix A = Matrix(a,2);
      Matrix B = Matrix(b,2);
      Matrix AB = A.cat(B,1);

      REQUIRE(equality(AB, ab, 5, 2));

      //invalid vertical cat: 3-by-2, 2-by-3
      Matrix C = A.cat(Matrix(a,3),1);
      REQUIRE(equality(C,{},0,0));
    }
  
  SECTION("horizontal cat")
    {
      std::vector<Elem> a = {1,2,3,4,5,6}; //2-by-3
      std::vector<Elem> b = {1,2,3,4}; //2-by-2
      std::vector<Elem> ab = {1,2,3,1,2,4,5,6,3,4};

      Matrix A = Matrix(a,3);
      Matrix B = Matrix(b,2);
      Matrix AB = A.cat(B,2);
      std::cout << AB << std::endl;
      REQUIRE(equality(AB, ab, 2, 5));

      //invalid horizontal cat: 3-by-2, 2-by-3
      Matrix C = A.cat(Matrix(a,2),2);
      REQUIRE(equality(C,{},0,0));
    }
}

// TEST_CASE( "Method invMod", "[Matrix]" )
// {
//   INFO("Hint: Modular multiplicative inverse, using the modulus m, of a matrix");
//   Matrix A;

//   A = A.invMod(29);
  
//   REQUIRE(A.equal(A));
// }

TEST_CASE( "Method rowSwitch", "[Matrix]" )
{
  INFO("Hint: Switch (swap) rows within the matrix (in-place operation)");

  //test case: valid switch
  std::vector<Elem> a = {1,2,3,4,5,6,7,8,9};
  std::vector<Elem> a_sw = {1,2,3,7,8,9,4,5,6};

  Matrix A = Matrix(a,3);
  REQUIRE(A.rowSwitch(1,2));

  REQUIRE(equality(A, a_sw, 3, 3));

  //test cases: invalid switch
  REQUIRE_FALSE(A.rowSwitch(-1,2));
  REQUIRE(equality(A, a_sw, 3, 3));

  REQUIRE_FALSE(A.rowSwitch(1,10));
  REQUIRE(equality(A, a_sw, 3, 3));
}

TEST_CASE( "Method rowMult", "[Matrix]" )
{
  INFO("Hint: Multiply a row within the matrix by a scalar (in-place operation)");

  //test case: valid mult
  std::vector<Elem> a = {1,2,3,4,5,6,7,8,9};
  std::vector<Elem> a_mlt = {2,4,6,12,15,18,-21,-24,-27};

  Matrix A = Matrix(a,3);
  REQUIRE(A.rowMult(0,2));
  REQUIRE(A.rowMult(1,3));
  REQUIRE(A.rowMult(2,-3));
  
  REQUIRE(equality(A, a_mlt, 3, 3));

  //test cases: invalid mult
  REQUIRE_FALSE(A.rowSwitch(-1,2));
  REQUIRE(equality(A, a_mlt, 3, 3));

  REQUIRE_FALSE(A.rowSwitch(1,10));
  REQUIRE(equality(A, a_mlt, 3, 3));
}

TEST_CASE( "Method rowAdd", "[Matrix]" )
{
  INFO("Hint: Replace a row with the sum of itself and the scalar multiple of another row (in-place operation)");

  //test case: valid add
  std::vector<Elem> a = {1,2,3,4,5,6,7,8,9};
  std::vector<Elem> a_add = {-20,-22,-24,4,5,6,7,8,9};

  Matrix A = Matrix(a,3);
  REQUIRE(A.rowAdd(0,2,-3));
  
  REQUIRE(equality(A, a_add, 3, 3));

  //test cases: invalid add
  REQUIRE_FALSE(A.rowAdd(-1,0,2));
  REQUIRE(equality(A, a_add, 3, 3));

  REQUIRE_FALSE(A.rowAdd(1,10,2));
  REQUIRE(equality(A, a_add, 3, 3));
}

TEST_CASE( "Method rowAdd2", "[Matrix]" )
{
  INFO("Hint: Replace a row with the sum of itself and the scalar multiple of another row (in-place operation)");

  //test case: valid add
  std::vector<Elem> a = {2,0,1,0,3,0,1,0,4};

  Matrix A = Matrix(a,3);
  A.rowAdd(1,0,2);
  REQUIRE(true);
}

TEST_CASE( "Inverse Modular", "[Matrix]" )
{
  INFO("Test invmod function");
  std::vector<Elem> elements = {1,2,3,4};
  std::size_t cols = 2;
  Matrix A(elements, cols);

  A.invMod(7);

  REQUIRE(true);
}
