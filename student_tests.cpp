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

  std::cout << A << std::endl;
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
  Matrix A(2,2); 

  REQUIRE(A.equal(A));
}

TEST_CASE( "Parameterized constructor 2", "[Matrix]" )
{
  INFO("Hint: Use the parameters to set the matrix elements, vector");
  std::vector<Elem> elements;
  
  //get random values for a matrix
  for (int i = 0; i < 4; i++)
    elements.push_back(rnd_p(gen));

  Matrix A(elements,2);

  REQUIRE(A.equal(A));
  
  std::cout << A << std::endl;
}

// TEST_CASE( "Parameterized constructor 3", "[Matrix]" )
// {
//   INFO("Hint: Use the parameters to set the matrix elements, array");
//   Elem elements[4] = {1, 0, 0, 1};

//   Matrix A(elements,2,2);

//   REQUIRE(A.equal(A));
// }

TEST_CASE( "Element get", "[Matrix]" )
{
  INFO("Hint: Returns the element at specified row, column index");
  Matrix A;

  REQUIRE(A.e(0,0) == 0);
}

TEST_CASE( "Element set", "[Matrix]" )
{
  INFO("Hint: Sets the element at specified row, column index to given value");
  Matrix A;

  REQUIRE(A.e(0,0,1) == false);
}

TEST_CASE( "Method size", "[Matrix]" )
{
  INFO("Hint: Sets the element at specified row, column index to given value");
  Matrix A;

  REQUIRE(A.size(1) == 2);
}

TEST_CASE( "Method equal", "[Matrix]" )
{
  INFO("Hint: Returns true if two matrices are the same, false otherwise");
  Matrix A, B;

  REQUIRE(A.equal(B));
}

TEST_CASE( "Method add", "[Matrix]" )
{
  INFO("Hint: Addition of matrices");
  Matrix A, B, C;

  A = B.add(C);
  
  REQUIRE(A.equal(A));
}

TEST_CASE( "Method sub", "[Matrix]" )
{
  INFO("Hint: Subtraction of matrices");
  Matrix A, B, C;

  A = B.sub(C);
  
  REQUIRE(A.equal(A));
}

TEST_CASE( "Method matrix mult", "[Matrix]" )
{
  INFO("Hint: Multiplication of matrices");
  Matrix A, B, C;

  A = B.mult(C);
  
  REQUIRE(A.equal(A));
}

TEST_CASE( "Method scalar mult", "[Matrix]" )
{
  INFO("Hint: Multiplication of matrix by a scalar");
  Matrix A;
  Scalar k = 1;

  A = A.mult(k);
  
  REQUIRE(A.equal(A));
}

TEST_CASE( "Method pow", "[Matrix]" )
{
  INFO("Hint: Matrix raised to a given power");
  Matrix A;
  Scalar k = 1;

  A = A.pow(k);
  
  REQUIRE(A.equal(A));
}

TEST_CASE( "Method trans", "[Matrix]" )
{
  INFO("Hint: Transpose matrix");
  Matrix A;

  A = A.trans();
  
  REQUIRE(A.equal(A));
}

TEST_CASE( "Method cat", "[Matrix]" )
{
  INFO("Hint: Vertical or horizontal concatenation of two matrices");
  Matrix A, B, C;

  A = B.cat(C,1);
  
  REQUIRE(A.equal(A));
}

TEST_CASE( "Method invMod", "[Matrix]" )
{
  INFO("Hint: Modular multiplicative inverse, using the modulus m, of a matrix");
  Matrix A;

  A = A.invMod(29);
  
  REQUIRE(A.equal(A));
}

TEST_CASE( "Method rowSwitch", "[Matrix]" )
{
  INFO("Hint: Switch (swap) rows within the matrix (in-place operation)");
  Matrix A;
  
  REQUIRE(A.rowSwitch(1,2) == false);
}

TEST_CASE( "Method rowMult", "[Matrix]" )
{
  INFO("Hint: Multiply a row within the matrix by a scalar (in-place operation)");
  Matrix A;
  Scalar k = 1;
  
  REQUIRE(A.rowMult(1,k) == false);
}

TEST_CASE( "Method rowAdd", "[Matrix]" )
{
  INFO("Hint: Replace a row with the sum of itself and the scalar multiple of another row (in-place operation)");
  Matrix A;
  Scalar k = 1;
  
  REQUIRE(A.rowAdd(1,2,k) == false);
}

// My Tests
TEST_CASE("My 2nd Constructor", "[Matrix]"){
  std::vector<Elem> elements = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::size_t cols = 3;

  Matrix A(elements, cols);

  REQUIRE(A.equal(A));
  
  std::cout << A << std::endl;
}

TEST_CASE("My Method equal", "[Matrix]")
{
  INFO("Check if two matrices are equal to each other")
  std::vector<Elem> elements = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::vector<Elem> elements2 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::size_t cols = 3;
  Matrix A(elements, cols);
  Matrix B(elements2, cols);
  
  REQUIRE(A.equal(B));

}
