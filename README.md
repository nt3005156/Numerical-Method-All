# Numerical-Method-All

Newton-Raphson Method Algorithm
ALGORITHM: Newton-Raphson Method for Finding Roots

INPUT: 
    - Function f(x)
    - Derivative f'(x)
    - Initial guess x₀
    - Tolerance ε (epsilon)
    - Maximum iterations N

OUTPUT:
    - Approximate root xₙ
    - Number of iterations

Step 1: START

Step 2: Set iteration count i = 0

Step 3: Input initial guess x₀

Step 4: Input tolerance value ε (e.g., 0.00001)

Step 5: Input maximum iterations N (e.g., 100)

Step 6: DO WHILE i < N
    Step 6.1: Calculate f(xᵢ) and f'(xᵢ)
    
    Step 6.2: IF |f'(xᵢ)| < ε THEN
        Print "Derivative is too small. Method fails."
        STOP
    END IF
    
    Step 6.3: Calculate next approximation:
        xᵢ₊₁ = xᵢ - f(xᵢ) / f'(xᵢ)
    
    Step 6.4: Calculate error = |xᵢ₊₁ - xᵢ|
    
    Step 6.5: Print iteration details: i, xᵢ, f(xᵢ), xᵢ₊₁, error
    
    Step 6.6: IF error < ε THEN
        Print "Root found:", xᵢ₊₁
        Print "Number of iterations:", i+1
        STOP
    END IF
    
    Step 6.7: Update xᵢ = xᵢ₊₁
    Step 6.8: i = i + 1

Step 7: END DO

Step 8: IF i = N THEN
    Print "Maximum iterations reached. Method may not converge."

Step 9: STOP

C program to find the root of a given equation using the Newton-Raphson method:

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAX_ITERATIONS 100
#define EPSILON 0.00001

// Function prototypes
double f(double x);      // The function f(x)
double f_prime(double x); // Derivative f'(x)

int main() {
    double x0, x1, error;
    int iteration = 0;
    
    printf("Newton-Raphson Method for Root Finding\n");
    printf("======================================\n");
    
    // Input initial guess
    printf("Enter initial guess: ");
    scanf("%lf", &x0);
    
    printf("\nIteration\t x0\t\t f(x0)\t\t x1\t\t Error\n");
    printf("----------------------------------------------------------------\n");
    
    // Newton-Raphson iteration
    do {
        // Check if derivative is zero
        if (fabs(f_prime(x0)) < EPSILON) {
            printf("Derivative is zero at x = %.6f. Method fails!\n", x0);
            exit(1);
        }
        
        // Apply Newton-Raphson formula
        x1 = x0 - f(x0) / f_prime(x0);
        
        // Calculate error
        error = fabs(x1 - x0);
        
        // Print iteration details
        printf("%d\t\t %.6f\t %.6f\t %.6f\t %.6f\n", 
               iteration + 1, x0, f(x0), x1, error);
        
        // Update for next iteration
        x0 = x1;
        iteration++;
        
        // Check for maximum iterations
        if (iteration > MAX_ITERATIONS) {
            printf("\nMaximum iterations reached. Method may not converge.\n");
            break;
        }
        
    } while (error > EPSILON);
    
    // Print the root
    printf("\n========================================\n");
    printf("Root found at x = %.6f\n", x1);
    printf("Function value at root: f(x) = %.6f\n", f(x1));
    printf("Number of iterations: %d\n", iteration);
    
    return 0;
}

// Example function: x^3 - 2x - 5 = 0
double f(double x) {
    return x*x*x - 2*x - 5;
}

// Derivative: 3x^2 - 2
double f_prime(double x) {
    return 3*x*x - 2;
}

Sample Output:
Newton-Raphson Method for Root Finding
======================================

Enter initial guess: 2

Iteration	 x0		 f(x0)		 x1		 Error
----------------------------------------------------------------
1		 2.000000	 -1.000000	 2.100000	 0.100000
2		 2.100000	 0.061000	 2.094568	 0.005432
3		 2.094568	 0.000202	 2.094551	 0.000017
4		 2.094551	 0.000000	 2.094551	 0.000000

========================================
Root found at x = 2.094551
Function value at root: f(x) = 0.000000
Number of iterations: 4


Advantages of Newton-Raphson Method:
Very fast convergence (quadratic)

Simple to implement

Usually converges in few iterations

Limitations:
Requires derivative calculation

May fail if derivative is zero

May not converge with poor initial guess

Does not guarantee convergence for all functions
