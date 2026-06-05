# Newton-Raphson Method — Root Finding in C

A clean implementation of the **Newton-Raphson Method** for finding roots of equations, written in C.

---

## What is the Newton-Raphson Method?

The Newton-Raphson method is a powerful numerical technique for finding roots (zeros) of a real-valued function. Starting from an initial guess, it iteratively refines the estimate using the formula:

```
x₁ = x₀ - f(x₀) / f'(x₀)
```

It converges **quadratically**, meaning the number of correct digits roughly doubles with each iteration — making it one of the fastest root-finding methods available.

---

## Algorithm

```
INPUT:  f(x), f'(x), initial guess x₀, tolerance ε, max iterations N
OUTPUT: Approximate root xₙ, number of iterations

1. Set i = 0
2. WHILE i < N:
   a. If |f'(xᵢ)| < ε → STOP (derivative too small)
   b. xᵢ₊₁ = xᵢ - f(xᵢ) / f'(xᵢ)
   c. error = |xᵢ₊₁ - xᵢ|
   d. If error < ε → ROOT FOUND
   e. xᵢ = xᵢ₊₁, i = i + 1
3. If i = N → "Max iterations reached"
```

---

## Code

```c
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAX_ITERATIONS 100
#define EPSILON 0.00001

double f(double x)       { return x*x*x - 2*x - 5; }
double f_prime(double x) { return 3*x*x - 2; }

int main() {
    double x0, x1, error;
    int iteration = 0;

    printf("Enter initial guess: ");
    scanf("%lf", &x0);

    printf("\nIteration\t x0\t\t f(x0)\t\t x1\t\t Error\n");
    printf("----------------------------------------------------------------\n");

    do {
        if (fabs(f_prime(x0)) < EPSILON) {
            printf("Derivative is zero. Method fails!\n");
            exit(1);
        }

        x1    = x0 - f(x0) / f_prime(x0);
        error = fabs(x1 - x0);

        printf("%d\t\t %.6f\t %.6f\t %.6f\t %.6f\n",
               iteration + 1, x0, f(x0), x1, error);

        x0 = x1;
        iteration++;

        if (iteration > MAX_ITERATIONS) {
            printf("Maximum iterations reached.\n");
            break;
        }
    } while (error > EPSILON);

    printf("\nRoot found at x = %.6f\n", x1);
    printf("f(x) at root   = %.6f\n", f(x1));
    printf("Iterations     = %d\n", iteration);

    return 0;
}
```

> The example solves **x³ - 2x - 5 = 0**. To use a different equation, just update `f(x)` and `f_prime(x)`.

---

## Sample Output

```
Enter initial guess: 2

Iteration    x0          f(x0)       x1          Error
----------------------------------------------------------------
1            2.000000    -1.000000   2.100000    0.100000
2            2.100000     0.061000   2.094568    0.005432
3            2.094568     0.000202   2.094551    0.000017
4            2.094551     0.000000   2.094551    0.000000

Root found at x = 2.094551
f(x) at root   = 0.000000
Iterations     = 4
```

---

## How to Compile and Run

```bash
gcc newton_raphson.c -o newton_raphson -lm
./newton_raphson
```

---

## Advantages

- Very fast convergence (quadratic)
- Simple to implement
- Usually converges in very few iterations

## Limitations

- Requires the derivative f'(x)
- Fails if f'(x) = 0 at any iteration point
- May not converge with a poor initial guess
- Does not guarantee convergence for all functions

---



========================================================================================
========================================================================================
========================================================================================
# Gauss Elimination Method — Solving Linear Equations in C

A clean implementation of the **Gauss Elimination Method** for solving systems of linear equations, written in C.

---

## What is the Gauss Elimination Method?

Gauss Elimination is a fundamental algorithm in linear algebra used to solve a system of n linear equations with n unknowns. It works in two phases:

1. **Forward Elimination** — transforms the coefficient matrix into an upper triangular form
2. **Back Substitution** — solves for unknowns starting from the last equation upward

For a system like:
```
2x₁ +  x₂ − x₃ =  8
−3x₁ − x₂ + 2x₃ = −11
−2x₁ + x₂ + 2x₃ = −3
```

We represent it as an augmented matrix `[A|b]` and reduce it step by step.

---

## Algorithm

```
INPUT:  n (number of equations), augmented matrix [A|b] of size n × (n+1)
OUTPUT: Solution vector x[0], x[1], ..., x[n-1]

FORWARD ELIMINATION:
For k = 0 to n-2:
    If A[k][k] == 0 → STOP (zero pivot)
    For i = k+1 to n-1:
        factor = A[i][k] / A[k][k]
        For j = k to n:
            A[i][j] = A[i][j] - factor × A[k][j]

BACK SUBSTITUTION:
x[n-1] = A[n-1][n] / A[n-1][n-1]
For i = n-2 down to 0:
    sum = A[i][n] - Σ (A[i][j] × x[j]) for j = i+1 to n-1
    x[i] = sum / A[i][i]
```

---

## Code

```c
#include <stdio.h>
#include <math.h>

#define MAX 10

int main() {
    float a[MAX][MAX+1], x[MAX];
    int n, i, j, k;
    float factor, sum;

    printf("Gauss Elimination Method\n");
    printf("========================\n");
    printf("Enter number of equations: ");
    scanf("%d", &n);

    printf("Enter augmented matrix [A|b] row by row:\n");
    for (i = 0; i < n; i++) {
        printf("Row %d: ", i+1);
        for (j = 0; j <= n; j++)
            scanf("%f", &a[i][j]);
    }

    /* Forward Elimination */
    for (k = 0; k < n-1; k++) {
        if (fabs(a[k][k]) < 1e-9) {
            printf("Zero pivot encountered. Method fails!\n");
            return 1;
        }
        for (i = k+1; i < n; i++) {
            factor = a[i][k] / a[k][k];
            for (j = k; j <= n; j++)
                a[i][j] -= factor * a[k][j];
        }
    }

    /* Back Substitution */
    x[n-1] = a[n-1][n] / a[n-1][n-1];
    for (i = n-2; i >= 0; i--) {
        sum = a[i][n];
        for (j = i+1; j < n; j++)
            sum -= a[i][j] * x[j];
        x[i] = sum / a[i][i];
    }

    /* Output */
    printf("\nSolution:\n");
    for (i = 0; i < n; i++)
        printf("x[%d] = %.4f\n", i+1, x[i]);

    return 0;
}
```

---

## Sample Output

**Input system:**
```
2   1  -1   8
-3  -1   2  -11
-2   1   2  -3
```

**Output:**
```
Solution:
x[1] =  2.0000
x[2] =  3.0000
x[3] = -1.0000
```

---

## How to Compile and Run

```bash
gcc gauss_elimination.c -o gauss_elimination -lm
./gauss_elimination
```

---

## How It Works (Visual)

```
Original Matrix         After Elimination       Back Substitution
[ 2   1  -1 |  8 ]     [ 2   1  -1 |  8  ]     x3 = -1
[-3  -1   2 | -11]  →  [ 0   0.5 0.5| 1  ]  →  x2 =  3
[-2   1   2 | -3 ]     [ 0   0   1  | -1 ]     x1 =  2
```

---

## Advantages

- Produces exact solution (not iterative)
- Works for any n×n non-singular system
- Simple and straightforward to implement
- Foundation for more advanced methods (LU decomposition, etc.)

## Limitations

- Fails when a pivot element is zero (fix: use partial pivoting)
- Numerical instability for large systems without pivoting
- Computationally O(n³) — slow for very large systems

---

## Topics

`numerical-methods` `linear-algebra` `gauss-elimination` `c-programming` `mathematics` `systems-of-equations`



`numerical-methods` `root-finding` `c-programming` `newton-raphson` `mathematics`
