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

## Topics

`numerical-methods` `root-finding` `c-programming` `newton-raphson` `mathematics`
