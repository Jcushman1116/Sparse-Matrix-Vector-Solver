# Sparse Matrix Data Structures and Matrix-Vector Product Algorithms

## Overview

Implements and evaluates six algorithms for computing sparse matrix-vector products
using four data structures: Compressed Sparse Row (CSR), Modified CSR, Ellpack-Itpack,
and Compressed Diagonal. Each routine exploits matrix sparsity to avoid unnecessary
computations. Algorithms are validated against MATLAB's built-in matrix-vector product
by comparing relative error and runtime across increasing dimensions.

## Algorithms

**Task 1a — COO to CSR** (`COO_to_CSR_Mv_mult`)
Converts an in row order coordinate (COO) matrix to CSR by building a row pointer
array IA, then computes the product using:

$$
z(i) = \sum_{k=IA(i)}^{IA(i+1)-1} AA(k) \cdot v(JA(k))
$$

**Task 1b — COO to Modified CSR** (`COO_to_modified_CSR_Mv_mult`)
Separates diagonal and off-diagonal elements into distinct arrays, builds a row
pointer for off-diagonal entries, concatenates into a modified storage structure,
then computes:

$$
z(i) = D(i) \cdot v(i) + \sum_{k=\text{rowptr}(i)}^{\text{rowptr}(i+1)-1} \text{offD}(k) \cdot v(\text{col}(k))
$$

**Task 1c — COO to Ellpack-Itpack** (`COO_to_ELL_Mv_mult`)
Stores nonzero values in a dense n × Nd matrix COEF with associated column indices
in JCOEF, where Nd is the maximum number of nonzeros in any row. Product computed as:

$$
z(i) = \sum_{j=1}^{N_d} \text{COEF}(i,j) \cdot v(\text{JCOEF}(i,j))
$$

**Task 2 — Unordered Row COO to Relaxed CSR** (`COO_to_relaxed_CSR_Mv_mult`)
Sorts unordered row indices, builds CSR, then adds elbow room padding after each
row to allow future element insertion without full restructuring. Product uses dot
product over each row's active entries.

**Task 3 — Symmetric CSR** (`Task3`)
Stores only the lower triangular portion of a symmetric matrix in CSR form.
During the product, each stored element $a_{ij}$ contributes to both $z(i)$
and $z(j)$, avoiding redundant storage of the upper triangle:

$$
z(i) \mathrel{+}= a_{ij} \cdot v(j), \quad z(j) \mathrel{+}= a_{ij} \cdot v(i), \quad i \neq j
$$

**Task 4 — Compressed Diagonal** (`task4`)
Detects the diagonal offset k, stores the main diagonal and sub/superdiagonals
into an n × 3 matrix Diag with offset array IOFF = [-k, 0, k], then computes:

$$
z(i) = \sum_{j=1}^{3} \text{Diag}(i,j) \cdot v(i + \text{IOFF}(j))
$$

Note: Task 4 produces correct results for k = 1 only. The algorithm breaks for
k > 1 due to a known indexing issue.

## Methodology

Each algorithm is validated against MATLAB's built-in matrix-vector product using
relative error:

$$
\text{Relative error} = \frac{\|x_{\text{comp}} - x_{\text{true}}\|}{\|x_{\text{true}}\|}
$$

Timing comparisons are run at n = 100, 500, and 1000. All working algorithms
produce relative error on the order of $10^{-16}$ to $10^{-17}$, consistent with
floating point machine precision. All custom implementations run slower than
MATLAB's built-in function across all tested dimensions, with the gap becoming
significant at n = 1000 for the symmetric CSR routine.

## Language

MATLAB

## How to Run

1. Ensure all `.m` files are in the same directory
2. Task 1: run `Task1_tester` — tests all three Task 1 routines against two
   predefined matrices and vectors
3. Tasks 2, 3, and 4: each has its own dedicated testing file — run the
   corresponding tester for each task
4. Each tester prints the output vector, relative error, and timing comparison
   against MATLAB's built-in function
