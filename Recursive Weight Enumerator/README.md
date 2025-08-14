# Polar Code's Recursive Weight Enumerator (WEF)

This repository provides MATLAB scripts to **compute the Weight Enumerator Function (WEF)** of Polar Codes using **recursive enumeration**.

These scripts have been implemented based on the recursive process demonstrated in the following paper and used to generate the results for the examples provided in this paper:

M. Rowshan, V-F Dragoi, ``A Tutorial on Weight Structure of Polar Codes," submitted to the Special Issue on Error-Correcting Codes of IEEE BITS magazine.

The WEF gives the number of codewords for each possible Hamming weight and is useful for understanding **minimum distance** and **error performance** of a code.

## Files

### `polar_code_recursive_wef.m`
Main script that:
- Constructs the polar code rate profile using **DEGA** (Density Evolution Gaussian Approximation).
- Identifies frozen and information bits.
- Enumerates all possible **cosets** up to the last frozen bit.
- Calls `CalcA.m` recursively to compute **partial WEFs**.
- Accumulates the results to get the **final WEF**.
- Saves intermediate and final results to:
  - `polar_wef_results.mat`
  - `polar_wef_results.txt`

### `CalcA.m`
Recursive function that:
- Implements the WEF computation based on the recursive structure of polar codes.
- Computes the polynomial coefficients where the index corresponds to the **Hamming weight**.
- Supports both **base case** (`N=1`) and recursive splitting of the code.

---

## Requirements
- MATLAB (tested on R2021a and later)
- Symbolic Math Toolbox (for `syms` and `poly2sym`)

---

## How to Run

1. Clone or download this repository.
2. Open MATLAB and set the working directory to the location of these files.
3. Run:
   ```matlab
   polar_code_recursive_wef
