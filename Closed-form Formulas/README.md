# Closed-form Weight Enumeration of Polar Codes

This MATLAB script computes the weight distribution of polar codes for weights less than twice the minimum distance $(<2 \cdot w_{\min})$ using group action on the maximum-degree monomials. It characterises codewords based on Type-I and Type-II structures and provides closed-form enumeration formulas, as described in the referenced works.

## Description
The script calculates the partial weight distribution of polar codes (or Reed-Muller codes, if configured) by analysing the information set $I$. It focuses on weights up to $2 \cdot w_{\min}$, where $w_{\min} = 2^{n-r}$ is the minimum distance, determined by the maximum degree $r$ of monomials in $I$. The script leverages permutation group action to count codewords of specific weights efficiently.

## Features
- Computes the maximum monomial degree $r$ and minimum distance $w_{\text{min}}$.
- Enumerates codeword multiplicities for weights $w_{\text{min}}, 1.5 \cdot w_{\text{min}}, 1.75 \cdot w_{\text{min}}, \ldots$.
- Distinguishes between Type-I and Type-II codeword structures.
- Supports both polar codes (using `rate_profile`) and Reed-Muller (RM) codes (using `RM_profile`).
- Includes input validation and error handling.

## Prerequisites
- MATLAB (tested with R2020a and later).
- No additional toolboxes required.

## Usage
1. **Clone the Repository**:
   ```bash
   git clone <repository-url>
   cd <repository-directory>
   ```

2. **Run the Script**:
   Open MATLAB, navigate to the script directory, and run:
   ```matlab
   weight_enum_script
   ```

3. **Configure Parameters**:
   Modify the system parameters in the script to suit your needs:
   ```matlab
   n = 6;              % log2(N), where N is the code length
   R = 0.5;            % Code rate (for polar codes)
   ...
   design_snr_db = 3;  % Design SNR for polar code reliability
   ...
   I = rate_profile(design_snr_db, N, K)';
   ```
   
4. **Optional: Reed-Muller Codes**:
   To enumerate weights for an RM code, uncomment the lines for `K` calculation and use `RM_profile`:
   ```matlab
   order = 3;          % Order of RM code (for RM codes)
   ...
   K = 0;
   for r = [0:order]
       K = K + nchoosek(n,r);
   end
   ...
   I = RM_profile(n, order)';
   ```

5. **Custom Information Set**:
   To use a specific information set $I$, modify the script to include your indices, e.g.:
   ```matlab
   I = [6, 9, 20, 26, 28, 37, 41, 42, 44, 50, 52]; % Example 
   ```

## Inputs
- `n`: Logarithm base-2 of the code length $(N = 2^n)$.
- `R`: Code rate (used to compute $K = N \cdot R$ for polar codes).
- `order`: Order of the RM code (used for `RM_profile`).
- `design_snr_db`: Design SNR in dB for polar code reliability calculation.
- `I`: Optional user-defined information set (vector of indices).

## Outputs
- `r`: Maximum degree of monomials in $I$, where minimum distance is $2^{n-r}$.
- `w`: Vector of weights $( w_{\text{min}}, 1.5 \cdot w_{\text{min}}, 1.75 \cdot w_{\text{min}}, \ldots$).
- `A_w`: Multiplicities of codewords for each weight in `w`.
- `A_typ`: Counts of Type I and Type II codewords.

## Example
For the illustrative example in ``A Tutorial on Weight Structure of Polar Code," submitted to IEEE BITS:
```matlab
I = [23, 26, 27, 28, 29, 30, 31, 37, 38, 39, 41, 42, 43, 44, 45, 46, 47, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63];
[r, w, A_w, A_typ] = weight_enum(I, n);
```
Output (example):
```
Code parameters: (64,32) with n=6
Maximum monomial degree: r = 3
Minimum distance: d_min = 2^(6-3) = 8
Weight distribution for w < 2*w_min:
  w_1 = 1.00*w_min = 8: 920 codewords
  w_2 = 1.50*w_min = 12: 25472 codewords
  w_3 = 1.75*w_min = 14: 32768 codewords
Type-I : 32768 codewords
Type-II: 25472 codewords
```

## References
1. V. -F. Drăgoi, M. Rowshan and J. Yuan, "On the Closed-Form Weight Enumeration of Polar Codes: 1.5d -Weight Codewords," in IEEE Transactions on Communications, vol. 72, no. 10, pp. 5972-5987, Oct. 2024, doi: 10.1109/TCOMM.2024.3394749.
2. M. Rowshan, V. Drăgoi and J. Yuan, "Weight Structure of Low/High-Rate Polar Codes and Its Applications," 2024 IEEE International Symposium on Information Theory (ISIT), Athens, Greece, 2024, pp. 2945-2950, doi: 10.1109/ISIT57864.2024.10619618.
3. V.-F. Dragoi and M. Rowshan, "On weight enumeration and structure characterization of polar codes via group actions," arXiv:2504.19544v2, 2025. [Online]. Available: https://arxiv.org/abs/2504.19544 (presented at ISIT'25).

## License
Copyright (c) 2025, Mohammad Rowshan and Vlad-Florin Dragoi. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the source code retains the above copyright notice.

## Contact
For questions or feedback, contact the authors: mrowshan at ieee dot org
