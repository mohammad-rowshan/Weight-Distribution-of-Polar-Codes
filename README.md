# Weight Distribution of Polar Codes

This repository contains two MATLAB scripts for weight enumeration and returning partial or the entire weight distribution (a.k.a. weight spectrum) of polar codes.

They are organised into two folders.
1. Closed-form Formulas:
This folder contains a script that uses closed-form formulas to compute the partial weight distribution of polar codes for weights less than twice the minimum distance $(<2 \cdot w_{\min})$ using group action on the maximum-degree monomials. Since this algorithm is based on closed-form formulas presented in [1-2] and further elaborated in [4], it is fast and can be used for short, medium, and long codes.

3. Recursive Weight Enumerator (WEF):
This folder provides MATLAB scripts to **compute the Weight Enumerator Function (WEF)** of Polar Codes using **recursive enumeration** and returns the complete weight distribution for short codes. This algorithm is a variation of the algorithm provided in [3] and explained in detail with an illustrative example in [4].


[1] V. -F. Drăgoi, M. Rowshan and J. Yuan, "On the Closed-Form Weight Enumeration of Polar Codes: 1.5d -Weight Codewords," in IEEE Transactions on Communications, vol. 72, no. 10, pp. 5972-5987, Oct. 2024, doi: 10.1109/TCOMM.2024.3394749.

[2] M. Rowshan, V. Drăgoi and J. Yuan, "Weight Structure of Low/High-Rate Polar Codes and Its Applications," 2024 IEEE International Symposium on Information Theory (ISIT), Athens, Greece, 2024, pp. 2945-2950, doi: 10.1109/ISIT57864.2024.10619618.

[3] H. Yao, A. Fazeli and A. Vardy, "A Deterministic Algorithm for Computing the Weight Distribution of Polar Code," in IEEE Transactions on Information Theory, vol. 70, no. 5, pp. 3175-3189, May 2024.

[4] M. Rowshan, V-F Dragoi, ``A Tutorial on Weight Structure of Polar Codes," submitted to the Special Issue on Error-Correcting Codes of IEEE BITS magazine.
