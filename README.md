# Weight Distribution of Polar Codes

This repository contains two algorithms (MATLAB scripts) for the enumeration of polar codes. These scripts are supplementary resources for the following paper:

M. Rowshan, V-F Dragoi, ``A Tutorial on Weight Structure of Polar Codes," submitted to the Special Issue on Error-Correcting Codes of IEEE BITS magazine.

They are organised into two folders.
1. Closed-form Formulas
This folder contains a script that uses closed-form formulas to compute the partial weight distribution of polar codes for weights less than twice the minimum distance $(<2 \cdot w_{\min})$ using group action on the maximum-degree monomials. This algorithm is fast and can be used for medium and long codes.
2. Recursive Weight Enumerator (WEF)
This folder provides MATLAB scripts to **compute the Weight Enumerator Function (WEF)** of Polar Codes using **recursive enumeration** and returns the complete weight distribution. The use of this script for code longer than 128 is not recommended due to excessive computational complexity.
