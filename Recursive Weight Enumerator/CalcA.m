% ========================================================================
% CalcA.m
% ========================================================================
% Recursive computation of polar code Weight Enumerator Function (WEF).
% Educational version with detailed explanatory comments.
%
% Author: Mohammad Rowshan
% Date: 2025-08-04
% License: © 2025 Mohammad Rowshan. All rights reserved.
%          Permission to use for educational purposes only.
%
% Description:
%   - Implements recursive decomposition of the polar code structure
%     to compute the WEF for a given coset defined by prefix bits u_0^i.
%   - Displays intermediate recursion steps for educational purposes.
%   - Returns polynomial coefficients where index k = number of codewords
%     with Hamming weight k.
%
% ========================================================================

function Ax = CalcA(N, u_i, i, depth)
    % Inputs:
    %   N      : Block length (power of 2)
    %   u_i    : Prefix binary vector (u_0, ..., u_{i}) of length i+1
    %   i      : Index defining the coset
    %   depth  : Recursion depth for indentation in educational printing
    %
    % Output:
    %   Ax     : Row vector of polynomial coefficients representing WEF

    if nargin < 4
        depth = 0; % recursion depth for indentation
    end
    indent = repmat(' ', 1, depth*4); % indentation for visualising recursion

    %fprintf('%s[N=%d] u_i=%s, i=%d\n', indent, N, mat2str(u_i), i);

    % --------------------------
    % Base case: N = 1
    % --------------------------
    if N == 1
        if u_i == 0
            Ax = 1;       % Polynomial: 1 (weight 0)
        else
            Ax = [0, 1];  % Polynomial: X (weight 1)
        end
        %fprintf('%s-> Base case polynomial: %s\n', indent, mat2str(Ax));
        return;
    end

    % --------------------------
    % Recursive case
    % --------------------------
    % Split prefix vector into even and odd positions
    u_even = u_i(1:2:end);
    u_odd  = u_i(2:2:end);

    % Case 1: (i-1) is even → upper branch of recursion
    if mod(i-1, 2) == 0
        % XOR even and odd positions
        u_xor = mod(u_even + u_odd, 2);

        % Recursively compute WEFs for both halves
        Ax_c1 = CalcA(N/2, u_xor, length(u_xor)-1, depth+1);
        Ax_c2 = CalcA(N/2, u_odd, length(u_odd)-1, depth+1);

        % Combine via polynomial convolution (weight addition)
        Ax = conv(Ax_c1, Ax_c2);

    % Case 2: (i-1) is odd → lower branch of recursion
    else
        % Two cases for the last bit: u_i = 0 and u_i = 1

        % Case: u_i = 0
        u_xor_0 = mod(u_even + [u_odd,0], 2);
        Ax0_c1  = CalcA(N/2, u_xor_0, length(u_xor_0)-1, depth+1);
        Ax0_c2  = CalcA(N/2, [u_odd,0], length([u_odd,0])-1, depth+1);
        Ax0     = conv(Ax0_c1, Ax0_c2);

        % Case: u_i = 1
        u_xor_1 = mod(u_even + [u_odd,1], 2);
        Ax1_c1  = CalcA(N/2, u_xor_1, length(u_xor_1)-1, depth+1);
        Ax1_c2  = CalcA(N/2, [u_odd,1], length([u_odd,1])-1, depth+1);
        Ax1     = conv(Ax1_c1, Ax1_c2);

        % Combine both cases by polynomial addition
        Ax = polyadd(Ax0, Ax1);
    end

    %fprintf('%s-> Combined polynomial: %s\n', indent, mat2str(Ax));
end

% --------------------------
% Helper: Polynomial addition
% --------------------------
function sum_poly = polyadd(poly1, poly2)
    % Adds two polynomials by aligning lengths
    max_len = max(length(poly1), length(poly2));
    poly1 = [poly1, zeros(1, max_len - length(poly1))];
    poly2 = [poly2, zeros(1, max_len - length(poly2))];
    sum_poly = poly1 + poly2;
end
