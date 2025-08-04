% Recursive weight enumeration of polar codes
% -------------------------------------------
% This script calculates the weight enumerator polynomial of a polar code
% using recursive enumeration of possible codewords up to the last frozen bit.
% The weight enumerator polynomial counts how many codewords have a given Hamming weight.

% --------------------------
% Code and design parameters
% --------------------------
N = 2^3;            % Block length (must be a power of two for polar codes)
K = 5;             % Number of information bits (code rate = K/N)

syms X;              % Symbol for polynomial representation

Ax_accum = 0;        % Accumulator for the sum of enumerator polynomials

design_snr_db = 0;   % Design SNR in dB (used for DEGA construction)

% --------------------------
% Construct reliability order and rate profile
% --------------------------
[I, orderRev] = construct_dega(design_snr_db, N, K); % Returns K most reliable sub-channel indices

rate_profile = zeros(1, N);   % Initial rate profile (0 = frozen, 1 = information)
rate_profile(I+1) = 1;        % Mark information bit positions

% --------------------------
% Identify range for enumeration
% --------------------------
last_frozen_index = find(rate_profile == 0, 1, 'last');  % Index of last frozen bit
i = last_frozen_index - 1;  % MATLAB is 1-based; subtract 1 for proper indexing

% Information bit positions up to the last frozen index
info_positions = find(rate_profile(1:last_frozen_index) == 1);

% Extract subvector of rate profile up to last frozen index
sub_RP = rate_profile(1:last_frozen_index);

num_info_bits = length(info_positions);  % Number of info bits in this subvector

% --------------------------
% Generate all possible information patterns
% --------------------------
% Each row in all_combos is a different info-bit pattern (0/1 assignment)
all_combos = dec2bin(0:(2^num_info_bits - 1)) - '0';

% --------------------------
% Loop over all possible combinations of info bits
% --------------------------
for combo = 1:size(all_combos, 1)

    % Make a copy of subvector for modification
    u_i = sub_RP;

    % Convert combo index to binary pattern of length num_info_bits
    bin_pattern = bitget(combo, num_info_bits:-1:1);

    % Assign binary pattern to corresponding information positions
    u_i(info_positions) = bin_pattern;
    
    % Calculate partial weight enumerator polynomial for current pattern
    Ax = CalcA(N, u_i, i);

    % Display the polynomial
    disp(sprintf('A_%d(u_0^%d=%s)(X):', N, i, mat2str(u_i)));
    syms X
    polynomial_expression = poly2sym(fliplr(Ax), X);
    disp(polynomial_expression);

    % Accumulate results
    Ax_accum = polyadd(Ax_accum, Ax);
end

% --------------------------
% Final accumulated enumerator polynomial
% --------------------------
syms X
polynomial_expression = poly2sym(fliplr(Ax_accum), X);
disp(polynomial_expression);


% =========================================================================
% Helper functions
% =========================================================================

% Add two polynomials (aligning their lengths)
function sum_poly = polyadd(poly1, poly2)
    % Align polynomials to the same length
    max_len = max(length(poly1), length(poly2));
    poly1 = [poly1, zeros(1, max_len - length(poly1))];
    poly2 = [poly2, zeros(1, max_len - length(poly2))];
    % Perform element-wise addition
    sum_poly = poly1 + poly2;
end

% Construct reliability sequence via Density Evolution Gaussian Approximation (DEGA)
function [I, orderREv] = construct_dega(design_snr_db, N, K)
    mllr = zeros(N,1); % Mean Log-Likelihood Ratios for sub-channels
    sigma_sq = 1/(2*K/N*power(10, design_snr_db/10));
    mllr(1) = 2/sigma_sq;

    % Recursively update MLLRs for all polarization levels
    for level = 1:log2(N)
        B = 2^level;
        for j = 1:B/2
            T = mllr(j);
            mllr(j) = calc_phi_inv(T); % Upper channel update
            mllr(B/2 + j) = 2 * T;     % Lower channel update
        end
    end
    
    % Natural to bit-reversed ordering
    for i = 0:N-1
        nat(i+1) = bitreversed(i, uint8(log2(N)));
    end
    
    % Create mask: [index, reliability, flag]
    mask = zeros(N,3);
    for i = 1:N
        mask(i,:) = [nat(i), mllr(i), 1];
    end
    
    % Sort by reliability
    mask = sortrows(mask, 2);
    
    % Mark N-K least reliable channels as frozen (flag=0)
    for i = 1:N-K
        mask(i,3) = 0;
    end
    
    % Extract reliability ordering
    order = mask(:,1);
    orderREv = order(end:-1:1);
    
    % Sort back to original index order
    mask = sortrows(mask, 1);
    I = find(mask(:,3) == 1) - 1;
end

% Bit-reversal function (same as bitrevorder in DSP toolbox)
function dec = bitreversed(i, n)
    dec = bin2dec(fliplr(dec2bin(i, n)));
end

% Piecewise linear approximation of φ⁻¹ function
function phi_inv = calc_phi_inv(x)
    if (x > 12)
        phi_inv = 0.9861 * x - 2.3152;
    elseif (x <= 12 && x > 3.5)
        phi_inv = x*(0.009005 * x + 0.7694) - 0.9507;
    elseif (x <= 3.5 && x > 1)
        phi_inv = x*(0.062883*x + 0.3678) - 0.1627;
    else
        phi_inv = x*(0.2202*x + 0.06448);
    end
end
