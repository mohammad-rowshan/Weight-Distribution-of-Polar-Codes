% ========================================================================
% polar_code_recursive_wef.m
% ========================================================================
% Recursive Weight Enumeration of Polar Codes
% Educational Implementation with Result Saving
%
% Author: Mohammad Rowshan
% Date: 2025-08-04
% License: © 2025 Mohammad Rowshan. All rights reserved.
%          Permission to use for educational purposes only.
%
% Description:
%   - Calculates the weight enumerator polynomial (WEF) of a polar code
%     using recursive enumeration of codewords up to the last frozen bit.
%   - Saves each coset’s partial WEF and the final accumulated WEF to file.
%   - Outputs results in .mat and .txt format for later analysis.
%
% ========================================================================
clc;
% --------------------------
% Code and design parameters
% --------------------------
N = 2^3;            % Block length (must be a power of two for polar codes)
K = 4;              % Number of information bits

syms X;              % Symbol for polynomial representation

Ax_accum = 0;        % Accumulator for total WEF

design_snr_db = 0;   % Design SNR in dB for DEGA construction

% --------------------------
% Construct reliability order and rate profile
% --------------------------
[I, orderRev] = construct_dega(design_snr_db, N, K);

rate_profile = zeros(1, N);   
rate_profile(I+1) = 1;        % Mark information bit positions

% --------------------------
% Identify range for enumeration
% --------------------------
last_frozen_index = find(rate_profile == 0, 1, 'last');
i = last_frozen_index - 1;  

info_positions = find(rate_profile(1:last_frozen_index) == 1);
sub_RP = rate_profile(1:last_frozen_index);
num_info_bits = length(info_positions);

% --------------------------
% Generate all information bit patterns
% --------------------------
all_combos = dec2bin(0:(2^num_info_bits - 1)) - '0';

% --------------------------
% Store results
% --------------------------
results = struct;
results.N = N;
results.K = K;
results.design_snr_db = design_snr_db;
results.rate_profile = rate_profile;
results.partial_polynomials = cell(size(all_combos,1), 1);

% --------------------------
% Loop over all combinations
% --------------------------
for combo = 1:size(all_combos, 1)

    % Make copy of subvector
    u_i = sub_RP;

    % Assign binary pattern
    bin_pattern = bitget(combo, num_info_bits:-1:1);
    u_i(info_positions) = bin_pattern;
    
    % Compute partial WEF
    Ax = CalcA(N, u_i, i);

    % Store partial result
    results.partial_polynomials{combo} = Ax;

    % Display
    fprintf('Coset %d / %d: u_i=%s\n', combo, size(all_combos,1), mat2str(u_i));
    disp(poly2sym(fliplr(Ax), X));

    % Accumulate into total WEF
    Ax_accum = polyadd(Ax_accum, Ax);
end

% --------------------------
% Final total WEF
% --------------------------
results.final_polynomial = Ax_accum;
results.final_polynomial_X = poly2sym(fliplr(Ax_accum), X);
disp('Final accumulated WEF:');
disp(poly2sym(fliplr(Ax_accum), X));

% --------------------------
% Save results
% --------------------------
save('polar_wef_results.mat', 'results');

fid = fopen('polar_wef_results.txt', 'w');
fprintf(fid, 'Polar Code WEF Enumeration\n');
fprintf(fid, 'N=%d, K=%d, Design SNR=%.2f dB\n\n', N, K, design_snr_db);

for c = 1:length(results.partial_polynomials)
    fprintf(fid, 'Coset %d: %s -> %s\n', c, ...
        mat2str(all_combos(c,:)), mat2str(results.partial_polynomials{c}));
end

fprintf(fid, '\nFinal WEF: %s\n', mat2str(results.final_polynomial));
fprintf(fid, '\nFinal WEF(X): %s\n', results.final_polynomial_X);
fclose(fid);

disp('Results saved to polar_wef_results.mat and polar_wef_results.txt');



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

