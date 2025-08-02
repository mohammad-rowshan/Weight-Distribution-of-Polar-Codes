N = 2^6;
K = N/2;

syms X;

Ax_accum = 0;

design_snr_db = 0;
[I, orderRev] = construct_dega(design_snr_db, N, K); % the indices of K information bits
I(18) = 37; I = sort(I); % modification
rate_profile = zeros(1,N);
rate_profile(I+1) = 1;

% --- Find last frozen position ---
last_frozen_index = find(rate_profile == 0, 1, 'last');
i = last_frozen_index-1; %Note that in MATLAB vectors are indexed from 1.

% --- Find indices of 1's from index 1 to last_frozen_index ---
info_positions = find(rate_profile(1:last_frozen_index) == 1);

sub_RP = rate_profile(1:last_frozen_index);

num_info_bits = length(info_positions);

% --- Generate all possible combinations of {0,1} for these positions ---
all_combos = dec2bin(0:(2^num_info_bits - 1)) - '0'; % Each row is a combination

% --- Loop over all combinations ---
for combo = 1:size(all_combos,1)

    % Make a copy of the subvector
    u_i = sub_RP;

    % Convert combination number to binary vector of length num_info_bits
    bin_pattern = bitget(combo, num_info_bits:-1:1);

    % Assign binary pattern to info positions
    u_i(info_positions) = bin_pattern;
    
    Ax = CalcA(N, u_i, i);

    disp(sprintf('A_%d(u_0^%d=%s)(X):', N, i, mat2str(u_i)));
    syms X
    polynomial_expression = poly2sym(fliplr(Ax), X);
    disp(polynomial_expression);

    Ax_accum = polyadd(Ax_accum , Ax);
end

syms X
polynomial_expression = poly2sym(fliplr(Ax_accum), X);
disp(polynomial_expression);





function sum_poly = polyadd(poly1, poly2)
    % Helper function to add two polynomials
    % Align polynomials to the same length
    max_len = max(length(poly1), length(poly2));
    poly1 = [poly1, zeros(1, max_len - length(poly1))];
    poly2 = [poly2, zeros(1, max_len - length(poly2))];
    sum_poly = poly1 + poly2;
end



% retruns Mean-LLRs (a measure for sub-channels' reliability) obtained from Density Evolution by Gaussian Approximation (DEGA) 
function [I, orderREv] =  construct_dega(design_snr_db, N, K)
    mllr = zeros(N,1);
    sigma_sq = 1/(2*K/N*power(10,design_snr_db/10));
    mllr(1) = 2/sigma_sq;
    for level = 1:log2(N)
        B = 2^level;
        for j = 1:B / 2
            T = mllr(j);
            mllr(j) = calc_phi_inv(T);
            mllr(B / 2 + j) = 2 * T;
        end
    end
    
    mask = zeros(N,3);
    for i = 0:N-1
        nat(i+1) = bitreversed(i,uint8(log2(N)));
    end
    %nat = bitrevorder(0:N-1);
    for i = 1:N
        mask(i,:) = [nat(i), mllr(i), 1];
    end
    % sort sub-channels by mllr
    mask = sortrows(mask,2); %direction: ascend (default)
    % set info bits to 1 for sub-channels with K largest mllr values
    for i = 1:N-K
        mask(i,3) = 0;
    end
    order = mask(:,1);
    orderREv = order(end:-1:1);
    % sort channels with respect to index (in bitreversal order; line 42
    mask = sortrows(mask,1); %direction: ascend (default)
    I = find(mask(:,3)==1)-1;
end

function dec = bitreversed(i,n) % bitrevorder() is in singal processing toolbox.
    dec = bin2dec(fliplr(dec2bin(i,n)));
end

% returns Phi inverse based on piece-wise linear approximation
function phi_inv = calc_phi_inv(x)
    if (x>12)
        phi_inv = 0.9861 * x - 2.3152;
    elseif (x<=12 && x>3.5)
        phi_inv = x*(0.009005 * x + 0.7694) - 0.9507;
    elseif (x<=3.5 && x>1)
        phi_inv = x*(0.062883*x + 0.3678)- 0.1627;
    else
        phi_inv = x*(0.2202*x + 0.06448);
    end
end


% returns the minimum distance of the code and the number of codewords with minumum distance (error coefficient)
function [dmin, A_dmin] = err_coeff(I,N)
    d = min(sum(dec2bin(I)-'0',2));
    dmin = 2^d; n = log2(N); A_dmin = 0;
    B = find(sum(dec2bin(I)-'0',2)==d);
    for i = B'
        Ki_size = n - d;
        for x = find(dec2bin(I(i),n)-'0'==1)
            %if x>1
                ii = dec2bin(bitxor(N-1,I(i)),n)-'0';
                Ki_size = Ki_size + sum(ii(1:x-1));
            %end
        end
        A_dmin = A_dmin + 2^Ki_size;
    end
end
