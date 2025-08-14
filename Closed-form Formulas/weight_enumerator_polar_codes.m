%% Weight Enumeration of Polar Codes: All weights less than 2w_min
% Copyright (c) 2025, Mohammad Rowshan and Vlad-Florin Dragoi
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification, 
% are permitted provided that the source code retains the above copyright notice.
%
% DESCRIPTION:
% This script computes the weight distribution of polar codes for weights
% less than twice the minimum distance (2*w_min) using group action theory.
% It characterizes codewords based on Type I and Type II structures and
% provides closed-form enumeration formulas.
%
% INPUTS:
%   n: log2(N) where N is the code length
%   R: Code rate (only used calculating code dimention K for rate_profile function: polar codes)
%   odrer: The order of RM code (only used for RM_profile function: RM codes)
% OUTPUTS:
%   r: Maximum degree of monomials (minimum distance = 2^(n-r))
%   w: Vector of weights [w_min, 1.5*w_min, 1.75*w_min, ...]
%   A_w: Multiplicities of codewords with weights w
%
% REFERENCES:
%   [1] V.-F. Dragoi, M. Rowshan, J. Yuan, "On the closed-form weight 
%       enumeration of polar codes: 1.5d-weight codewords," IEEE Trans. 
%       Commun., 2024.
%   [2] M. Rowshan, V.-F. Dragoi, J. Yuan, "Weight structure of low/high-rate 
%       polar codes and its applications," ISIT 2024.
%   [3] V.-F. Dragoi, M. Rowshan, "On Weight Enumeration and Structure
%       Characterization of Polar Codes via Group Actions," 2025, 
%       http://arxiv.org/abs/2504.19544v2 (presented at ISIT'25)


%% Clear workspace and initialize parameters
clear; % Clear workspace
% clc; % Clear command window (commented out for debugging)

%% System parameters
n = 6;              % log2(N) where N is the code length
R = 0.5;            % Code rate of the polar code
order = 3;          % The order of RM code (unused for polar codes). 
design_snr_db = 3;  % design SNR for finding the information set of the polar code (N,K)
N = 2^n;            % Code length
K = N * R;          % Code dimention: the number of information bits
% Uncomment the following lines (and line 56) if you want to enumerate the RM code (order,n).  
% K = 0;
% for r = [0:order]
%     K = K + nchoosek(n,r);
% end

fprintf('Code parameters: (%d,%d) with n=%d\n', N, K, n);


%% Finding the information set I and calling the weight enumerator
I = rate_profile(design_snr_db, N, K)'; % the indices of K most relaible bit-channels
%I = RM_profile(n,order)'; 
% The information set used for the illustrative example in the tutorial paper: 
%I = [23, 26, 27, 28, 29, 30, 31, 37, 38, 39, 41, 42, 43, 44, 45, 46, 47, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63];


[r,w,A_w,A_typ] = weight_enum(sort(I),n);
% Note that Type-I codewords contributes to the weights larger than
% 1.5w_min becuase of mu>=3 condition.

%% Display results
fprintf('Maximum monomial degree: r = %d\n', r);
fprintf('Minimum distance: d_min = 2^(%d-%d) = %d\n', n, r, 2^(n-r));
fprintf('Weight distribution for w < 2*w_min:\n');
for mu = 1:length(w)
    weight_val = w(mu);
    multiplicity = A_w(mu);
    fprintf('  w_%d = %.2f*w_min = %d: %d codewords\n', mu, w(mu)/w(1), weight_val, multiplicity);
end
fprintf('Type-I : %d codewords\n', A_typ(1));
fprintf('Type-II: %d codewords\n', A_typ(2));


%% Weight enumeration for polar codes for weights less than 2*w_min using group action theory.
function [r,w,A_w,A_typ] = weight_enum(I,n)
    % INPUTS:
    %   I: Sorted vector of information bit indices
    %   n: Log2 of code length (code length = 2^n)
    %
    % OUTPUTS:
    %   r: Maximum degree of monomials in I
    %   w: Vector of weights [w_min, 1.5*w_min, 1.75*w_min, 1.875*w_min]
    %   A_w: Multiplicities of codewords with weights w
    
    %% Input validation
    if isempty(I) || ~isvector(I) || any(I ~= floor(I)) || any(I < 0)
        error('I must be a non-empty vector of non-negative integers');
    end
    
    if ~isscalar(n) || n ~= floor(n) || n <= 0
        error('n must be a positive integer');
    end
    
    if max(I) >= 2^n
        warning('Some indices in I exceed code length 2^%d', n);
    end

    %% Find maximum degree of monomials in information set
    r = max(sum(~(dec2bin(I)-'0'),2));

    % Get monomials of maximum degree in information set
    Ir = I(find(sum(~(dec2bin(I)-'0'),2)==r)); % The ascending order in I is needed for conditions in computing alpha

    % Get all monomials of maximum degree (for finding those not in I for Type-I: subtype B1)
    All = [0:2^n-1];
    Ir_all = All(find(sum(~(dec2bin(All)-'0'),2)==r)); % All monomials of degree r
    Ir_not = setdiff(Ir_all, I); % Monomials of degree r not in I

    %% Initialize output variables
    mu_max = max(floor((n-r+2)/2), min(n-r, r));  % Limit mu_max based on constraints
    w = zeros(1, mu_max);
    A_w = zeros(1, mu_max);
    
    %% Pre-compute weight values
    for mu = 1:mu_max
        w(mu) = 2^(n-r) * (2 - 1/2^(mu-1));
    end
    
    %% Main computation loop
    % Track counts for Type I (suptypes A1 and A2) and Type II
    A_typ = [0, 0];  
    A1 = []; A2 = [];
    
    % Store monomial indices for recursive processing
    fi = zeros(n-r, r);  
   
    for i = Ir
        % Process Type II codewords (minimum weight contributions)
        mu = 1;
        f = find(~(reverse(dec2bin(i,n))-'0'));
        fi(1,:) = f;

        A_w(mu) = A_w(mu) + 2^(r+lambda(f,f)); % Counting w_min codewords

        Ir_sub = Ir(Ir>i);
        mu = mu + 1;

        % Check pairs for higher weight codewords
        for j = [i Ir_sub]
            g = find(~(reverse(dec2bin(j,n))-'0'));
            h = intersect(f,g);
            foh = setdiff(f,h);  % foh: f over h : f\h 
            goh = setdiff(g,h); 

            %% Type II codewords: length(h) == r-mu
            if length(h) == r-mu 
                % Compute alpha adjustment for overlapping terms
                foh(1,:) = setdiff(f,h); foh(2,:) = setdiff(g,h); fi(2,:) = g; 
                alpha_sum = alpha(foh,mu); %alpha(foh(1,:),foh(2,:)); %(foh(2)>goh(2) & goh(2)>foh(1)) + (goh(1)>foh(1));
                lambda_sum = r-mu + lambda(h,h) + length(foh(1,:))+lambda(f,foh(1,:)) + length(foh(2,:))+lambda(g,foh(2,:));
                
                A_fg = 2^(lambda_sum - alpha_sum);
                A_w(mu) = A_w(mu) + A_fg;  % 1.5*w_min codewords with mu=2
                A_typ(2) = A_typ(2) + A_fg; 
                
                % Recursive processing for higher mu
                if 2*(mu+1) <= n-r+2
                    A_w1 = A_w(mu+1);
                    [A_w, w,A2] = weights_typeII(h,fi,foh,j,Ir_sub,mu,A_w,w,lambda_sum,r,n,A2);
                    A_typ(2) = A_typ(2) + A_w(mu+1)-A_w1; % subtracting accumulated count (A_w(mu+1)) from previous A_w1
                end
            end

            %% Type I codewords (subtypes A1 and A2)
            for m = [mu:n-r] % mu is changing to consider the number of independent linear forms required for the pair of monomials
                if (m >= 3) && (r >= m) % Type I conditions
                    A_w0 = A_w(m);
                    [A_w, w] = weights_typeI(h,f,g,goh,m,A_w,w,r,n);
                    A_typ(1) = A_typ(1) + A_w(m)-A_w0; % subtracting accumulated count (A_w(mu+1)) from previous A_w1
                    if A_w(m)-A_w0 >0
                        A1 = [A1; [f,g,A_w(m)-A_w0]];
                    end
                end
            end
        end
    end

    % Add %% Type I codewords (subtype B1)
    % Loop over f in Ir_not (monomials of degree r not in I)
    for i_not = Ir_not
        %f_idx = find(All == i_not);
        f = find(~(reverse(dec2bin(i_not,n))-'0'));
        
        % For each possible mu (starting from 3 for Type I)
        for mu = [3:min(r, n-r)]  
            % Generate all possible factorizations f = hh* with deg(h) = r-mu
            if r-mu > 0
                h_candidates = nchoosek(f, r-mu);  % All possible h subsets of f with size r-mu
            else
                % Special case when r-mu = 0 (h is empty)
                h_candidates = [];
            end
            
            % Handle special case when r-mu = 0
            if r == mu
                h_candidates = [h_candidates; []];  % Add empty set as a valid h
            end
            
            for h_idx = 1:size(h_candidates,1)
                % Get h (could be empty)
                if r-mu > 0
                    h = h_candidates(h_idx,:);
                else
                    h = [];
                end
                
                % Calculate h* = f/h
                if isempty(h)
                    h_star = f;
                    lambda_h = 0;
                else
                    h_star = setdiff(f, h);
                    lambda_h = lambda(h,h);
                end
                
                % Check condition: hs* divides the complement of f
                f_complement = setdiff(1:n, f);
                [c, satisfied] = find_dec_compl(f_complement, h_star);
                
                if satisfied
                    % Implement the subtype B1 counting formula
                    A_fg = 2^(r+mu-1+lambda_h+lambda(f,h_star)) * hj2c_permut(h_star, c);
                    
                    % Add to the appropriate weight counter
                    %if mu <= mu_max
                        A_w(mu) = A_w(mu) + A_fg;
                        A_typ(1) = A_typ(1) + A_fg;  % Add to Type I count
                    %end
                end
            end
        end
    end    
end

%% Computing the contribution of Type II codewords with mu > 2 recursively.
function [A_w, w,A2] = weights_typeII(h,fi,foh,i,Ir_sub,mu,A_w,w,lambda_sum_init,r,n,A2)
    Ir_sub = Ir_sub(Ir_sub>i);
    mu = mu + 1;
    for j = Ir_sub
        foh(mu,:) = [0,0];
        g = find(~(reverse(dec2bin(j,n))-'0'));
        if cnt_h(h,fi,foh,g,mu) == mu
            lambda_sum = lambda_sum_init;
            foh(mu,:) = setdiff(g,h);
            fi(mu,:) = g;
            a = alpha(foh,mu); 
            lambda_sum = lambda_sum + length(foh(mu,:))+lambda(g,foh(mu,:));
            A_fg = 2^(lambda_sum - a);
            A_w(mu) = A_w(mu) + A_fg; 
            A2 = [A2; [fi(1,:),fi(2,:),g,h,A_fg]];
            if mu+1 <= n-r
                [A_w, w] = weights_typeII(h,fi,foh,j,Ir_sub,mu,A_w,w,lambda_sum,r,n,A2);
            end
        end
    end
end

%% Computing the contribution of Type I (subtypes A1 and A2) codewords.
function [A_w, w] = weights_typeI(h,f,g,gd,mu,A_w,w,r,n)
    if length(h) < r
        fgc = setdiff([1:n],union(f,g));
        foh = setdiff(f,h);
        goh = setdiff(g,h);
        H = nchoosek(h,r-mu); 
        for hs_idx = 1:size(H,1) 
            hs = H(hs_idx,:);
            if r-mu == 0 
                hs = [];
            end
            hj = setdiff(h,hs);
            [c,satisfied] = find_dec_compl(fgc,hj);
            foh_dec = find_dec_compl(foh,hj);
            if satisfied
                A_fg =  2^(length(hs)+lambda(hs,hs)+length([foh,hj])+lambda(hs,[foh,hj])+length(goh)+lambda(h,goh))*hj2foh_permut(hj,foh_dec)*hj2c_permut(hj,c);
                A_w(mu) = A_w(mu) + A_fg;
            end
        end
    elseif length(h) == r
        fgc = setdiff([1:n],union(f,g)); %complement of f.g
        foh = setdiff(f,h);
        goh = setdiff(g,h);
        H = nchoosek(h,r-mu); 
        for hs_idx = 1:size(H,1) 
            hs = H(hs_idx,:);
            if r-mu == 0 
                hs = [];
            end
            hj = setdiff(h,hs);
            [c,satisfied] = find_dec_compl(fgc,hj); % checking free variables
            foh_dec = find_dec_compl(foh,hj); % checking the variables of the other 
            if satisfied
                A_fg = 2^(length(hs)+lambda(hs,hs)+length([foh,hj])+lambda(hs,[foh,hj])+length(goh)+lambda(h,goh)-1)*hj2foh_permut(hj,foh_dec)*hj2c_permut(hj,c);
                A_w(mu) = A_w(mu) + A_fg;
            end
        end
    end
end

%%
% find_dec_compl(a, h) finds elements in a that are less than elements in h, checking the decreasing condition.
function [compl,satisfied] = find_dec_compl(a,h)
    h = sort(h); a = sort(a); % not needed as all already sorted
    compl = [];
    compl_h_el = {};
    for i = [1:length(h)]
        compl = union(compl, a(find(a<h(i)))); % the ouptut of compl as a column vector but it works as we need only the length
        compl_h_el = [compl_h_el; [a(find(a<h(i)))]];
    end
    satisfied = length(compl_h_el) >= length(h);
end
% hj2c_permut(h, c) calculates the permutation factor for hj joining c
function multip = hj2c_permut(h,c)
    c_sets = set1_lt_set2(c,h);
    multip = 1;
    for i = [1:length(c_sets)]
        multip = multip * (2^length(c_sets{i})-2^(i-1)); %length(c) >> length(h)
    end
end
% hj2foh_permut(h, di) calculates the number of linear independent forms for hj joining goh 
% (making them different with hj joining foh).
function multip = hj2foh_permut(h,di) 
    p_sets = set1_lt_set2(di,h);
    multip = 1;
    for i = [1:length(p_sets)]
        multip = multip * 2^(1+length(p_sets{i})); % +1 representing the translation
    end
end
% set1_lt_set2(a, h) creates a cell array where each cell contains elements of a that are 
% less than the corresponding element in h.
function subsets = set1_lt_set2(a,h)
    subsets = {}; % a cell
    for i = [1:length(h)]
        subsets = [subsets; {[a(find(a<h(i)))]}];
    end
end
% lambda(f, g) computes the sum of lengths of set differences used in the weight enumeration formulas.
function orbit =  lambda(f,g)
    orbit = 0;
    for i = g
        orbit = orbit + length(setdiff(setdiff([1:i-1],g),f));
    end
end
% alpha(foh, mu) computes the alpha adjustment needed for Type II codeword enumeration.
function a = alpha(foh, mu)
    C = nchoosek(1:mu, 2);
    a = 0;
    for m = 1:size(C,1)
        c = C(m,:);
        a = a + (foh(c(1),2) > foh(c(2),2) & foh(c(2),2) > foh(c(1),1)) + ...
                (foh(c(2),1) > foh(c(1),1));
    end
end
% cnt_h(h, fi, foh, g, mu) counts the number of times h intersects with previous fi values.
function cnt = cnt_h(h, fi, foh, g, mu)
    cnt = 1;
    for m = 1:mu-1
        hj = intersect(g, fi(m,:));
        if isequal(hj, h)
            cnt = cnt + 1;
        end
    end
end


%%
% retruns the indices of K-most relaible bit-channels, set I.
% The measure for bit-channels' reliability is Mean-LLRs obtained from Density Evolution 
% by Gaussian Approximation (DEGA) by Chung et. al. and Trifonov, Algorithm by Vangala et. al.
function I = rate_profile(design_snr_db, N, K)
    mllr = zeros(N,1);
    sigma_sq = 1/(2*K/N*power(10,design_snr_db/10));
    mllr(1) = 2/sigma_sq;
    
    for level = 1:log2(N)
        B = 2^level;
        for j = 1:B/2
            T = mllr(j);
            mllr(j) = calc_phi_inv(T);
            mllr(B/2 + j) = 2 * T;
        end
    end
    
    mask = zeros(N,3);
    for i = 0:N-1
        nat(i+1) = bitreversed(i, uint8(log2(N)));
    end
    
    for i = 1:N
        mask(i,:) = [nat(i), mllr(i), 1];
    end
    
    % sort sub-channels by mllr
    mask = sortrows(mask, 2); % direction: ascend (default)
    
    % set info bits to 1 for sub-channels with K largest mllr values
    for i = 1:N-K
        mask(i,3) = 0;
    end
    
    % sort channels with respect to index (in bitreversal order)
    mask = sortrows(mask, 1); % direction: ascend (default)
    I = find(mask(:,3)==1)-1;
end
% bitreversed(i, n) computes the bit-reversed value of i % with n bits.
function dec = bitreversed(i, n)
    dec = bin2dec(fliplr(dec2bin(i, n)));
end

% returns Phi inverse based on piece-wise linear approximation, by Trifonov
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

% returns RM profile 
function I =  RM_profile(m,r)
    N=[0:2^m-1];
    I = [];
     for order = 0: r
        I = [I N(find(sum(~(dec2bin(N)-'0'),2)==order))];
     end
     I = sort(I)';
end
