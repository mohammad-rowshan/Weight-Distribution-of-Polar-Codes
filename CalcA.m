function Ax = CalcA(N, u_i, i)
    % CalcA:
    %   Implements an algorithm for recursively computing the
    %   Weight Enumerator Function (WEF) A_N^{(i)}(u_{0}^{i})(X)
    %   for a given coset defined by prefix u_0^i.
    %
    % Inputs:
    %   N   : Block length (must be a power of 2 for polar codes)
    %   u_i : Prefix binary vector (u_0, u_1, ..., u_{i}) of length i+1
    %   i   : Index defining the coset (position up to which u_i is fixed)
    %
    % Output:
    %   Ax  : Row vector of polynomial coefficients representing the WEF.
    %         Index k corresponds to the number of codewords with weight k.
    %
    % Notes:
    %   - This function recursively decomposes the code into two halves
    %     based on the polar transform structure.
    %   - In the base case N = 1, the WEF is either 1 (weight 0) or X (weight 1).
    %   - Convolution of polynomials represents weight addition from independent halves.

    % --------------------------
    % Base case: N = 1
    % --------------------------
    if N == 1
        if u_i == 0
            % If bit is 0 → only one codeword of weight 0
            Ax = 1; % Polynomial: 1
        else
            % If bit is 1 → only one codeword of weight 1
            Ax = [0, 1]; % Polynomial: X
        end
        return;
    end
    
    % --------------------------
    % Recursive case
    % --------------------------
    % Split vector u_i into even and odd indexed components
    u_even = u_i(1:2:end); % Positions: u_0, u_2, ...
    u_odd  = u_i(2:2:end); % Positions: u_1, u_3, ...

    % --------------------------
    % Case 1: (i-1) is even
    % --------------------------
    % Means current node is in the 'upper' branch of polar recursion
    if mod(i-1, 2) == 0
        % Compute XOR of even and odd parts (u_even ⊕ u_odd)
        u_xor = mod(u_even + u_odd, 2);

        % Recursively compute WEF for both halves (size N/2)
        Ax_c1 = CalcA(N/2, u_xor, floor(length(u_xor)) - 1);
        Ax_c2 = CalcA(N/2, u_odd, floor(length(u_odd)) - 1);

        % Convolution of polynomials corresponds to combining weights
        Ax = conv(Ax_c1, Ax_c2);

    % --------------------------
    % Case 2: (i-1) is odd
    % --------------------------
    % Means current node is in the 'lower' branch of polar recursion
    else
        % Create two cases for u_{i}:
        %   - When u_{i} = 0
        %   - When u_{i} = 1
        % We need both because the WEF sums over both possibilities.

        % Case: u_i = 0
        u_xor_0 = mod(u_even + [u_odd, 0], 2);  % even ⊕ (odd ⊕ 0)
        Ax0_c1  = CalcA(N/2, u_xor_0, floor(length(u_xor_0)) - 1);
        Ax0_c2  = CalcA(N/2, [u_odd, 0], floor(length([u_odd, 0])) - 1);
        Ax0     = conv(Ax0_c1, Ax0_c2);

        % Case: u_i = 1
        u_xor_1 = mod(u_even + [u_odd, 1], 2);  % even ⊕ (odd ⊕ 1)
        Ax1_c1  = CalcA(N/2, u_xor_1, floor(length(u_xor_1)) - 1);
        Ax1_c2  = CalcA(N/2, [u_odd, 1], floor(length([u_odd, 1])) - 1);
        Ax1     = conv(Ax1_c1, Ax1_c2);

        % Add the two polynomials for the final WEF
        Ax = polyadd(Ax0, Ax1);
    end

    % Uncomment for debugging:
    % disp(sprintf('A_%d(u_0^%d=%s)(X):', N, i, mat2str(u_i)));
    % syms X
    % polynomial_expression = poly2sym(fliplr(Ax), X);
    % disp(polynomial_expression);
end

% ------------------------------------------------
% Helper function to add two polynomials
% ------------------------------------------------
function sum_poly = polyadd(poly1, poly2)
    % Align lengths by padding with zeros
    max_len = max(length(poly1), length(poly2));
    poly1 = [poly1, zeros(1, max_len - length(poly1))];
    poly2 = [poly2, zeros(1, max_len - length(poly2))];

    % Add element-wise
    sum_poly = poly1 + poly2;
end
