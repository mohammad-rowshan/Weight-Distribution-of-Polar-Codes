function Ax = CalcA(N, u_i, i)
    % CalcA: Implements Algorithm 1 to compute WEFs A_N^{(i)}(u_{i-1}, 0)(X) and A_N^{(i)}(u_{i-1}, 1)(X)
    % Inputs:
    %   N: Block length (power of 2)
    %   u_i: Prefixed binary vector (u_0, u_1, ..., u_{i}) of length i+1
    %   i: Index defining the coset
    % Outputs:
    %   Ax: Coefficients of WEF A_N^{(i)}(u_{i})(X)
    
    % Base case: N = 1
    if N == 1
        if u_i == 0
            Ax = 1; % Polynomial 1 (weight 0)
        else
            Ax = [0, 1]; % Polynomial X (weight 1)
        end
        return;
    end
    
    % Compute subvectors
    if mod(i-1, 2) == 0
        % Even index: Use u_{i-1, even} ⊕ u_{i-1, odd} and u_{i-1, odd}
        u_even = u_i(1:2:end); % (u_0, u_2, ...)
        u_odd = u_i(2:2:end); % (u_1, u_3, ...)
        u_xor = mod(u_even + u_odd, 2); % u_{i-1, even} ⊕ u_{i-1, odd}
        
        % Recursive calls with N/2
        Ax_c1 = CalcA(N/2, u_xor, floor(length(u_xor))-1);
        Ax_c2 = CalcA(N/2, u_odd, floor(length(u_odd))-1);
        Ax = conv(Ax_c1, Ax_c2);
    else
        u_even = u_i(1:2:end); % (u_0, u_2, ...)
        u_odd = u_i(2:2:end); % (u_1, u_3, ...)
        u_xor_0 = mod(u_even + [u_odd,0], 2); % u_{i, even} ⊕ (u_{i, odd} ⊕ 0)
        u_xor_1 = mod(u_even + [u_odd,1], 2); % u_{i, even} ⊕ (u_{i, odd} ⊕ 1)
        
        % Recursive calls with N/2
        Ax0_c1 = CalcA(N/2, u_xor_0, floor(length(u_xor_0))-1);
        Ax0_c2 = CalcA(N/2, [u_odd,0], floor(length([u_odd,0]))-1);
        Ax0 = conv(Ax0_c1, Ax0_c2);
        Ax1_c1 = CalcA(N/2, u_xor_1, floor(length(u_xor_1))-1);
        Ax1_c2 = CalcA(N/2, [u_odd,1], floor(length([u_odd,1]))-1);
        Ax1 = conv(Ax1_c1, Ax1_c2);

        Ax = polyadd(Ax0, Ax1);
    end
%     disp(sprintf('A_%d(u_0^%d=%s)(X):', N, i, mat2str(u_i)));
%     syms X
%     polynomial_expression = poly2sym(fliplr(Ax), X);
%     disp(polynomial_expression);
end

function sum_poly = polyadd(poly1, poly2)
    % Helper function to add two polynomials
    % Align polynomials to the same length
    max_len = max(length(poly1), length(poly2));
    poly1 = [poly1, zeros(1, max_len - length(poly1))];
    poly2 = [poly2, zeros(1, max_len - length(poly2))];
    sum_poly = poly1 + poly2;
end