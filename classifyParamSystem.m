function [singular_vals, inf_vals, inc_vals] = classifyParamSystem(A_sym, b_sym, a_sym)
% classifyParamSystem  Classify a parameterized linear system A(a)x = b(a)
%
% [singular_vals, inf_vals, inc_vals] = classifyParamSystem(A_sym,b_sym,a_sym)
%
% Inputs:
%   A_sym  : n×n symbolic matrix, entries possibly in terms of a_sym
%   b_sym  : n×1 symbolic vector, entries possibly in terms of a_sym
%   a_sym  : symbolic variable (e.g. a_sym = sym('a'))
%
% Outputs:
%   singular_vals : real values of a_sym where det(A_sym)==0
%   inf_vals      : subset of singular_vals where infinitely many solutions
%   inc_vals      : subset of singular_vals where the system is inconsistent
%
% Usage example:
%   syms a
%   A = [2,2*a,3,2; 1,a-1,1-a,a+3; 3,-1,4-a,a+2; -2,2-2*a,-3,a-1];
%   b = [1;0;1;-1];
%   [S,I,In] = classifyParamSystem(A,b,a);

    %% 1) compute the determinant
    detA = simplify( det(A_sym) );

    % --- special case: determinant is identically zero ----------------
    if isequal(detA, sym(0))
        % matrix is singular for every a
        singular_vals = []; 
        inf_vals     = []; 
        inc_vals     = [];
        fprintf('det(A) = 0: matrix is singular for all values of %s.\n', char(a_sym));
        return
    end
    % ------------------------------------------------------------------

    %% 2) find all isolated a where det(A_sym) == 0
    p         = sym2poly(detA);             % polynomial coefficients
    rts_all   = roots(double(p));           % numeric roots
    real_rts  = real(rts_all(abs(imag(rts_all))<1e-8));
    singular_vals = unique(real_rts);

    %% 3) classify each singular a
    n         = size(A_sym,1);
    inf_vals  = [];
    inc_vals  = [];

    for ai = singular_vals.'
        % substitute and convert to numeric
        Ai = double( subs(A_sym, a_sym, ai) );
        bi = double( subs(b_sym, a_sym, ai) );

        rA  = rank(Ai);
        rAb = rank([Ai, bi]);

        if rA == rAb && rA < n
            % underdetermined → infinitely many solutions
            inf_vals(end+1) = ai;      %#ok<AGROW>
        elseif rA < rAb
            % inconsistent
            inc_vals(end+1) = ai;      %#ok<AGROW>
        end
    end

    %% 4) display a summary
    fprintf('\nParameter values where det(A)=0: %s\n', mat2str(singular_vals.',6));
    if ~isempty(inf_vals)
        fprintf(' → Infinitely many solutions at a = %s\n', mat2str(inf_vals,6));
    end
    if ~isempty(inc_vals)
        fprintf(' → Inconsistent at a = %s\n',   mat2str(inc_vals,6));
    end
    if isempty(inf_vals) && isempty(inc_vals)
        fprintf(' → For all roots of det(A)=0, the system has no special cases (unexpected).\n');
    end
end
