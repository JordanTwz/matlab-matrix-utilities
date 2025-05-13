% Compute an orthogonal basis U for the columns of A, exactly.
function U = gram_schmidt(A)
    if ~isa(A,'sym')
        A = sym(A); % lift to symbolic
    end
    [n, m] = size(A);
    U = sym(zeros(n, m));
    for j = 1:m
        u = A(:, j);
        for i = 1:j-1
            coeff = (U(:,i).' * A(:,j)) / (U(:,i).' * U(:,i));
            u = u - coeff * U(:,i);
        end
        U(:, j) = u;
    end
end
