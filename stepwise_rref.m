function stepwise_rref(A)
    % Step-by-step RREF transformation of matrix A
    [m, n] = size(A);
    fprintf('Original Matrix:\n');
    disp(A);

    pivotRow = 1;
    
    for col = 1:n
        % Find pivot in the current column
        [~, maxRow] = max(abs(A(pivotRow:m, col)));
        maxRow = maxRow + pivotRow - 1;

        % Skip column if pivot is zero
        if A(maxRow, col) == 0
            continue;
        end

        % Swap rows if necessary
        if maxRow ~= pivotRow
            A([pivotRow, maxRow], :) = A([maxRow, pivotRow], :);
            fprintf('Swapped row %d with row %d:\n', pivotRow, maxRow);
            disp(A);
        end

        % Scale pivot row to make pivot = 1
        pivotValue = A(pivotRow, col);
        A(pivotRow, :) = A(pivotRow, :) / pivotValue;
        fprintf('Divided row %d by %.4f to make pivot 1:\n', pivotRow, pivotValue);
        disp(A);

        % Eliminate all other entries in the pivot column
        for r = 1:m
            if r ~= pivotRow
                factor = A(r, col);
                A(r, :) = A(r, :) - factor * A(pivotRow, :);
                fprintf('Row %d - (%.4f * Row %d):\n', r, factor, pivotRow);
                disp(A);
            end
        end
        
        pivotRow = pivotRow + 1; % Move to next row
        if pivotRow > m
            break;
        end
    end

    fprintf('Final RREF Matrix:\n');
    disp(A);
end
