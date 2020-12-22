function [U, V] = aca_2d(m, n, Aelem, Arow, Acol)
%ACA_2D 

rowindex = 1;
row = Arow(rowindex);

[~, colindex] = max(abs(row));

U = zeros(m, 0);
V = zeros(m, 0);

for j = 1 : 6
    pivot = row(colindex);
    col = Acol(colindex) - U * V(colindex, :)';
    
    fprintf('Pivot (%d, %d): %e\n', rowindex, colindex, pivot);
    
    U = [ U, col / pivot ];
    V = [ V, row ];
    
    % Select the next pivot
    [~, nrowindex] = max(abs(...
        col([ 1 : rowindex - 1, rowindex + 1 : end ])) ...
    );
    
    if nrowindex >= rowindex
        nrowindex = nrowindex + 1;
    end
    
    rowindex = nrowindex;
    row = Arow(rowindex) - ( U(rowindex, :) * V' )';
    
    [~, colindex] = max(abs(row));
end

end

