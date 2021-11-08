function [U, V] = aca_2d(m, n, Aelem, Arow, Acol, tol)
%ACA_2D 

rowindex = 1;
row = Arow(rowindex); row = row(:);

[~, colindex] = max(abs(row));

U = zeros(m, 0);
V = zeros(m, 0);

nsamples = 10;

for j = 1 : 200
    pivot = row(colindex);
    
    fprintf('Pivot (%d, %d): %e\n', rowindex, colindex, pivot);    
    
    if abs(pivot) <= tol
        ii = randi(m, 1, nsamples); jj = randi(n, 1, nsamples);
        
        for k = 1 : nsamples
            xx = Aelem(ii(k), jj(k)) - U(ii(k), :) * V(jj(k), :)';
            if abs(xx) > abs(pivot)
                pivot = xx;
                colindex = jj(k);
                rowindex = ii(k);
            end
        end
        
        if abs(pivot) <= tol
            return;
        else
            fprintf('New pivot (%d, %d): %e\n', rowindex, colindex, pivot);    
        end
    end
    
    col = Acol(colindex) - U * V(colindex, :)';
    
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
    row = Arow(rowindex); 
    row = row(:); 
    row = row - ( U(rowindex, :) * V' )';
    
    [~, colindex] = max(abs(row));
end

end

