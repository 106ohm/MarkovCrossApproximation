function U = aca_nd(n, Afiber, tol)
%
% n = [ n(1), ..., n(d) ]
% Aelem(i) -> element A(i(1), ..., i(d));
% Afiber(if, i) -> A(i(1), ..., i(if-1), :, i(if+1), ..., i(d))
%
% U = { U1, ..., Uk }
%
% A = U1(:,1) x ... x Ud(:,1) + ... + U1(:,k) x ... x Ud(:,k)

strategy = 'random';

d = length(n);
U = cell(1, d);

for j = 1 : d
    U{j} = zeros(n(j), 0);
end

i = ones(1, d);
i = find_pivot(Afiber, i, n, strategy);

maxit = 100;

for s = 1 : maxit
    
    for j = 1 : d
        fibers{j} = ( Afiber(j, i) - aca_eval_fiber(U, j, i) );
    end
    
    pivot = fibers{1}(i(1));
    
    fprintf('Step %d, |pivot| = %e\n', s, abs(pivot));
    
    if abs(pivot) <= tol
        return;
    end
    
    % Select the next pivot. 
    U{1} = [ U{1} , fibers{1} ];
    for j = 2 : d
        U{j} = [ U{j}, fibers{j} / pivot ];
    end
    
    i = find_pivot(@(j,i) Afiber(j, i) - aca_eval_fiber(U, j, i), i, n, strategy);
end

end

