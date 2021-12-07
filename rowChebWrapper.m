function row = rowChebWrapper(f, pts)
%

n = 15;
a = min(pts);
b = max(pts);

x = chebpts(n, [a, b]);
fx = x;

for j = 1 : length(x)
    fx(j) = f(x(j));
end

res = inf;

while res > 100 * eps * length(x)
    n = 2 * n - 1;
    
    % Evaluate the function at the new points
    x = chebpts(n, [a, b]);
    nfx = zeros(length(x)-2,1);
    nfx(2:2:end) = fx(2:end-1);
    fx = [fx(1) ; nfx ; fx(end)];
    for j = 2 : 2 : length(x)
        fx(j) = f(x(j));
    end
    
    sz = length(fx) - 1;
    v = [fx(end:-1:1) ; fx(2:sz) ];
    c = real(fft(v)) / (sz);
    c = [c(1) / 2 ; c(2:sz) ; c(sz+1) / 2 ];

    % Estimate the residual
    mp = (length(c) - 1) / 2;
    res = norm(c(mp+1:end)) / norm(c);
    
    % fprintf('Deg = %d, Res = %e\n', (length(c) - 1) / 2, res);
end

% Use the Chebyshev coefficients to evaluate the function at the prescribed
% points; we do it the Clenshaw way.
row = clenshaw(c, 2 * (pts - a) / (b-a) - 1);


end


