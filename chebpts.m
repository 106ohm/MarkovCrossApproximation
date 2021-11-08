function [x, w] = chebpts(n, domain)
%CHEBPTS Get Chebyshev points and integration weights

if ~exist('domain', 'var')
    domain = [-1, 1];
end

% Compute the Chebyshev points for the domain [-1, 1], and rescale them to
% fit the given domain. 
x = cos( (n-1:-1:0) / (n-1) * pi).';
x = domain(1) + (x + 1) / 2 * (domain(2) - domain(1));

if nargout >= 2
    w = ones(1, n) * pi / n;
end

end