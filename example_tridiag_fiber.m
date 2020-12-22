function v = example_tridiag_fiber(j, i, t, lambda1, lambda2, m, Q, pi0)
%EXAMPLE_10_4_FIBER 

% comp = 4;
r = ones(m, 1); r(end) = 0;

switch j
    case 1
        % fiber in time
        w = KolmogorovODE(Q([ lambda1(i(2)), lambda2(i(3)) ], m), pi0, t);
        v = w * r;
    case 2
        % theta 1
        v = zeros(length(lambda1), 1);
        for jj = 1 : length(v)
            w = KolmogorovPoint(Q([ lambda1(jj), lambda2(i(3)) ], m), pi0, t(i(1)));
            v(jj) = dot(r, w);
        end
    case 3
        % theta 2
        v = zeros(length(lambda2), 1);
        for jj = 1 : length(v)
            w = KolmogorovPoint(Q([ lambda1(i(2)), lambda2(jj) ], m), pi0, t(i(1)));
            v(jj) = dot(r, w);
        end
end
    

end

