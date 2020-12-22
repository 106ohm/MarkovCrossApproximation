function v = example_10_12_fiber(j, i, t, lambda, Q, pi0)
%EXAMPLE_10_12_FIBER 

r = [ 1 ; 1 ; 0 ];

switch j
    case 1
        % fiber in time
        w = KolmogorovODE(Q(lambda(i(2))), pi0, t);
        v = w * r;
    case 2
        % theta 1
        v = zeros(length(lambda), 1);
        for jj = 1 : length(v)
            w = KolmogorovPoint(Q(lambda(jj)), pi0, t(i(1)));
            v(jj) = dot(r, w);
        end
    case 3
        % component of pi
        v = KolmogorovPoint(Q(lambda(i(2))), pi0, t(i(1)));
end
    

end

