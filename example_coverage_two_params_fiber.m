function v = example_battery_two_params_fiber(j, i, t, lambda1, lambda2, pi0, Q, r, rr, pp, kk)

switch j
    case 1
        % fiber in time
        w = KolmogorovODE(Q(lambda1(i(2)), lambda2(i(3))), pi0, t);
        % w = KolmogorovIntegralODE(Q(lambda1(i(2)), lambda2(i(3))), pi0, t);
        v = w * r;
    case 2
        % theta 1
        v = rowChebWrapper(@(x) dot(r, KolmogorovPoint2(Q(x, lambda2(i(3))), pi0, t(i(1)), rr, pp, kk)), lambda1)';
        % v = w' * r;
    case 3
        % theta 2
        v = rowChebWrapper(@(x) dot(r, KolmogorovPoint2(Q(lambda1(i(2)), x), pi0, t(i(1)), rr, pp, kk)), lambda2)';
        % v = w' * r;
end

end

