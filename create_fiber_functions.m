function Afiber = create_fiber_functions(Q, intervals, pi0, r, tol, kind)
%CREATE_FIBER_FUNCTIONS

t = intervals{1};

[rr,pp,kk] = rational(min(14, 5 + ceil(-log10(tol))), 1);

Afiber = @(j,i) fiber_evaluator(j, i, t, intervals(2:end), pi0, Q, r, rr, pp, kk, kind);

end

function v = fiber_evaluator(j, i, t, theta, pi0, Q, r, rr, pp, kk, kind)

if j == 1
    % fiber in time
    
    % Construct the arguments where Q is evaluated
    args = cell(1, length(theta));
    for k = 1 : length(args)
        args{k} = theta{k}(i(k+1));
    end
    
    switch kind
        case 'instantaneous'
            w = KolmogorovODE(Q(args{:}), pi0, t);
        case 'accumulated'
            w = KolmogorovIntegralODE(Q(args{:}), pi0, t);
    end
    
    v = w * r;
else
    args = cell(1, length(theta));
    for k = 1 : length(args)
        args{k} = theta{k}(i(k+1));
    end
    
    % We note that the parameter with index j is actually args{j-1},
    % because we have removed the time. 
    switch kind
        case 'instantaneous'
            v = rowChebWrapper(@(x) dot(r, KolmogorovPoint(...
                Q(args{1:j-2}, x, args{j:end}), pi0, t(i(1)), rr, pp, kk)), ...
                theta{j-1})';
        case 'accumulated'
            v = rowChebWrapper(@(x) dot(r, KolmogorovIntegralPoint(...
                Q(args{1:j-2}, x, args{j:end}), pi0, t(i(1)), rr, pp, kk)), ...
                theta{j-1})';
    end
end

end
