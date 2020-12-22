function i = find_pivot(Afiber, i, n, strategy)
%FIND_PIVOT 

d = length(n);

switch strategy
    case 'random'
        for j = 1 : d
            i(j) = randi(n(j), 1, 1);
        end
        
    case 'fiber'
        i(d) = randi(n(d), 1, 1);
        rounds = 1;
        for s = 1 : rounds
            for j = 1 : d
                f = Afiber(j, i);
                [~, i(j)] = max(abs(f));
            end
        end
        
    otherwise
        error('Unsupported selection strategy');
end

end

