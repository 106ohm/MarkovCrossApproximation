function v = aca_eval_fiber(U, j, i)
%ACA_EVAL_FIBER 

n = arrayfun(@(i) size(U{i}, 1), 1 : length(U));
d = length(n);

v = zeros(n(j), 1);

for s = 1 : size(U{1}, 2)
    w = U{j}(:,s);
    for jj = 1 : d
        if jj ~= j
            w = w * U{jj}(i(jj), s);
        end
    end
    
    v = v + w;
end

end

