function T = cpd2full(U)
%CPD2FULL 

d = length(U);
n = arrayfun(@(i) size(U{i}, 1), 1 : d);
k = size(U{1}, 2);

T = zeros(n);

for s = 1 : k
    w = U{1}(:, s);
    for j = 2 : d
        w = kron(U{j}(:,s), w);
    end
    
    T = T + reshape(w, n);
end

end

