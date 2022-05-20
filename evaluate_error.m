function err = evaluate_error(U, Aelem)

npoints = 1000;
err = 0;

n = arrayfun(@(i) size(U{i}, 1), 1 : length(U));
d = length(n);

i = zeros(1, d);

for k = 1 : npoints
    for j = 1 : d
        i(j) = randi(n(j), 1, 1);
    end

    appr = 0;
    for l = 1 : size(U{1}, 2)
        newfactor = 1;
        for ll = 1 : d
            newfactor = newfactor * U{ll}(i(ll), l);
        end

        appr = appr + newfactor;
    end

    err = max(err, abs(Aelem(i) - appr));
end

end