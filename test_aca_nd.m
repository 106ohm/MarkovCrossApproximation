k = 2;
n = [ 40, 25, 32, 40 ];
d = length(n);

UU = cell(1, d);

for j = 1 : d
    UU{j} = rand(n(j), k);
end

T = cpd2full(UU);

Afiber = @(j,i) aca_eval_fiber(UU, j, i);

U = aca_nd(n, Afiber, 1e-6);