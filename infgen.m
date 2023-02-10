%% Infinitesimal Generator Matrix
function Q = infgen(R2, f, nreplicas, varargin)
% birth-death process plus random deaths in group
br = 1;
dr = 1;

% seed = 1;
% rng(seed);

R1 = diag(dr*ones(nreplicas,1),1) + diag(br*ones(nreplicas,1),-1);

% R2 = sprand(nreplicas+1,nreplicas+1,density);
% R2 = R2-diag(diag(R2));
[i,j,~] = find(R2);
for k = 1 : length(i)
    R2(i(k),j(k)) = f(varargin{mod(k, length(varargin)) + 1});
end

R = R1 + R2;

Q = (R-diag(R*ones(nreplicas+1,1)));
end
