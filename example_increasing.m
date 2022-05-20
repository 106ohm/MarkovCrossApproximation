functions = { ...
    @(p) p, ...
    @(p) 1 / p, ...
    @(p) p + sign(p - 0.35) * 0.2 ...
};

% KK is functions{KK}
% [NPARAMS, KK, ERR, STEPS, TOL]
data = zeros(0, 4);

for nparams = [2, 3, 4]
    for kk = 1 : length(functions)
        for tol = [ 1e-2, 1e-4, 1e-6, 1e-8 ]
            f = functions{kk};
            
            %% Model parameters
            nreplicas = 4;
            
            nstates = nreplicas+1; % labelled as n, n-1,..., 1 and 0 (failed system state)
            
            pi0 = zeros(nstates,1);
            pi0(1) = 1;
            
            % theta is a vector of parameters, length(theta)=nparams, passed 
            % entry by entry
            density = 0.8;
            R2 = sprand(nreplicas+1,nreplicas+1,density);
            R2 = R2-diag(diag(R2));
            Q = @(varargin)infgen(R2, f, nreplicas, varargin{:});
            
            %% Model variable parameters (lambda and c)
            tf = 10*nreplicas;
            low = 0.25;
            up = 0.5;
            
            npoints = 50;
            t = linspace(0, tf, npoints);
            % lambda = repmat(linspace(low, up, npoints),nparams,1);
            n = [ length(t), npoints*ones(1,nparams)];
            
            intervals = {};
            intervals{1} = t;
            for p = 1 : nparams
                intervals{p+1} = linspace(low, up, npoints);
            end
            
            npar = prod(n);
            
            %% Reliability at time t
            kind = 'instantaneous';
            r = [ones(nreplicas,1);0];
            
            [Afiber, Aelem] = create_fiber_functions(Q, intervals, pi0, r, tol, kind);
            
            % ACA
            timer_aca = tic;
            U = aca_nd(n, Afiber, tol);
            timer_aca = toc(timer_aca);
            
            % Reference solution
            timer_reference = tic;
            if nparams < 6 && false
                [RR, err] = create_reference_approximation(Afiber, intervals, U);
            else
                err = evaluate_error(U, Aelem);
                % RR = 0; err = 0;
            end
            timer_reference = toc(timer_reference);
            
            fprintf('NPARAMS = %d\n', nparams);
            fprintf('function number %d\n', kk);
    
            fprintf('Reliability, the full tensor has %d entries\n', npar);
            
            fprintf('Time for ACA: %f seconds\n', timer_aca);
            fprintf('Time for explicit computation: %f seconds\n', timer_reference);
            fprintf('Acceleration factor: %2.3fx\n', timer_reference / timer_aca);
            fprintf('Relative error from the reference solution: %e\n', err);
    
            data = [ data ; [ nparams, kk, err, size(U{1}, 2), tol] ];
            
            % save(strcat('example_increasing_RR_',int2str(nparams),'.mat'), "RR");
            % norm(RR(:) - mean(RR(:))) / norm(RR(:))
            
            % % plot Reliability at time tend
            % tend = round(npoints/5);
            % grid = zeros(npoints, npoints);
            % grid(:,:)=RR(tend,:,:);
            % imagesc(lambda,c,grid);
            % colorbar;
            % title('Reliability at $t_{\mathrm{end}}=2$ years','interpreter','latex');
            % xlabel('$\lambda$','interpreter','latex'); 
            % ylabel('c'); 
        end
    end
end

dlmwrite('example_increasing.dat', data, '\t');

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
