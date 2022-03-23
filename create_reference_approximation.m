function [RR, err] = create_reference_approximation(Afiber, intervals, U)
    err = 0;
    n = cellfun(@(i) length(i), intervals);
    RR = zeros(n);
    
    t = intervals{1};
    indices = cell(1, length(intervals)-1);
    for i = 1 : length(indices)
        indices{i} = 1 : length(intervals{i+1});
    end
    [indices{:}] = ndgrid(indices{1:end});
    
    for i = 1 : numel(indices{1})
        v = cellfun(@(x) x(i), indices);
        cv = num2cell(v);
        RR(:, cv{:}) = Afiber(1, [1, v]);

        % Construct the approximate sol
        as = zeros(length(t), 1);
        for k = 1 : length(t)
            as(k) = 0;
            for kk = 1 : size(U{1}, 2)
                xx = U{1}(k, kk);
                for jj = 1 : length(indices)
                    xx = xx * U{jj+1}(indices{jj}(i), kk);
                end
                as(k) = as(k) + xx;
            end
        end

        err = hypot(err, norm( RR(:, cv{:}) - as, 'fro' ) );
    end
   
    err = err / norm(RR(:));
end