function [x_sol, niters] = PCG(A, b)
    
    % set up problem parameters
    % L = sparse(ichol(sparse(A)));
    % L = sparse( ichol(sparse(A), struct('type','ict','droptol',1e-2,'michol','off')));
    L = sparse( ichol(sparse(A), struct('type','ict','droptol',1e-9,'michol','off')));
    x = ones(size(A, 1), 1);
    k = 0;

    tic

    r0 = b - A*x;
    z0 = L' \ (L \ r0);
    p = z0;

    res = norm(r0);

    while res > 1e-6 * norm(b)

        q = A*p;
        alpha = (r0' * z0) / (p' * q);
        x = x + alpha*p;
        r = r0 - alpha*q;

        z = L' \ (L \ r);
        gam = (r' * z) / (r0' * z0);
        p = z + gam*p;

        r0 = r;
        z0 = z;
        res = norm(r);
        k = k + 1;    
    end
    
    toc
    niters = k;
    x_sol = x;
end


