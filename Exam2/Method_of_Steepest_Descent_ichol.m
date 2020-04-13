function [x_sol, niters] = Method_of_Steepest_Descent_ichol(A, b)
    
    % set up problem parameters
    k = 0;
    % L = sparse(ichol(sparse(A)));
    L = sparse( ichol(sparse(A), struct('type','ict','droptol',1e-9,'michol','off')));
    x = ones(size(A, 1), 1);

    tic
    r = b - A*x;

    while norm(r) > 1e-6 * norm(b)
        
        p = L' \ (L \ r);
        q = A * p;
        alpha = (p' * r) / (p' * q);
        x = x + alpha * p;
        r = r - alpha * q;
        k = k + 1;

    end
    toc
    niters = k;
    x_sol = x;
end


