function [x_sol, niters] = Method_of_Steepest_Descent(A, b)

    % set up problem parameters
    k = 0;
    x = ones(size(A, 1), 1);

    tic
    
    r = b - A*x;
    while norm(r) > 1e-6 * norm(b)
        p = r;
        q = A * p;
        alpha = (p'*r) / (p' * q);
        x = x + alpha * p;
        r = r - alpha * q;
        k = k + 1;
    end
    toc
    niters = k;
    x_sol = x;
end