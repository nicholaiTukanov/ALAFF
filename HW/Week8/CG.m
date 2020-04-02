% name: Nicholai Tukanov
% eid:  nst383

close all

m = 100;

A = generateSPDmatrix(m);
L = generatePreconditioner(A);
M = L * L';
b = ones(m,1);

OCG(A, b);
PCG(A, b, M);
pcg(A, b, 1e-13, 500, M);

% Generate a n x n SPD matrix
% Code was taken from StackOverflow
function A = generateSPDmatrix(n)
    A = rand(n,n);
    A = 0.5*(A+A');
    A = A + n*eye(n);
end

% Generate a preconditioner for a given A using incomplete cholesky 
function L = generatePreconditioner(A)
    L = sparse(ichol(sparse(A)));
end

% Ordinary Conjugate Gradient
% Code is based off algorithm given by 8.3.5.2
function x_end = OCG(A, b)
    disp("Ordinary Conjugate Gradient Method"+ newline)
    k = 0;

    tic
    x = zeros(size(A, 1), 1);
    r = b;
    res = r' * r;
    while res > 1e-13
        if k == 0
            p = r;
        else
            gam = (-p'*A*r) / (p'*A*p);
            p = r + gam*p;
        end

        alpha = (p' * r) / (p' * A * p);
        x = x + alpha * p;
        r = r - alpha * A * p;
        res = r' * r;
        k = k + 1;
        
    end
    k
    res
    toc
    x_end = x;
end

% Preconditioned Conjugate Gradient
% Used Conjugate Gradient wiki page as a source
function x_end = PCG(A, b, M)
    disp("Preconditioned Conjugate Gradient Method"+ newline);
    k = 0;

    tic
    x = zeros(size(A, 1), 1);
    r = b;
    res = r' * r;
    while res > 1e-13
        z = M^-1*r;
        if k == 0
            p = z;
        else
            gam = (-p'*z) / (p'*A*p);
            p = z + gam*p;
        end

        alpha = (p' * r) / (p' * A * p);
        x = x + alpha * p;
        r = r - alpha * A * p;
        res = r' * r;
        k = k + 1;
    end
    k
    res
    toc
    x_end = x;

end