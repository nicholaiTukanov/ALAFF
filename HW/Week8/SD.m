close all

m = 100;
n = 10;

A = generateSPDmatrix(m);
L = sparse(ichol(sparse(A)));
M = L * L';
b = ones(m,1);
x = zeros(m,1);

OSD(A, b, x);
PSD(A, b, x, M);

function A = generateSPDmatrix(n)
    A = rand(n,n);
    A = 0.5*(A+A');
    A = A + n*eye(n);
end

function x_end = OSD(A, b, x_start)
    k = 0;
    tol = 10^(-6);
    tic
    r_cur = b - A*x_start;
    x_cur = x_start;
    res = r_cur' * r_cur;
    while res > tol;
        p = r_cur;
        q = A * p;
        alpha = (res) / (p' * q);
        x_cur = x_cur + alpha * p;
        r_cur = r_cur - alpha * q;
        res = r_cur' * r_cur;
        k = k + 1;
        % next = input( 'press RETURN to continue' );
    end
    toc
    k
    x_end = x_cur;
end

% function [x_end, r_end] = TSD(A, b, x_start, M):
%     r_cur = b - A*x_start;
%     x_cur = x_start;
    
%     while r_cur ~= 0
        
%         p = r_cur;
%         q = A * p;
%         alpha = (p' * r_cur) / (p' * q);
%         x_cur = x_cur + alpha * p;
%         r_cur = r_cur - alpha * q;

%     end

%     x_end = x_cur;
%     r_end = r_cur;

function x_end = PSD(A, b, x_start, M)
    r_cur = b - A*x_start;
    x_cur = x_start;
    k = 0;
    while r_cur ~= 0
        p = M^-1 * r_cur;
        q = A * p;
        alpha = (p' * r_cur) / (p' * q);
        x_cur = x_cur + alpha * p;
        r_cur = r_cur - alpha * q;
        k = k + 1;
    end
    x_end = x_cur;
    k
end