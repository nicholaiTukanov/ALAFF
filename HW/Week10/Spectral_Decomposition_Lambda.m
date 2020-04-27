function Lambda = Spectral_Decomposition_Lambda( T )

    tic
    n = size(T,1);
    Lambda = [];

    for k = n:-1:3
        itr = 0;
        while abs( T(k,k-1) ) > 1e-10
            T = Francis_Step( T );
            itr = itr + 1;
        end
        Lambda = [ Lambda T(k,k) ];
        T = T(1:k-1,1:k-1);
    end
    Lambda = diag([ Lambda eig(T)'])
    toc

end