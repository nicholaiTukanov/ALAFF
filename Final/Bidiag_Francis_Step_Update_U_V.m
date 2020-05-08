function [U, B_next, V] = Bidiag_Francis_Step(U_A, B, V_A)


    [m, n] = size(B);
    [m_u,n_u] = size(U_A);
    U = U_A;
    V = V_A;

    % first G comes from B' * B (the values are given in the notes) 
    G = Givens_rotation([B(1,1)^2 - B(m-1,m)^2 - B(m,m)^2  
                            B(1,2)*B(1,1)    ]);
    B(:,1:2) = B(:,1:2) * G;

    blkdiag(G, eye(m_u-2))
    V = V * blkdiag(G, eye(m_u-2));

    F = Givens_rotation( [B(1,1),B(2,1)]' );
    B(1:2,:) = F' * B(1:2,:);

    U = U * blkdiag(F, eye(m_u-2));

    % Only loop if the matrix is bigger than a 2x2
    if m > 2
        for i=1:m-3

            G = Givens_rotation( [B(i,i+1) B(i,i+2)]' );
            B(:,i+1:i+2) = B(:,i+1:i+2) * G;

            blkdiag(eye(i), G, eye(m_u-i-2))
            V = V * blkdiag(eye(i), G, eye(m_u-i-2));

            F = Givens_rotation( [B(i+1,i+1), B(i+2,i+1)]');
            B(i+1:i+2,:) = F' * B(i+1:i+2,:);

            U = U * blkdiag(eye(i), F, eye(m_u-i-2));

        end

        G = Givens_rotation( [B(m-2,m-1), B(m-2,m)]' );
        B(:,m-1:m) = B(:,m-1:m) * G;

        blkdiag(eye(m_u-3), G, eye(1))
        V = V * blkdiag(eye(m_u-3), G, eye(1));

        F = Givens_rotation( [B(m-1,m-1), B(m,m-1)]');
        B(m-1:m,:) = F' * B(m-1:m,:);

        U = U * blkdiag(eye(m_u-3), F, eye(1));

    end

    B_next = B;

end