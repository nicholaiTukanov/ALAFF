format long

% Create a random m x m tridiagonal matrix.  Even though we will only
% update the lower triangular part, we create the symmetric matrix.
T = rand( m, m );

% Extract upper triangle plus first subdiagonal.
T = triu( T, -1 );

% Set strictly upper triangular part to transpose of strictly lower
% triangular part
T = tril( T ) + tril( T,-1 )'