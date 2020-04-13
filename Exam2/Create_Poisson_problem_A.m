function [ A ] = Create_Poisson_problem_A( N )

  % Create the archtypical matrix A for an N x N Poisson problem (5-point
  % stencil.

  % Set the diagonal

  D = diag(ones(1,N^2)*4);
  
  % Set the entries of the first sub and super diagonals
  sub_super_d = [];

  for k=1:N^2-1
    if mod(k,N) == 0
      sub_super_d = [sub_super_d, 0];
    else
      sub_super_d = [sub_super_d, 1];
    end
  end

  L = diag(sub_super_d, -1) + diag(sub_super_d, 1);

  % Set the other off-diagonal entries

  off_diag = diag(ones(1, N^2 - N), -N) + diag(ones(1, N^2 - N), N);

  A = (D - L - off_diag);
  
end





