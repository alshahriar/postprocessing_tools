function [psi, sigma, phi] = opt_forcingD(A, Q, omg);

[F, flag] = chol(Q');
Fi = F\eye(size(F,1));

[psi, sigma, phi] = svd(F*inv(-1i*omg*eye(size(A,1))-A)*Fi, 0);
psi = Fi*psi;
phi = Fi*phi;


end