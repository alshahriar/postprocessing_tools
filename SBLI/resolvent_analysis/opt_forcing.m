function [psi, sigma, phi] = opt_forcing(lambda,V,omg)
% [psi, sigma, phi] = opt_forcing(lambda,V,Q,omg)

At = diag(lambda);
% Q = V'*Q*V;
Qt = V'*V;

[psi, sigma, phi] = opt_forcingD(At, Qt, omg);

psi = V*psi;
phi = V*phi;
sigma = diag(sigma);

end