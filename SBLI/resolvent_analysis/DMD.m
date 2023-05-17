% DMD
function [lambda, psi, phi, b] = DMD(X,Y,dt,r)

nt = size(X,2);
[U,s,V]= svd(X, 0);
r = min([r, size(U,2)]);
Ur = U(:,1:r);
Sr = s(1:r,1:r);
Vr = V(:,1:r);

At = (Ur'*Y)*(Vr/Sr);

[rho, W, Wadj] = eigen_dual(At, eye(r));

psi = Y*(Vr/Sr*W);
phi = Ur*Wadj;

for i=1:r
    psi(:,i) = psi(:,i)/sqrt(psi(:,i)'*psi(:,i));
    phi(:,i) = phi(:,i)/sqrt(phi(:,i)'*phi(:,i));
    psi(:,i) = psi(:,i)/(phi(:,i)'*psi(:,i));
end

b = psi\X(:,1);
lambda = log(rho)/dt;

end