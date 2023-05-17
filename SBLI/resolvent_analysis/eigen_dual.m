function [lambda, V, W] = eigen_dual(A,Q)

Aadj = adjointQ(A, Q);

[V, lambda] = eig(A);
[W, lambdab] = eig(Aadj);

lambda=diag(lambda); lambdab=diag(lambdab);

[~, sort_lambda] = sort(-real(log(lambda)));
[~, sort_lambdab] = sort(-real(log(lambdab)));

V=V(:,sort_lambda); lambda=lambda(sort_lambda);
W=W(:,sort_lambdab); lambdab=lambdab(sort_lambdab);

for i = 1:size(V, 2)
    V(:,i) = V(:,i)/(W(:,i)'*Q*V(:,i));
end

end