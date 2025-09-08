import casadi.*

n_obj = 200;                % number of variables
w = SX.sym('w', n_obj, 1);   % decision variables
    r = w(1:end-1).*w(2:end);  % pairwise products
    f = sum(r.^2) + sum(r.^4); % L2 + L4 contributions
H = hessian(f, w);                 % symbolic Hessian
H_fun = Function('H_fun', {w}, {H});

% evaluate at a dummy point
w0 = 2*ones(n_obj,1);
H_val = full(H_fun(w0));
max(H_val(:))

% plot sparsity pattern
figure; spy(H_val);
