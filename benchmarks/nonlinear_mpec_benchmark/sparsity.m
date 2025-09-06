import casadi.*

n_obj = 100;                % number of variables
w = SX.sym('w', n_obj, 1);   % decision variables

% case 'CHEBYQAD'
  % % Parameters
  %   H = 1/(n_obj+1);
  %   H2 = H^2;
  %   KAPPA = 1.0;
  %   KAPPAH2 = KAPPA*H2;
  % 
  %   % Residuals (gradient-like)
  %   r = SX.zeros(n_obj,1);
  %   r(1)     = w(1) - 0;               % G(0)
  %   for i=1:n_obj-1
  %       r(i+1) = w(i) - w(i+1);        % G(i)
  %   end
  %   r = r - 2*H2;                       % L(i)
  %   r(n_obj) = r(n_obj) - (1 + 2*H2);   % L(N)
  % 
  %   % Function
  %   f = 0.5 * r.' * r;


f = 0;
for i = 1:n_obj
    % finite-difference like group (difference with next)
    if i < n_obj
        q = w(i+1) - w(i);
    else
        q = w(i);
    end
    % element contribution (scaled sin/cos)
    f = f + 0.5 * q^2 + 100*sin(0.01*w(i)) + 1e8*cos(w(i));
end


H = hessian(f, w);                 % symbolic Hessian
H_fun = Function('H_fun', {w}, {H});

% evaluate at a dummy point
w0 = 1*ones(n_obj,1);
H_val = full(H_fun(w0));

% plot sparsity pattern
figure; spy(H_val);
%%
import casadi.*

N = 10;             % adjust N as needed
KAPPA = 1.0;
OBJSCALE = 1e8;

w = SX.sym('w', N, 1);

% Define individual contributions
S = 100.0*sin(0.01*w);
C = OBJSCALE*cos(w);

% Define groups (HALFL2 structure)
f = 0.5*OBJSCALE*sum((diff([0; w])).^2);  % finite-difference quadratic
f = f + sum(C) + sum(S);                  % add individual element contributions

% Hessian
H_f = hessian(f, w);

% Function to evaluate and plot
H_fun = Function('H_fun', {w}, {H_f});
H_val = full(H_fun(zeros(N,1)));

figure; spy(H_val); title('Hessian sparsity of FLETBV3M');
