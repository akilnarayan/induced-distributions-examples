function[c] = compute_projection_error(f, indices, recurrence, p);
% Computes projection coefficients for the function f associated to the
% multivariate indices "indices" using p-point Gauss quadrature rules in each
% dimension.
%
% recurrence is a function handle with the calling syntax recurrence(N) that
% returns the first N recurrence coefficients for each univariate family.

d = size(indices, 2);
N = size(indices, 1);

c = zeros([N 1]);

% Proceed dimension-by-dimension
[a,b] = recurrence(p+1);
[x,w] = gauss_quadrature(a, b, p);

siz = p*ones([d-1 1]);
ind = 1;

subs = cell([1 d-1]);

% Sum over dimension 1, iterating over fixed values of other dimensions
while ind < p^(d-1)

  % Get inices for other dimensions
  [subs{:}] = ind2sub(siz, ind);
  xx = [x repmat(x(cell2mat(subs)).', [p 1])];

  V = mpoly_eval(xx, indices, recurrence);
  c = c + prod(w(cell2mat(subs))) * V'*(w.*f(xx));

  ind = ind + 1;

end
