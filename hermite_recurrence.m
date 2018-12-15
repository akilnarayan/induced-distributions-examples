function[a,b] = hermite_recurrence(N, rho)
% [a,b] = hermite_recurrence(N)
%
%     Returns the recurrence coefficients with indices n for the Hermite
%     polynomials, orthogonal under the weight
%
%      w(x) = C x^(2*rho) * exp(-x^2)
%
%     where C is a normalization constant so that w is a probability density on
%     R.

N = max(N(:));
n = (1:N).' - 1;

a = zeros(size(n));
b = zeros(size(n));

neq0 = (n==0);
nodd = boolean(mod(n,2));

% Un-normalized
%b(neq0) = gamma(rho+1/2);

% Normalized
b(neq0) = 1;

b(~neq0) = 1/2*(n(~neq0));
b(nodd) = b(nodd) + rho;
