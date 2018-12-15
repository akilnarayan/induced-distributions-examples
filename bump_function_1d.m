function[y] = bump_function(x)
% Evaluates the bump function defined as
%
%       {  exp(-1/(1-x^2)),   |x| < 1
%   y = {
%       {  0                  |x| >= 1
%
% This function is smooth, but not analytic.

flags = (abs(x) < 1);
y = zeros(size(x));

y(flags) = exp(1)*exp(-1./(1-x(flags).^2));
