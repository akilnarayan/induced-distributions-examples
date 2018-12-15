% Generates polynomial approximations to a one-dimensional function, and computes errors.

clear
close all

addpath('induced-distributions');

% Define polynomial family:
alph = 0;
bet = 0;

% Define function:
f = @(xx) bump_function_1d(1.5*xx - 0.2);

% Sizes of Gaussian interpolation grids to plot
Ns = [5 10 20];
M = 2e3;
[a,b] = jacobi_recurrence(max(max(Ns),M)+1, alph, bet);

% Plotting grid
[xp,wp] = gauss_quadrature(a, b, M);

c = zeros([max(Ns) length(Ns)]);

counter = 1;

Vp = poly_eval(a, b, xp, max(Ns)-1);

for N = Ns

  [x,w] = gauss_quadrature(a,b,N);

  c(1:N,counter) = poly_eval(a, b, x, N-1)'*(w.*f(x));

  counter = counter + 1;

end

f_approx = Vp*c;

% Generate L^2 errors for different N's
Nerrors = 5:100;
Vp = poly_eval(a, b, xp, max(Nerrors)-1);
errors = zeros([numel(Nerrors) 1]);
fNerrors = zeros([numel(Nerrors) 1]);
coeff_errors = zeros([numel(Nerrors) 1]);
counter = 1;
for N = Nerrors
  [x,w] = gauss_quadrature(a, b, N);

  c_error = poly_eval(a, b, x, N-1)'*(w.*f(x));
  errors(counter) = sqrt(sum(wp.*(f(xp) - Vp(:,1:N)*c_error).^2));

  c_fN = poly_eval(a, b, xp, N-1)'*(wp.*f(xp));
  fNerrors(counter) = sqrt(sum(wp.*(f(xp) - Vp(:,1:N)*c_fN).^2));

  coeff_errors(counter) = norm(c_fN - c_error);
  counter = counter + 1;
end

%%%%%%%%%%% 
% Plotting

styles = {'b--', 'r', 'k:'};
axis_props = {'fontsize', 16, 'fontweight', 'b'};
set_latex = @(hh) set(hh, 'interpreter', 'latex');

assert(length(Ns) == length(styles));

figure; 
subplot(1,3,1);
hold on;
counter = 1;
legend_texts = {};
for N = Ns
  plot(xp, f_approx(:,counter), styles{counter}, 'linewidth', 3); hold on;
  legend_texts{counter} = sprintf('$N$ = %d', N);
  counter = counter + 1;
end
legend_texts{counter} = 'Exact function';
plot(xp, f(xp), 'k');
set_latex(legend(gca, legend_texts, 'location', 'northwest'));
legend boxoff

set(gca, axis_props{:});
set_latex(xlabel('$x$'));
set_latex(title('Bump function approximation'));

subplot(1,3,2);
semilogy(Nerrors, errors, 'b', Nerrors, fNerrors, 'k:', 'linewidth', 3);
set_latex(xlabel('$M = N$'));
set_latex(ylabel('Error'));
set_latex(title('Error using $N$-point Gaussian quadrature'));
set(gca, axis_props{:});
set_latex(legend(gca, {'$\|f - g_N\|$', '$\|f - f_N\|$'}));
legend boxoff

subplot(1,3,3);
plot(Nerrors, coeff_errors./fNerrors, 'b', 'linewidth', 3);
set_latex(xlabel('$M = N$'));
set_latex(title('Error'));
set(gca, axis_props{:});
set_latex(legend(gca, {'$\|g_N - f_N\| / \| f - f_N \|$'}));
legend boxoff

