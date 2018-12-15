% Generates plots of induced distribution pdfs from various 1d distributions

clear
close all 

addpath('induced-distributions');

% Define subspace
N = 20;

% Visualization parameters
M = 301;
x = zeros([M 3]);
w = zeros([M 3]);
rho = zeros([M 3]);
rhoinf = zeros([M 3]);

%%%% w: Uniform
id = 1;
x(:,id) = linspace(-1, 1, M).';
w(:,id) = 0.5 * ones(size(w(:,id)));
rhoinf(:,id) = 1/pi * (1 - x(:,id).^2).^(-1/2);
% rhoinf is singular at boundaries, so just truncate there
rhoinf(abs(x(:,id))==1) = max(rhoinf(isfinite(rhoinf(:,id)),id));

% Define polynomial family:
alph = 0;
bet = 0;
[a,b] = jacobi_recurrence(N, alph, bet);
b(1) = 1; % Ensure orthonormality under a probability distribution
V = poly_eval(a, b, x(:,id), N-1);
rho(:,id) = 1/N * sum(V.^2, 2) .* w(:,id);

%%%% w: One-sided exponential 
id = 2;
x(:,id) = linspace(0, 4.5*N, M).';
w(:,id) = exp(-x(:,id));
flags = x(:,id) <= 4*N;
rhoinf(flags,id) = 1/(2*pi*N) * sqrt((4*N - x(flags,id))./x(flags,id));

% Define polynomial family:
alph = 1;
p = 0;
[a,b] = hfreud_recurrence(N, alph, p);
b(1) = 1; % Ensure orthonormality under a probability distribution
V = poly_eval(a, b, x(:,id), N-1);
rho(:,id) = 1/N * sum(V.^2, 2) .* w(:,id);

%%%% w: Gaussian
id = 3;
x(:,id) = linspace(-sqrt(2.5*N), sqrt(2.5*N), M).';
w(:,id) = 1/sqrt(pi) .* exp(-x(:,id).^2);
flags = abs(x(:,id)) <= sqrt(2*N);
rhoinf(flags,id) = 1/(pi*N) * sqrt(2*N - x(flags,id).^2);

% Define polynomial family:
p = 0;
[a,b] = hermite_recurrence(N, p);
b(1) = 1; % Ensure orthonormality under a probability distribution
V = poly_eval(a, b, x(:,id), N-1);
rho(:,id) = 1/N * sum(V.^2, 2) .* w(:,id);


%%%%% Visualization
styles = {'b--', 'r', 'k:'};
lineprops = {'linewidth', 3};
axis_props = {'fontsize', 16, 'fontweight', 'b'};
set_latex = @(hh) set(hh, 'interpreter', 'latex');
ymaxes = [1.5, 1, 0.6];
ymins = [0, 1e-10, 0];
scales = {'linear', 'log', 'linear'};
titles = {'$w$ uniform on $[-1,1]$', '$w$ exponential on $[0, \infty)$', '$w$ Gaussian on $(-\infty,\infty)$'};

for q = 1:3 

  subplot(1,3,q);
  if strcmpi(scales{q},'linear')
    set(plot(x(:,q), w(:,q), styles{1}), lineprops{:});
    hold on;
    set(plot(x(:,q), rho(:,q), styles{2}), lineprops{:});
    set(plot(x(:,q), rhoinf(:,q), styles{3}), lineprops{:}, 'linewidth', 2);
  elseif strcmpi(scales{q}, 'log');
    set(semilogy(x(:,q), w(:,q), styles{1}), lineprops{:});
    hold on;
    set(semilogy(x(:,q), rho(:,q), styles{2}), lineprops{:});
    set(semilogy(x(:,q), rhoinf(:,q), styles{3}), lineprops{:}, 'linewidth', 2);
  else
    assert false
  end
  set_latex(xlabel('$x$'));
  axis([min(x(:,q)), max(x(:,q)), ymins(q), ymaxes(q)]);
  set(gca, axis_props{:});
  set_latex(title(titles{q}));

  if q == 1
    set(legend('$w$', '$\rho$', '$\rho_\infty$'), 'interpreter', 'latex', axis_props{:});
    legend boxoff;
  end

end
