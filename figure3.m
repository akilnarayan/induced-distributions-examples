% Generates polynomial approximations to a 2-dimensional function, and computes errors.

clear
close all

addpath('induced-distributions');
% Define polynomial family:
rho = 0;
recurrence = @(NN) hermite_recurrence(NN, rho);
d = 2;

% Define approximation space
kmax = 25;
kchoice = 20;
maxindices = total_degree_indices(d, kmax);
indices = total_degree_indices(d, kchoice);
N = size(maxindices, 1);
sampling_ratio = 4;

% Initialized induced sampler from the above
univ_inv = @(uu,nn) fidistinv_freud(uu, nn, 2,rho);
idist_sampler = @(MM) idist_mixture_sampling(MM, maxindices, univ_inv);

% Define function:
center = [0.2 -0.1];
assert(d == numel(center))
%f = @(xx) bump_function_1d(sqrt(sum(xx/(sqrt(kmax)) - repmat(center, [size(xx,1) 1]), 2).^2));
f = @(xx) bump_function_1d(sqrt(sum(0.25*xx - repmat(center, [size(xx,1) 1]), 2).^2));
f2 = @(xx) f(xx).^2;

% Use p-point tensorized Gauss rule to compute "exact" projection coefficients
p = 400;
c_exact = compute_projection_error(f, maxindices, recurrence, p);
f_l2norm_squared = compute_projection_error(f2, zeros([1 d]), recurrence, p)/mpoly_eval(zeros([1 d]), zeros([1 d]), recurrence);

% Number of batches
B = 100;

% First test: for various N, choose M = k*N samples, and do least-squares
proj_errors = zeros([kmax 1]);
sampling_errors = zeros([kmax B]);
Ns = zeros([kmax 1]);
for k = 1:kmax
  indices = total_degree_indices(d, k);
  N = size(indices, 1);
  Ns(k) = N;
  proj_errors(k) = sqrt(f_l2norm_squared - norm(c_exact(1:N)).^2);

  % Now construct B least-squares estimators
  for b = 1:B
    M = ceil(sampling_ratio * N);
    x = randn([M d]) * 1/sqrt(2);
    V = mpoly_eval(x, indices, recurrence);
    w = ones([M 1]);
    V = repmat(w, [1 size(V,2)]).*V;
    sampling_errors(k, b) = norm(c_exact(1:N) - V\(w.*f(x)));
  end
end

%% Second test: fix N and construct errors for various M
% Number of samples
indices = total_degree_indices(d, kchoice);
N = size(indices, 1);
Ms = round(linspace(N, 10*N, 100));

% L2 errors
errors_classical = zeros([numel(Ms) B]);
errors_induced = zeros([numel(Ms) B]);

% Compute estimators
Mcounter = 1;
for M = Ms

  for b = 1:B

    % Gaussian sampling on R^d, variance 1/2
    x = randn([M d]) * 1/sqrt(2);
    V = mpoly_eval(x, indices, recurrence);
    w = ones([M 1]);
    V = repmat(w, [1 size(V,2)]).*V;
    errors_classical(Mcounter, b) = norm(c_exact(1:N) - V\(w.*f(x)));

    % Inducd sampling on R^d
    errors_induced(Mcounter, b) = norm(c_exact(1:N) - V\(w.*f(x)));
    x = idist_sampler(M);
    V = mpoly_eval(x, indices, recurrence);
    w = 1./sqrt(sum(V.^2, 2));
    V = repmat(w, [1 size(V,2)]).*V;
    errors_induced(Mcounter, b) = norm(c_exact(1:N) - V\(w.*f(x)));
   
  end

  Mcounter = Mcounter + 1;
  fprintf('M = %d finished\n', M);

end


% Normalize errors (make relative error)
errors_induced = errors_induced/norm(c_exact);
errors_classical = errors_classical/norm(c_exact);

%%% Visualization
axis_props = {'fontsize', 16, 'fontweight', 'b'};
set_latex = @(hh) set(hh, 'interpreter', 'latex');

Msplot = Ms/size(total_degree_indices(d, kchoice),1);

figure;
subplot(1,3,1);
alh1 = semilogy(Ns, sampling_errors, 'b.'); hold on;
alh2 = semilogy(Ns, proj_errors, 'k', 'linewidth', 3);
N = size(total_degree_indices(d,kchoice), 1);
plot([N N], [1e-4, 1e12], 'k:', 'linewidth', 2);
set_latex(xlabel('$N$'));
set_latex(legend([alh1(1), alh2], '$\|f - g_N\|$', '$\|f - f_N\|$'));
set_latex(title('Gaussian measure ($w$) sampling'));
legend boxoff
set(gca, axis_props{:});
axis([0 max(Ns) 1e-4 1e12]);

subplot(1,3,2);
%semilogy(repmat(Msplot', [1 B]), errors_classical, 'b.'); 
semilogy(repmat(Msplot', [1 B]), errors_classical./proj_errors(kchoice), 'b.'); 
hold on;
%alh = plot([Msplot(1) Msplot(end)], [proj_errors(kchoice) proj_errors(kchoice)], 'k:', 'linewidth', 3);
alh = plot([Msplot(1) Msplot(end)], [1 1], 'k:', 'linewidth', 3);
set_latex(xlabel('$M/N$'));
h = ylabel('$\eta_N$');
set(h, 'rotation', 0);
set_latex(h);
set_latex(title('Gaussian measure ($w$) sampling'));
%set_latex(legend(alh, 'Best approximation error'));
%legend boxoff
set(gca, axis_props{:});
ymax = axis;
ymax = max(ymax(3:4));


subplot(1,3,3);
%semilogy(repmat(Msplot', [1 B]), errors_induced, 'r.');
alh = plot([Msplot(1) Msplot(end)], [1 1], 'k:', 'linewidth', 3);
set(gca, 'yscale', 'log');
hold on;
semilogy(repmat(Msplot', [1 B]), errors_induced./proj_errors(kchoice), 'r.');
%alh = plot([Msplot(1) Msplot(end)], [proj_errors(kchoice) proj_errors(kchoice)], 'k:', 'linewidth', 3);
set_latex(xlabel('$M/N$'));
h = ylabel('$\eta_N$');
set(h, 'rotation', 0);
set_latex(h);
set_latex(title('Induced measure ($\rho$) sampling'));
%set_latex(legend(alh, 'Best approximation error'));
%legend boxoff
set(gca, axis_props{:});
ymin = axis;
ymin = min(ymin(3:4));
temp = axis;
axis([Msplot(1) Msplot(end) ymin ymax]);

subplot(1,3,3);
temp = axis;
axis([Msplot(1) Msplot(end) 1e-1 1e16]);
subplot(1,3,2);
temp = axis;
axis([Msplot(1) Msplot(end) 1e-1 1e16]);
