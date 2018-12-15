% Generates polynomial approximations to an n-dimensional function, and
% computes errors.
% Computes hyperbolic cross approximations

clear
close all

cd('/Users/akil/work/papers/mine/2018-lsq-induced/code');

addpath('induced-distributions');
% Define polynomial family:
rho = 0;
recurrence = @(NN) hermite_recurrence(NN, rho);
ds = [4 8 20];
%ds = [1 2 3];

% Define approximation space
%ks = [50 30 15];
ks = [20 10 5];

% Ratios of sample counts
ratio_min = 1;
ratio_max = 5;
total_ratios = 100;

B = 100; % Number of trials

% Initialized induced sampler from the above
univ_inv = @(uu,nn) fidistinv_freud(uu, nn, 2,rho);

% Define function:
%center = [0.2 -0.1 -0.2 0.1];
%assert(d == numel(center))
f = @(xx) prod(exp(-xx.^2./repmat(1:size(xx,2), [size(xx,1) 1])),2);
f1d = @(xx,dd) exp(-xx.^2/dd);

% L2 errors
errors_classical = zeros([total_ratios B length(ds)]);
errors_induced = zeros([total_ratios B length(ds)]);

for dcounter = 1:numel(ds)
  d = ds(dcounter);

  indices = hyperbolic_cross_indices(d, ks(dcounter));
  idist_sampler = @(MM) idist_mixture_sampling(MM, indices, univ_inv);
  N = size(indices, 1);

  % Use p-point tensorized Gauss rule to compute "exact" projection coefficients
  p = 300;
  [a,b] = recurrence(p+1); b(1) = 1;
  [xx,ww] = gauss_quadrature(a, b, p);
  c_exact = ones([size(indices,1) 1]);
  f_l2norm_squared = 1;

  % Since f is a product of 1d functions, can do this dimension-by-dimension
  f1ds = zeros([p d]);
  for qd = 1:d
    f1ds(:,qd) = f1d(xx,qd);

    f_l2norm_squared = f_l2norm_squared * (ww.'*f1d(xx,qd).^2);
  end
  V = poly_eval(a, b, xx, max(indices(:)));
  for qq = 1:size(indices,1)
    for qd = 1:d
      c_exact(qq) = c_exact(qq) * (ww.'*(f1ds(:,qd).*V(:,1+indices(qq,qd))));
    end
  end

  %f_l2norm_squared = compute_projection_error(f2, zeros([1 d]), recurrence, p)/mpoly_eval(zeros([1 d]), zeros([1 d]), recurrence);

  proj_error(dcounter) = sqrt(f_l2norm_squared - norm(c_exact).^2);

  % Number of samples
  Ms = round(linspace(N*ratio_min, N*ratio_max, total_ratios));

  % Compute estimators
  Mcounter = 1;
  for M = Ms

    for b = 1:B

      % Gaussian sampling on R^d, variance 1/2
      x = randn([M d]) * 1/sqrt(2);
      V = mpoly_eval(x, indices, recurrence);
      w = ones([size(x,1) 1]);
      V = repmat(w, [1 size(V,2)]).*V;
      errors_classical(Mcounter, b,dcounter) = norm(c_exact - V\(w.*f(x)));

      % Inducd sampling on R^d
      x = idist_sampler(M);
      V = mpoly_eval(x, indices, recurrence);
      w = 1./sqrt(sum(V.^2, 2));
      V = repmat(w, [1 size(V,2)]).*V;
      errors_induced(Mcounter, b,dcounter) = norm(c_exact - V\(w.*f(x)));
     
    end

    Mcounter = Mcounter + 1;
    fprintf('M = %d finished\n', M);

  end

  % Normalize errors (make relative error)
  errors_induced(:,:,dcounter) = errors_induced(:,:,dcounter)/proj_error(dcounter);
  errors_classical(:,:,dcounter) = errors_classical(:,:,dcounter)/proj_error(dcounter);

end

clear V x

save test_approx_hermite_nd.mat
load test_approx_hermite_nd.mat

visualize = false;

%%%% Visualization
if visualize

  axis_props = {'fontsize', 16, 'fontweight', 'b'};
  set_latex = @(hh) set(hh, 'interpreter', 'latex');

  Msplot = Ms/N;
  ymin = Inf;
  ymax = 0;
  marksize = 10;

  figure;
  for dcounter = 1:length(ds);
    subplot(1,length(ds),dcounter);
    alh1 = semilogy(repmat(Msplot', [1 B]), errors_classical(:,:,dcounter), 'b.', 'markersize', marksize); 
    hold on;
    alh2 = semilogy(repmat(Msplot', [1 B]), errors_induced(:,:,dcounter), 'r.', 'markersize', marksize); 
    plot([Msplot(1) Msplot(end)], [1 1], 'k:', 'linewidth', 3);
    set_latex(xlabel('$M/N$'));
    set(ylabel('$\eta_N$', 'interpreter', 'latex', 'rotation', 0));
    set_latex(title(['$w$ Gaussian in $d = ' num2str(ds(dcounter)) '$ dimensions']));
    if dcounter == 1
      set_latex(legend([alh1(1),alh2(1)], 'Standard sampling', 'Induced distribution sampling'));
      legend boxoff;
    end
    set(gca, axis_props{:});
    temp = axis;
    ymax = max([ymax temp(3:4)]);
    ymin = min([ymin temp(3:4)]);
  end

  for dcounter = 1:length(ds);
    subplot(1,length(ds),dcounter);
    axis([ratio_min ratio_max ymin, ymax]);
  end

end
