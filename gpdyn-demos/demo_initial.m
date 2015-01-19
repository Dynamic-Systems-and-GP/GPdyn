%% generate data
% training
n = 20;
rand('state',18);
randn('state',20);

cov = {'covSum', {'covSEiso','covNoise'}};

loghyper = [log(1.0); log(1.0); log(0.1)];

inf = @infExact;
lik = @likGauss;
mean = @meanZero;


x = 15*(rand(n,1)-0.5);
y = chol(feval(covfunc{:}, loghyper, x))'*randn(n,1);
% valid
xstar = linspace(-7.5, 7.5, 201)';

covfunc = {'covSEiso'};

%% gp - logarithmic scale
% train
hyp0 = gp_initial([], inf, mean, cov, lik, x, y);
hyp = minimize(hyp0, @gpx, -100, inf, mean, cov, lik, x, y);
% predict
[mu, S2] = gpx(hyp, inf, mean, cov, lik, x, y, xstar);
S2 = S2 - exp(2*hyp.lik);
% plot
figure;
f = [mu+2*sqrt(S2);flipdim(mu-2*sqrt(S2),1)];
fill([xstar; flipdim(xstar,1)], f, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
hold on
plot(xstar,mu,'k-','LineWidth',2);
plot(x, y, 'k+', 'MarkerSize', 17);

%% gp - linear scale
% train
bounds = [exp(-8) exp(7)];
hyp0 = gp_initial(bounds, inf, mean, cov, lik, x, y, 0); % exp(-8) ~= 0.0003, exp(7) ~= 1096.6
hyp = minimize(hyp0, @gpx, -100, inf, mean, cov, lik, x, y);
% predict
[mu, S2] = gpx(hyp, inf, mean, cov, lik, x, y, xstar);
S2 = S2 - exp(2*hyp.lik);
% plot
figure;
f = [mu+2*sqrt(S2);flipdim(mu-2*sqrt(S2),1)];
fill([xstar; flipdim(xstar,1)], f, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
hold on
plot(xstar,mu,'k-','LineWidth',2);
plot(x, y, 'k+', 'MarkerSize', 17);