clear all;
close all;
clf;

% This is the linear reference model used for the state predictor.
Am = [0 1; -1 -1.4];
b = [0 1]';
Gamma = 100000;

% This is the P matrix used to ensure that the adaptive updates happen
% in the right direction.  The point-wise normalization doesn't change
% its nature, I think SLICOT just goes the other direction when it
% solves the lyapunov equation.
P = lyap(Am, eye(rank(Am))) .* [1 -1; -1 1];

% ksp0 = 1.4 * sqrt(Gamma) - 1;

% Using Ksp here yields modified prediction error dynamics:
%
% xtildedot = (Am - Ksp) * xtilde + b * sigmatilde
%
% Ksp lets us move the poles of this system to get better damping
% on parameter estimates.
Ksp = [0 0; 0 10];

s = tf('s');
f = @(G, A, K, b, P) G * (inv(s * eye(rank(A)) - A + K)*b)' * P * b;

% Loop shaping
F = f(Gamma, Am, Ksp, b, P);
adap = F / (s + F);

% No loop shaping
FF = f(Gamma, Am, zeros(rank(Am)), b, P);
adap0 = FF / (s + FF);

% This shows the transfer function from sigma to sigmahat with the
% loop-shaping matrix Ksp.
figure(1);
bode(adap);

% This shows the transfer function from sigma to sigmahat without the
% loop-shaping matrix Ksp in place.  Note that it's underdamped.
figure(2);
bode(adap0);
