function [x y] = fn_plot_covariance_2D(mu,Sigma,steps,M)

[U,L] = eig(Sigma);
% L: eigenvalue diagonal matrix
% U: eigen vector matrix, each column is an eigenvector

if nargin < 4
    M = 2; % choose your own N
end

if nargin < 3
    steps = 100;
end

% For N standard deviations spread of data, the radii of the eliipsoid will
% be given by N*SQRT(eigenvalues).
radii = M*sqrt(diag(L));

% plot the ellipse
degree = atan(U(2)/U(1))*180/pi;

[x y] = fn_calculateEllipse(mu(1),mu(2),radii(1),radii(2), degree, steps);