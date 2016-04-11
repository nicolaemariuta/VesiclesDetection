function [ RX, RY, RZ, EX,EY,EZ ] = fitElipsoid( data, N )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Calculate the eigenvectors and eigenvalues
covariance = cov(data);
%covariance = covariance *10;

mu = mean(data);

[eigenvec, eigenval ] = eig(covariance);
% eigenvec: eigenvalue diagonal matrix
% eigenval: eigen vector matrix, each column is an eigenvector

% For N standard deviations spread of data, the radii of the eliipsoid will
% be given by N*SQRT(eigenvalues).
%N = 1; 
radii = N*sqrt(diag(eigenval));


% generate coordinates for "unrotated" ellipsoid
RX = radii(1);
RY = radii(2);
RZ = radii(3);
[xc,yc,zc] = ellipsoid(0,0,0,radii(1),radii(2),radii(3));
s = size(xc);




X = eigenvec*[xc(:),yc(:),zc(:)]';
EX = reshape(X(1,:),s)+mu(1);
EY = reshape(X(2,:),s)+mu(2);
EZ = reshape(X(3,:),s)+mu(3);

% rotate data with orientation matrix eigenvec and center mu
%a = kron(eigenvec(:,1),xc);
%b = kron(eigenvec(:,2),yc);
%c = kron(eigenvec(:,3),zc);

%data = a+b+c; n = size(data,2);

%the rotated ellipsoid
%EX = data(1:n,:)+mu(1);
%EY = data(n+1:2*n,:)+mu(2);
%EZ = data(2*n+1:end,:)+mu(3);

end

