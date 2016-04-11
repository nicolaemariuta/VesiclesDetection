% plot 3D ellipsoid
% developed from the original demo by Rajiv Singh
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/42966
% 5 Dec, 2002 13:44:34
% Example data (Cov=covariance,mu=mean) is included.

Cov = [10 1.8
       1.8 2];
mu = [10 2];

% steps = 100; % pick the number of ellipse point yourself 
% M = 2; % pick the Cov scale by yourself
% [x y] = fn_plot_covariance_2D(mu,Cov,steps,M);

% by default M = 1 and steps = 100;
[x y] = fn_plot_covariance_2D(mu,Cov);

d = mvnrnd(mu, Cov, 100);
figure; plot(d(:,1), d(:,2),'r*'); hold on;
plot(x,y,'b-'); daspect([1 1 1]);
