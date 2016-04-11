%read images and build block
N = 256;
istart = 170;
iend = 239;
I = zeros(N,N,iend-istart+1);

for i = istart:iend
    tmp = imread(sprintf('Images/Syn1sec%d.tif',i));
    I(:,:,1+i-istart) = imresize(tmp,[N,N]);
end

I = correctSlice(I,37);


%coordinates for point inside the vescicle
x = 102;
y = 120;
z = 35;
%set of points for calculating covariance
data = [102 121 35; 203 122 6; 104 220 67; 5 121 38; 301 20 44; 102, 120, 35];

% Create some random data
s = [2 2 2];
x = randn(160,1);
y1 = normrnd(s(1).*x,1);
y2 = normrnd(s(2).*x,1);
y3 = normrnd(s(3).*x,1);
data = [y1 y2 y3];






%generate ellipsoid
[ex,ey,ez] = fitElipsoid(data);




%figure;
%surf(ex, ey, ez);
%axis equal;



%display slices
cd imStacks;
OrthoSlicer3d(I);
hold on; surf(data);surf(ex,ey,ez); hold off;

figure;
% now plot the rotated ellipse
% sc = surf(x,y,z); shading interp; colormap copper
h = surfl(ex, ey, ez); colormap copper
title('actual ellipsoid represented by mu and Cov')
axis equal
alpha(0.7)