cd ..
%set parameters
N = 256;
r = 5;
sigma = 2;
T = 20;
alpha = .01;
beta = .05;
istart = 170;
iend = 239;
in = 10;



%read images and build block
I = zeros(N,N,iend-istart+1);


for i = istart:iend
    tmp = imread(sprintf('Images/Syn1sec%d.tif',i));
    I(:,:,1+i-istart) = imresize(tmp,[N,N]);
end


%make correction
a=sum(sum(I(:,:,36)));
b=sum(sum(I(:,:,38)));
c=sum(sum(I(:,:,37)));
I(:,:,37)=I(:,:,37)*(a+b)/2/c;

FI = fftn(I);
if(false)
    I1 = real( ifftn(scalen(FI,sigma,[1,0,0])));
    I2 = real(ifftn(scalen(FI,sigma,[0,1,0])));
    I3 = real(ifftn(scalen(FI,sigma,[0,0,1])));
    f = sqrt(I1.^2+I2.^2+I3.^2);
else
    I11 = real(ifftn(scalen(FI,sigma,[2,0,0])));
    I22 = real(ifftn(scalen(FI,sigma,[0,2,0])));
    I33 = real(ifftn(scalen(FI,sigma,[0,0,2])));
    f = I11+I22+I33;
end
f = f/max(max(max(f)));

%display slices
%cd imStacks;
%OrthoSlicer3d(f);



%try region filling
%coordinates for point inside the vescicle
x = 102;
y = 120;
z = 35;

A = [];
th = 0.2;

for ix = -20:19
    for iy = -20:19
        for iz = -20:19
            if f(x+ix,y+iy,z+iz) <0.3
                A = [A,[102,120,35]];
            else A = [A,[x+ix,y+iy,z+iz]]; 
            end
        end
        
    end
end

v = diag(cov(A));
%v = v.*5;

[ex, ey, ez] = ellipsoid(x,y,z,v(1),v(2),v(3));

%figure;
%surf(ex, ey, ez);
%axis equal;



display slices
cd imStacks;
OrthoSlicer3d(f);
hold on; surf(ex, ey, ez); hold off;












%region growing attempt
%coordinates for point inside the vescicle
%x = 102;
%y = 121;
%z = 70;

%cd ..;
%cd regionGrowing;
%[P,J] = regionGrowing(I, [x,y,z]);
%cd ..;

%cd imStacks;
%OrthoSlicer3d(J);

