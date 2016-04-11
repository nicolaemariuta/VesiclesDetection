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

%read images
I = zeros(N,N,iend-istart+1);

for i = istart:iend
    tmp = imread(sprintf('Images/Syn1sec%d.tif',i));
    I(:,:,1+i-istart) = imresize(tmp,[N,N]);
end
figure(1)
imshow3D(I);
title('1 x 2')
figure(2)
J = permute(I,[3,1,2]);
imshow3D(J);
title('3 x 1')
figure(3)
K = permute(I,[2,3,1]);
imshow3D(K);
title('2 x 3')

% Prepare the Image term.  The function scale is a convolution with a
% Gaussian and may be found on my homepage.
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



% Setup the figure
iborder = ceil(sigma);
imagesc(I(:,:,iborder)); colormap(gray);
drawnow;
for in=iborder:(size(I,3)-iborder)
    fn = f(:,:,in);
    In = I(:,:,in);
    
    for j = 1:3
        figure;
        
        [x,y] = ginput(1);
        t = linspace(0,2*pi,T);
        X = r*[cos(t);sin(t)]+[x;y]*ones(1,length(t));
        hold on; plot(x,y,'r+'); hold off;
        imagesc(fn);
        hold on; plot(X(1,:),X(2,:),'r-'); hold off
        title(sprintf('Image %d',in+istart));
        drawnow;
        
        Y=flipud(X);
        for i = 1:100
            Y = snake(Y',alpha,beta,fn,100,100)';
            imagesc(fn);
            hold on; plot(X(1,:),X(2,:),'r-'); hold off
            hold on; plot(Y(2,:),Y(1,:),'go-'); hold off
            title(sprintf('Image %d',in+istart));
            drawnow;
        end
        
        imagesc(In);colormap(gray);
        hold on; plot(Y(2,:),Y(1,:),'go-'); hold off
        title(sprintf('Image %d',in+istart));
        
    end
    
end



