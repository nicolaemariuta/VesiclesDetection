N = 256;
r = 5;
sigma = 3;
T = 15;
alpha = .001;
beta = .01;
istart = 170;
iend = 239;
in = 10;
I = zeros(N,N,iend-istart+1);
for i = istart:iend
    tmp = imread(sprintf('Images/Syn1sec%d.tif',i));
    I(:,:,1+i-istart) = imresize(tmp,[N,N]);
end

I = correctSlice(I,37);


for i=1:size(I,3)
    imagesc(I(:,:,i)); colormap(gray);
    drawnow;
end

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

iborder = ceil(sigma);
iborder = iborder+30;
imagesc(I(:,:,iborder)); colormap(gray);
drawnow;

for in=iborder:(iborder)

    
    for j = 1:1
        [x,y] = ginput(1);
        
        fn = f(round(y),:,:);
        In = I(round(y),:,:);
        
    %    imagesc(fn);colormap(gray);
   %     figure;
        
        fn = reshape(fn,[size(fn,2),size(fn,3)]);
        In = reshape(In,[size(In,2),size(In,3)]);
    
        figure;
        
        t = linspace(0,2*pi,T);
        X = r*[cos(t);sin(t)]+[iborder;x]*ones(1,length(t));
        hold on; plot(iborder,y,'r+'); hold off;
        imagesc(fn);colormap(gray);
        hold on; plot(X(1,:),X(2,:),'r-'); hold off
        title(sprintf('Image %d',in+istart));
        drawnow;
        
        Y=flipud(X);
        
        figure;
        %horizontal surface
        for i = 1:200
            dY = Y(:,2:end)-Y(:,1:end-1);
            l = [0,cumsum(sqrt(sum(dY.^2,1)))];
            Y = [interp1(l,Y(1,:),linspace(min(l),max(l),max(l)-min(l)+1)); interp1(l,Y(2,:),linspace(min(l),max(l),max(l)-min(l)+1))];
            Y(:,end+1)=Y(:,1);
            Y = snake(Y',alpha,beta,fn,100,100)';
            imagesc(fn);colormap(gray);
            hold on; plot(X(1,:),X(2,:),'r-'); hold off
            hold on; plot(Y(2,:),Y(1,:),'g.-'); hold off
            title(sprintf('Image %d',in+istart));
            drawnow;
        end
        Z = transpose(Y);
     %  Z = [Z(:,2),Z(:,1)];
        data1 = zeros(size(Z,1),1)+(round(x));
        data1 =[data1,Z];
        
        
        imagesc(In);
        hold on; plot(Y(2,:),Y(1,:),'g.-'); hold off
        title(sprintf('Image %d',in+istart));
        
        
   
        
        
        %generate ellipsoid
        [~,~,~,ex,ey,ez] = fitElipsoid(data1);
       
        %display slices
        figure;
        cd imStacks;
        OrthoSlicer3d(I);
        hold on; surf(ex,ey,ez); hold off;
       % cd ..;
       
       figure;
% now plot the rotated ellipse
% sc = surf(x,y,z); shading interp; colormap copper
%h = surfl(ez, ey, ex); colormap copper
title('actual ellipsoid represented by mu and Cov')
axis equal

     
    end
   
end