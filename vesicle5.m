N = 256;
r = 5;
sigma = 3;
T = 15;
alpha = .001;
beta = .01;
istart = 170;
iend = 239;
in = 40;
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
iborder = iborder+33;
imagesc(I(:,:,iborder)); colormap(gray);
drawnow;
x = 0;
y = 0;
fn = zeros(0);
In = zeros(0);
dataEllipsoid = zeros(0);

for in=1:3
  
    %init for horizontal slice
    if(in == 1)
        [x,y] = ginput(1);
        fn = f(:,:,iborder);
        In = I(:,:,iborder);
        
        t = linspace(0,2*pi,T);
        X = r*[cos(t);sin(t)]+[x;y]*ones(1,length(t));
        hold on; plot(x,y,'r+'); hold off;
        imagesc(fn);
        hold on; plot(X(1,:),X(2,:),'r-'); hold off
        title(sprintf('Image %d',iborder+istart));
        drawnow;
        
        Y=flipud(X);
        
    end
    
     %init for vertical slice
     if(in == 2)
      %  [x,y] = ginput(1);
        
        fn = f(:,round(x),:);
        In = I(:,round(x),:);
        
        fn = reshape(fn,[size(fn,1),size(fn,3)]);
        In = reshape(In,[size(In,1),size(In,3)]);
    
        figure;
        
        t = linspace(0,2*pi,T);
        X = r*[cos(t);sin(t)]+[iborder;y]*ones(1,length(t));
        hold on; plot(iborder,y,'r+'); hold off;
        imagesc(fn);colormap(gray);
        hold on; plot(X(1,:),X(2,:),'r-'); hold off
        title(sprintf('Image %d',iborder+istart));
        drawnow;
        
        Y=flipud(X);
     end
 
    
      
     %init for second vertical slice
     if(in == 3)
    %  [x,y] = ginput(1);
        
        fn = f(round(y),:,:);
        In = I(round(y),:,:);
        
        fn = reshape(fn,[size(fn,2),size(fn,3)]);
        In = reshape(In,[size(In,2),size(In,3)]);
    
        figure;
        
        t = linspace(0,2*pi,T);
        X = r*[cos(t);sin(t)]+[iborder;x]*ones(1,length(t));
        hold on; plot(x,iborder,'r+'); hold off;
        imagesc(fn);colormap(gray);
        hold on; plot(X(1,:),X(2,:),'r-'); hold off
        title(sprintf('Image %d',iborder+istart));
        drawnow;
        
        Y=flipud(X);
     end
     
     
     
        
   
        
       %apply snake algorithm
        for i = 1:200
            dY = Y(:,2:end)-Y(:,1:end-1);
            l = [0,cumsum(sqrt(sum(dY.^2,1)))];
            Y = [interp1(l,Y(1,:),linspace(min(l),max(l),max(l)-min(l)+1)); interp1(l,Y(2,:),linspace(min(l),max(l),max(l)-min(l)+1))];
            Y(:,end+1)=Y(:,1);
            Y = snake(Y',alpha,beta,fn,100,100)';
            imagesc(fn);
            hold on; plot(X(1,:),X(2,:),'r-'); hold off
            hold on; plot(Y(2,:),Y(1,:),'g.-'); hold off
            title(sprintf('Image %d',iborder+istart));
            drawnow;
        end
        
        %add curve points to data for building the ellipsoid
        Z = transpose(Y);
        
  if(in == 1)     
        Z = [Z(:,2),Z(:,1)];
        data1 = zeros(size(Z,1),1)+(iborder);
        data1 =[Z,data1];
        
        dataEllipsoid = [data1;dataEllipsoid];
  end
  
  
    if(in == 2)     
      
        data2 = zeros(size(Z,1),1)+(round(x));
        data2 =[data2,Z(:,1),Z(:,2)];
        
         dataEllipsoid = [data2;dataEllipsoid];    
        
    end
  
      if(in == 3)     
     %   Z = [Z(:,2),Z(:,1)];
        data3 = zeros(size(Z,1),1)+(round(y));
        data3 =[Z(:,1),data3,Z(:,2)];

        
        dataEllipsoid = [data3;dataEllipsoid];    
      end
  
    
        imagesc(In);
        hold on; plot(Y(2,:),Y(1,:),'g.-'); hold off
        title(sprintf('Image %d',in+istart));
        
        
  
end
        
        
        %generate ellipsoid
        [~,~,~,ex,ey,ez] = fitElipsoid(dataEllipsoid);
       
         hold on; plot(ex,ey,'b.-'); hold off
        
        %display slices
        figure;
        cd imStacks;
        OrthoSlicer3d(I);
        hold on; surf(ex,ey,ez); hold off;
     %   hold on; surf(data1); hold off;
       % cd ..;
       
       figure;
% now plot the rotated ellipse
% sc = surf(x,y,z); shading interp; colormap copper
h = surfl(ez, ey, ex); colormap copper
title('actual ellipsoid represented by mu and Cov')
axis equal




%read image energy in the ellipsoid's points


ImEnergy = zeros(size(I))+1;
ImVesicle = zeros(size(I))+255;

for p = 1 : size(ez)
    ImEnergy(round(ex(p)),round(ey(p)),round(ez(p))) = f(round(ex(p)),round(ey(p)),round(ez(p))) ;
    ImVesicle(round(ex(p)),round(ey(p)),round(ez(p)))  = I(round(ex(p)),round(ey(p)),round(ez(p))) ;
end

figure;
OrthoSlicer3d(ImEnergy);
figure;
OrthoSlicer3d(ImVesicle);
   
%use interpn instead of round and 