clearvars;
N = 256;
r = 5;
sigma = 3;
T = 15;
alpha = .001;
beta = .01;
istart = 170;
iend = 239;
in = 10;
deviation = 0.008;
I = zeros(N,N,iend-istart+1);
for i = istart:iend
    tmp = imread(sprintf('Images/Syn1sec%d.tif',i));
    I(:,:,1+i-istart) = imresize(tmp,[N,N]);
end

figure(1);
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
iborder = iborder + 30;
figure(1);
imagesc(I(:,:,iborder)); colormap(gray);
drawnow;
[x,y] = ginput(1);
t = linspace(0,2*pi,T);
Y = flipud(r*[cos(t);sin(t)]+[x;y]*ones(1,length(t)));
firstTime = true;
%Data = cell(size(I,3)-iborder,1);
%ind = iborder+[0:size(Data,1)-1];

dataEllipsoid = zeros(0);

prevCurve = zeros(0);


%move up
for j=1:100                           %length(ind) %(size(I,3)-iborder)
    in = iborder+j-1;
    
    if(in > size(I,3))
        break;
    end
    
    fn = f(:,:,in);
    In = I(:,:,in);
    
    imagesc(fn);
    Yold = Y;
    if(firstTime)
        M = 200;
        firstTime=false;
    else
        M = 50;
    end
    for i = 1:M
        dY = Y(:,2:end)-Y(:,1:end-1);
        l = [0,cumsum(sqrt(sum(dY.^2,1)))];
        Y = [interp1(l,Y(1,:),linspace(min(l),max(l),max(l)-min(l)+1)); interp1(l,Y(2,:),linspace(min(l),max(l),max(l)-min(l)+1))];
        Y(:,end+1)=Y(:,1);
        Y = snake(Y',alpha,beta,fn,100,100)';
        imagesc(fn);
        hold on; plot(Yold(2,:),Yold(1,:),'r-'); hold off
        hold on; plot(Y(2,:),Y(1,:),'g.-'); hold off
        title(sprintf('Image %d',in+istart));
        drawnow;
    end
    Data{j} = Y;
    
    imagesc(In);
    hold on; plot(Y(2,:),Y(1,:),'g.-'); hold off
    title(sprintf('Image %d',in+istart));
    
    C = cov(Y');
    [E,V] = eig(C);
    M = mean(Y,2);
    t = linspace(0,2*pi,36);
    Z = [sqrt(2*V(1,1))*cos(t);sqrt(2*V(2,2))*sin(t)];
    Z = E'*Z+M*ones(1,length(t));
    hold on; plot(Z(2,:),Z(1,:),'b-','linewidth',3); hold off
    
        %create ellipsoid
    Z = transpose(Y);
    Z = [Z(:,2),Z(:,1)];
    data1 = zeros(size(Z,1),1)+(in);
    data1 =[Z,data1];
    
    
   %calculate difference between points of prev curve and next curve
   
   if(abs(prevCurve) ~= 0)
      diff = abs(sum(sum(prevCurve)) - sum(sum(Y)));
       
      if(diff/sum(sum(Y)) > deviation)
           break;
       end
   end

    dataEllipsoid = [data1;dataEllipsoid];
    
    prevCurve = Y;
end


t = linspace(0,2*pi,T);
Y = flipud(r*[cos(t);sin(t)]+[x;y]*ones(1,length(t)));

prevCurve = zeros(0);

%move down
for j=1:100                           %length(ind) %(size(I,3)-iborder)
    in = iborder-j-1;
    fn = f(:,:,in);
    In = I(:,:,in);
    
        if(in == 0)
            break;
        end
    
    
    
    imagesc(fn);
    Yold = Y;
    if(firstTime)
        M = 200;
        firstTime=false;
    else
        M = 50;
    end
    for i = 1:M
        dY = Y(:,2:end)-Y(:,1:end-1);
        l = [0,cumsum(sqrt(sum(dY.^2,1)))];
        Y = [interp1(l,Y(1,:),linspace(min(l),max(l),max(l)-min(l)+1)); interp1(l,Y(2,:),linspace(min(l),max(l),max(l)-min(l)+1))];
        Y(:,end+1)=Y(:,1);
        Y = snake(Y',alpha,beta,fn,100,100)';
        imagesc(fn);
        hold on; plot(Yold(2,:),Yold(1,:),'r-'); hold off
        hold on; plot(Y(2,:),Y(1,:),'g.-'); hold off
        title(sprintf('Image %d',in+istart));
        drawnow;
    end
    Data{j} = Y;
    
    imagesc(In);
    hold on; plot(Y(2,:),Y(1,:),'g.-'); hold off
    title(sprintf('Image %d',in+istart));
    
    C = cov(Y');
    [E,V] = eig(C);
    M = mean(Y,2);
    t = linspace(0,2*pi,36);
    Z = [sqrt(2*V(1,1))*cos(t);sqrt(2*V(2,2))*sin(t)];
    Z = E'*Z+M*ones(1,length(t));
    hold on; plot(Z(2,:),Z(1,:),'b-','linewidth',3); hold off
    
    %create ellipsoid
    Z = transpose(Y);
    Z = [Z(:,2),Z(:,1)];
    data1 = zeros(size(Z,1),1)+(in);
    data1 =[Z,data1];
    
   %calculate difference between points of prev curve and next curve
   
     if(abs(prevCurve) ~= 0)
      diff = abs(sum(prevCurve) - sum(Y));
       
       
       if(diff/sum(Y) > deviation)
           break;
       end
   end
    dataEllipsoid = [data1;dataEllipsoid];
    
    prevCurve = Y;
end








%plot curves
figure(2);  hold on; for i = 1:length(Data); Y = Data{i}; plot3(Y(2,:),Y(1,:),i*ones(size(Y,2)),'g.-'); end; hold off



 %generate ellipsoid
 [~, ~, ~, ex,ey,ez] = fitElipsoid(dataEllipsoid);
    
 %display slices
 figure;
 cd imStacks;
 OrthoSlicer3d(I);
 hold on; surf(ex,ey,ez); hold off;
     %   hold on; surf(data1); hold off;
       % cd ..;

