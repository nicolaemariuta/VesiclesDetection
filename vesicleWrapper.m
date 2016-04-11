clearvars;

N = 256;
r = 5;
sigma = 2;
T = 15;
alpha = .001;
beta = .005;
istart = 170;
iend = 239;
in = 40;
numError = 0;
snakeRuns = 170;
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
iborder = iborder+5;

drawnow;
x = 0;
y = 0;
fn = zeros(0);
In = zeros(0);
numErrorPrev = 0;

gridx = 1:size(I,1);
gridy = 1:size(I,2);
gridz = 1:size(I,3);

allEllipsoids = zeros(0);
vesicleClicks = zeros(0);


%gather clicks from the middle of vesicles
%click outside image to move to next slice or stop
for click = 1:1000
   
    imagesc(I(:,:,iborder)); colormap(gray);
    title(sprintf('Image %d',iborder+istart));
    [x,y] = ginput(1);
    z = iborder;
    
    if(x<0||y<0)
        if(iborder == (iend-istart-5))
            break;
        else
            iborder = iborder+1;
        end
    else
      
    found = false;
    %check if the vesicle is not already clicked    
    for index = 1 : size(vesicleClicks,1)
         point = vesicleClicks(index,:); 
         if ((x-point(1) <10) && (y-point(2) <10) && (z-point(3) <5))
           found = true;     
         end
    end
     
    if found == false
        vesicleClicks = [x y z; vesicleClicks];
    end
    
    end
    
end





%display slices
%figure;
%cd imStacks;
%OrthoSlicer3d(I);
%cd ..;

averageEnergy = zeros(0);

%alpha = .001;

nrVesicles = 0;

snakeRuns = 210;




for run = 1 : 1

statisticEnergy = zeros(0);

for index = 1 : size(vesicleClicks,1)
  
    point = vesicleClicks(index,:);
    x = point(1);
    y = point(2);
    z = point(3);
    dataEllipsoid = zeros(0);
    
    for in=1:3
  
    %init for horizontal slice
    if(in == 1)
     
        fn = f(:,:,z);
        In = I(:,:,z);
        
        t = linspace(0,2*pi,T);
        X = r*[cos(t);sin(t)]+[x;y]*ones(1,length(t));
      %  hold on; plot(x,y,'r+'); hold off;
      %  imagesc(fn);
      %  hold on; plot(X(1,:),X(2,:),'r-'); hold off
     %   title(sprintf('Image %d',z+istart));
      %  drawnow;
        
        Y=flipud(X);
        
    end
    
     %init for vertical slice
     if(in == 2)
      %  [x,y] = ginput(1);
        
        fn = f(:,round(x),:);
        In = I(:,round(x),:);
        
        fn = reshape(fn,[size(fn,1),size(fn,3)]);
        In = reshape(In,[size(In,1),size(In,3)]);
    
       
        
        t = linspace(0,2*pi,T);
        X = r*[cos(t);sin(t)]+[z;y]*ones(1,length(t));
     %   hold on; plot(z,y,'r+'); hold off;
     %   imagesc(fn);colormap(gray);
     %   hold on; plot(X(1,:),X(2,:),'r-'); hold off
     %   title(sprintf('Image %d',z+istart));
      %  drawnow;
        
        Y=flipud(X);
     end
 
    
      
     %init for second vertical slice
     if(in == 3)
    %  [x,y] = ginput(1);
        
        fn = f(round(y),:,:);
        In = I(round(y),:,:);
        
        fn = reshape(fn,[size(fn,2),size(fn,3)]);
        In = reshape(In,[size(In,2),size(In,3)]);
    
    
        
        t = linspace(0,2*pi,T);
        X = r*[cos(t);sin(t)]+[z;x]*ones(1,length(t));
    %    hold on; plot(x,z,'r+'); hold off;
    %    imagesc(fn);colormap(gray);
    %   hold on; plot(X(1,:),X(2,:),'r-'); hold off
    %    title(sprintf('Image %d',z+istart));
    %    drawnow;
        
        Y=flipud(X);
     end
     
     
     
        
   
        try
         %apply snake algorithm
        for i = 1:snakeRuns
            dY = Y(:,2:end)-Y(:,1:end-1);
            l = [0,cumsum(sqrt(sum(dY.^2,1)))];
            Y = [interp1(l,Y(1,:),linspace(min(l),max(l),max(l)-min(l)+1)); interp1(l,Y(2,:),linspace(min(l),max(l),max(l)-min(l)+1))];
            Y(:,end+1)=Y(:,1);
            Y = snake(Y',alpha,beta,fn,100,100)';
   %         imagesc(fn);
   %        hold on; plot(X(1,:),X(2,:),'r-'); hold off
    %        hold on; plot(Y(2,:),Y(1,:),'g.-'); hold off
  %         title(sprintf('Image %d',iborder+istart));
   %         drawnow;
        end
        catch err
            numError = numError +1;
        end
        
        %add curve points to data for building the ellipsoid
   
        Z = transpose(Y);
        
  if(in == 1)     
        Z = [Z(:,2),Z(:,1)];
        data1 = zeros(size(Z,1),1)+(z);
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
  
    
  %      imagesc(In);
   %     hold on; plot(Y(2,:),Y(1,:),'g.-'); hold off
   %     title(sprintf('Image %d',in+istart));
        
        
  
    end
    
    if(numError ==numErrorPrev)
        %generate ellipsoid, check if it is not already existent and add to
         %collection of ellipsoids
        [~,~,~,ex,ey,ez] = fitElipsoid(dataEllipsoid,2);
        E = [ex,ey,ez];
    
     %    sumE = sum(E);
   
      %   found = false;
   
      %   for check = 1:size(allEllipsoids,1)
      %      Ec = allEllipsoids(check,:);
      %       if sum(E)-sum(sumE) <20
      %          found = true;
      %       end
      %   end
   
       
  
     %       if found==false
              
        %       hold on; surf(ex,ey,ez); hold off;
         %      drawnow;
               
               
               %calculate total energy of current ellipsoid in all its
               %points and store the value into array
             try
               totalEnergy = 0;
               
               for i = 1:size(ex,1)
                   for j = 1:size(ex,2)
                   totalEnergy = totalEnergy + f(round(ex(i,j)),round(ey(i,j)),round(ez(i,j)));
                   end
               end
               

%                 for i = 1:size(dataEllipsoid)
%                     point = dataEllipsoid(i,:);
%                     totalEnergy = totalEnergy + f(round(point(1)),round(point(2)),round(point(3)));
%                 end
%  nrVesicle = nrVesicles+1;  
%                
                nrVesicle = nrVesicles+1;  
               statisticEnergy = [(totalEnergy); statisticEnergy ];
             catch
             end
               
            
               
    %       end
    end
            numErrorPrev = numError;
            
         
       
    
end
   break;
%calculate averga energy of all the ellipsoids to compare results between
%runs
%average = sum(statisticEnergy)/nrVesicles;
%averageEnergy = [ellipsoidDeviation average; averageEnergy];

 %ellipsoidDeviation =ellipsoidDeviation + 0.4;

            figure;
            hold on;plot(statisticEnergy); hold off;
            xlabel('ellipsoidIndex'), ylabel('energy');
            drawnow;
            
         
 
 
 
end
        
    


 %  cd imStacks;  
       
    %     hold on; plot(ex,ey,'b.-'); hold off
   

      
     %   hold on; surf(data1); hold off;
       % cd ..;
       
    %   figure;
% now plot the rotated ellipse
% sc = surf(x,y,z); shading interp; colormap copper
%h = surfl(ez, ey, ex); colormap copper
%title('actual ellipsoid represented by mu and Cov')
%axis equal

  
   
