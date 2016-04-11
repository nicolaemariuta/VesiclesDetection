%read images and build block
N = 256;
istart = 170;
iend = 239;
I = zeros(N,N,iend-istart+1);


for i = istart:iend
    tmp = imread(sprintf('Images/Syn1sec%d.tif',i));
    I(:,:,1+i-istart) = imresize(tmp,[N,N]);
end


%display slices
cd imStacks;
OrthoSlicer3d(I);


while (true)  %Infinite loop to read clicked pixel
    [xVector,yVector,zVector] = ginput(1);  
    x = int64(xVector(1));  %coordinates of clicked pixel
    y = int64(yVector(1));
    z = int64(zVector(1));
    disp(sprintf('X=%u', y));  %display coordinates in console
    disp(sprintf('Y=%u', x));
    disp(sprintf('Z=%u', z));
    OrthoSlicer3d(I),title(sprintf('Clicked (X=%u,Y=%u,Z=%u)',y,x,z)); 
   
end