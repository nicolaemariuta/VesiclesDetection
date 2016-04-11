I = dicomread('volumedata.tif');
%addpath('D:\master\project\imStacks');
%Image = squeeze('volumedata.tif');
imshow3D(I);
figure; for i=1:1065;I=imread('volumedata.tif',i); imagesc(I); colormap(gray); axis equal tight; title(sprintf('image %d',i)); drawnow; end