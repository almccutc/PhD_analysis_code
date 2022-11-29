clc;
clear all;
workingDir = tempname;
mkdir(workingDir)
mkdir(workingDir,'images')

shuttleVideo = VideoReader('/Users/almccutc/Desktop/Hero8Test_Laser5_autosettings.mov');

ii = 1;
%

while hasFrame(shuttleVideo)
   img = readFrame(shuttleVideo);
   filename = [sprintf('%05d',ii) '.tif'];
   fullname = fullfile(workingDir,'images',filename);
   imwrite(img,filename)    % Write out to a JPEG file (img1.jpg, img2.jpg, etc.)
   ii = ii+1;
end
%%
i = imread('/Users/almccutc/Documents/MATLAB/PIVanalysis/00010.tif');

I = rgb2gray(i);
imshow(I)

J = imadjust(I,[0 0.7],[]);

imshowpair(I,J,'montage') %shows images side by side

