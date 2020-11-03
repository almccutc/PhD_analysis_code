clc;
clear all;
workingDir = tempname;
mkdir(workingDir)
mkdir(workingDir,'images')


%Trial 2 calib
%shuttleVideo = VideoReader('slomo_1584120680_17-31-20.mov');

%Trial 3
%shuttleVideo = VideoReader('slomo_1584634793_16-19-53.mov');
%shuttleVideo = VideoReader('slomo_1584634967_16-22-47.mov');
%shuttleVideo = VideoReader('slomo_1584635130_16-25-30.mov');
%shuttleVideo = VideoReader('slomo_1584635367_16-29-27.mov');
shuttleVideo = VideoReader('032420_AM_01_80_120_720_10min_4.5V_0010.mov');

% 11 cm FOV

% 10 V 5% jets on
%shuttleVideo = VideoReader('slomo_1584047722_21-15-22.mov');
%shuttleVideo = VideoReader('slomo_1584047863_21-17-43.mov');
%shuttleVideo = VideoReader('slomo_1584048004_21-20-04.mov');

% 11 V 5% jets on
%shuttleVideo = VideoReader('slomo_1584048150_21-22-30.mov');
%shuttleVideo = VideoReader('slomo_1584048292_21-24-52.mov');
%shuttleVideo = VideoReader('slomo_1584048440_21-27-20.mov');

% 12 V 5% jets on
%shuttleVideo = VideoReader('');
%shuttleVideo = VideoReader('');
%shuttleVideo = VideoReader('');

% Calibration
%shuttleVideo = VideoReader('');

% 20 cm FOV

% 10 V 10% jets on
%shuttleVideo = VideoReader('slomo_1583528920_21-08-40.mov');
%shuttleVideo = VideoReader('slomo_1583529221_21-13-41.mov');
%shuttleVideo = VideoReader('slomo_1583529502_21-18-22.mov');

% 10 V 5% jets on
%shuttleVideo = VideoReader('slomo_1583952163_18-42-43.mov');
%shuttleVideo = VideoReader('slomo_1583952372_18-46-12.mov');
%shuttleVideo = VideoReader('slomo_1583952590_18-49-50.mov');

% 11 V
%shuttleVideo = VideoReader('slomo_1583529812_21-23-32.mov');
%shuttleVideo = VideoReader('slomo_1583530099_21-28-19.mov');
%shuttleVideo = VideoReader('slomo_1583530378_21-32-58.mov');

% 12 V
%shuttleVideo = VideoReader('slomo_1583530759_21-39-19.mov');
%shuttleVideo = VideoReader('slomo_1583531038_21-43-58.mov');
%shuttleVideo = VideoReader('slomo_1583531314_21-48-34.mov');

% Calibration
%shuttleVideo = VideoReader('calib_slomo_1583765584_14-53-04.mov');
%shuttleVideo = VideoReader('slomo_1583518450_18-14-10.mov');

%ii = 1;
%ii = 1250;
%ii =2499;
%ii = 3748;
%ii = 4997;
%ii = 6246;
%ii = 7495;
%ii = 8744;
%ii = 9993;
%ii = 11242;


while hasFrame(shuttleVideo)
   img = readFrame(shuttleVideo);
   filename = [sprintf('%05d',ii) '.tif'];
   fullname = fullfile(workingDir,'images',filename);
   imwrite(img,filename)    % Write out to a JPEG file (img1.jpg, img2.jpg, etc.)
   ii = ii+1;
end

%search for the long fullname variable on the computer to find the images


