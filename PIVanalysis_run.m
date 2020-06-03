%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script for testing PIVanalysis and methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc
delete all

% Load data. specific for PIVLAB mat files

load('F:\032420_AM_0.6_40_105_630_10min_4.5V_Correct.mat','u_original','v_original') %Change% 
%load('032420_AM_0.6_40_105_630_10min_4.5V_Correct.mat','u_original','v_original') %Change% 
%load('032420_AM_01_80_120_720_10min_4.5V_Correct.mat','u_original','v_original')
%load('4Vworkspace.mat','Utotal','Wtotal');
%load('6Vworkspace.mat','Utotal','Wtotal');

% u_o=permute(Utotal,[2 1 3]);
% v_o=permute(Wtotal,[2 1 3]);

 u_o = cell2mat(u_original); %No filter is done within PIVLAB
 v_o = cell2mat(v_original); %These values have NaNs in them
 
 % 
% % Rotate Matrix
 nlay = length(u_original);
 [r,c] = size(u_o);
% 
 u_o   = permute(reshape(u_o',[c,r/nlay,nlay]),[2,1,3]);
 v_o   = permute(reshape(v_o',[c,r/nlay,nlay]),[2,1,3]);

%save ('0.6_40_105_630_10min_4.5V.mat', 'u_o', 'v_o')
%save ('0.6_20_105_630_10min_4.5V.mat', 'u_o', 'v_o')
%save ('0.6_60_105_630_10min_4.5V.mat', 'u_o', 'v_o')
%save ('0.6_80_110_660_10min_4.5V.mat', 'u_o', 'v_o')
%save ('01_20_115_690_10min_4.5V.mat', 'u_o', 'v_o')
%save ('01_40_120_720_10min_4.5V.mat', 'u_o', 'v_o')
%save ('01_60_120_720_10min_4.5V.mat', 'u_o', 'v_o')
%save ('J_16_6Volts.mat', 'u_o', 'v_o')

beep
%% Bring in data, the is for a 3D double array, mine is 53 (height) by 79 (length) by 10,500 (in time)
clear all

%load ('J_16_4Volts.mat', 'u_o', 'v_o')
%load ('J_16_6Volts.mat', 'u_o', 'v_o')

%load ('0.6_40_105_630_10min_4.5V.mat', 'u_o', 'v_o')
%load ('0.6_20_105_630_10min_4.5V.mat', 'u_o', 'v_o')
%load ('0.6_60_105_630_10min_4.5V.mat', 'u_o', 'v_o')
%load ('0.6_80_110_660_10min_4.5V.mat', 'u_o', 'v_o')
load ('01_20_115_690_10min_4.5V.mat', 'u_o', 'v_o')
%load ('01_40_120_720_10min_4.5V.mat', 'u_o', 'v_o')
%load ('01_60_120_720_10min_4.5V.mat', 'u_o', 'v_o')
%load ('01_80_120_720_10min_4.5V.mat', 'u_o', 'v_o')

beep

%%
clc
% 

PIVanalysis_test=PIVanalysis(); 
%PIVanalysis_test.Main();

PIVanalysis_test.checkHistogram;

% % this will initialize PIVanalysis_test with all properties and methods defined in
% % PIVanalysis with all default definitions. it will also run the object 
% % initialization method.


      
      
      
      
      
      