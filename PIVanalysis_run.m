%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script for testing PIVanalysis and methods


%save ('0.6_40_105_630_10min_4.5V.mat', 'u_o', 'v_o')
%save ('0.6_20_105_630_10min_4.5V.mat', 'u_o', 'v_o')
%save ('0.6_60_105_630_10min_4.5V.mat', 'u_o', 'v_o')
%save ('0.6_80_110_660_10min_4.5V.mat', 'u_o', 'v_o')
%save ('01_20_115_690_10min_4.5V.mat', 'u_o', 'v_o')
%save ('01_40_120_720_10min_4.5V.mat', 'u_o', 'v_o')
%save ('01_60_120_720_10min_4.5V.mat', 'u_o', 'v_o')

beep
%%
clear all

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
clear, clc                                                                                                                  
close all


%make sure in the correct folder
PIVanalysis_test = PIVanalysis(); 
%PIVanalysis_test.Main();

PIVanalysis_test.allfunctions;

% % this will initialize PIVanalysis
%%sis_test with all properties and methods defined in
% % PIVanalysis with all default definitions. it will also run the object 
% % initialization method.


      
      
      
      
      
      