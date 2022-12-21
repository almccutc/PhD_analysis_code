 clear all
close all
clc

% BEFORE CONNECTINGj
% 1) Make sure all boards are off; switch is in middle position
% 2) Turn on board power supply - blue lights only should turn on
% 3) Run zeros on matlab by uncommenting Asum = [0,0,0,...] and commenting
% out loader and Asum(nr,:) = zeros;
%This is the code developed in Matlab to control an array of 256jets
%via the arduino with a stored matrix of the random jet states
%in order to zero all jets,can simply replace Asumas:
%Asum = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

%Jet states are supposed to update every 0.1s

%
%52 - left
%blue 53, middle 
%51 - right,  clock/latch/data = blue/red/grey

%for single jet tests, set on time%
load 1jettest_1secon_5secoff
Asum=onejet__6secon_5_4secoff_36secdelay;
%Asum=onejets_1secon_5secoff;
%Asum=onejet_1_4secon_4_6secoff;

% load identity_test32
%load j20_1s_20p_10min_limit10_lastrow0 %0.1089 t_avg with 0.095 pause
%load j20_1s_40p_5min_limit10_lastrow0 %0.1075 t_avg with 0.095 pause
%load j20_1s_60p_5min_limit10_lastrow0 %0.1056 t_avg with 0.095 pause
%load j20_1s_80p_5min_limit10_lastrow0 %0.1056 t_avg with 0.095 pause

%load j20_1_6s_20p_5min_limit10_lastrow0 %0.1074 t_avg with 0.095 pause
%load j20_1_6s_400p_5min_limit10_lastrow0 %0.1074 t_avg with 0.095 pause
%load j20_1_6s_60p_5min_limit10_lastrow0 %0.1075 t_avg with 0.095 pause
%load j20_1_6s_80p_5min_limit10_lastrow0 %0.1056 t_avg with 0.095 pause

%load j20_1_6s_20p_15min_limit10_lastrow0 %0.1075 t_avg with 0.095 pause
%load j20_1_6s_40p_15min_limit10_lastrow0 %0.1070 t_avg with 0.095 pause
%load j20_1_6s_60p_15min_limit10_lastrow0 %0.1076 t_avg with 0.095 pause
%load j20_1_6s_80p_15min_limit10_lastrow0 %0.1053 t_avg with 0.095 pause

%load j20_1s_20p_15min_limit10_lastrow0 %0.1076 t_avg with 0.095 pause
%load j20_1s_40p_15min_limit10_lastrow0 %0.1075 t_avg with 0.095 pause
%load j20_1s_60p_15min_limit10_lastrow0 %0.1075 t_avg with 0.095 pause
%load j20_1s_80p_15min_limit10_lastrow0 %0.1075 t_avg with 0.095 pause



%figure(1);
%imagesc(Asum);
%colorbar;

delete(instrfindall);


clear s
fclose('all')
% s=serial('/dev/tty.usbmodem411','BaudRate',115200)
s=serial('COM4','BaudRate',115200)   %COM6/COM8, 8 for board 9
fopen(s)

Azero=Asum;
[nr nc]=size(Azero);
kko=nr; %nr change to 1000 for test

Azero(kko,:) = zeros(1, nc);
timestep=0;
tstart=tic

while timestep < kko

rx = fgets(s);
timestep = timestep+1;

fwrite(s,Azero(timestep,1))
fwrite(s,Azero(timestep,2))
fwrite(s,Azero(timestep,3))
fwrite(s,Azero(timestep,4))
% fwrite(s,Azero(timestep,5))
% fwrite(s,Azero(timestep,6))
% fwrite(s,Azero(timestep,7))
% fwrite(s,Azero(timestep,8))
% fwrite(s,Azero(timestep,9));
% fwrite(s,Azero(timestep,10));
% fwrite(s,Azero(timestep,11));
% fwrite(s,Azero(timestep,12));
% fwrite(s,Azero(timestep,13));
% fwrite(s,Azero(timestep,14));
% fwrite(s,Azero(timestep,15));
% fwrite(s,Azero(timestep,16));
% fwrite(s,Azero(timestep,17));
% fwrite(s,Azero(timestep,18));
% fwrite(s,Azero(timestep,19));
% fwrite(s,Azero(timestep,20));
% fwrite(s,Azero(timestep,21));
% fwrite(s,Azero(timestep,22));
% fwrite(s,Azero(timestep,23));
% fwrite(s,Azero(timestep,24));
% fwrite(s,Azero(timestep,25));
% fwrite(s,Azero(timestep,26));
% fwrite(s,Azero(timestep,27));
% fwrite(s,Azero(timestep,28));
% fwrite(s,Azero(timestep,29));
% fwrite(s,Azero(timestep,30));
% fwrite(s,Azero(timestep,31));
% fwrite(s,Azero(timestep,32));

 pause(0.095); %adjusted to meet Fs=0.1Hz 
%0.03; Junior 0.04, Aubrey 0.06, Hannah T1: 0.06 T2: 0.06 T3: 0.06
% pause;
% % For single jet test
% pause % spacebar for test

end  

tstop=toc(tstart)
t_avg=tstop./kko
timestep=0;

fclose(s);
delete(instrfindall);