clear
close all
% saved should then be inserted into Limit_ jets_replacements.m

N = 32; % Number of jets; edit pump_update_independent to reflect correct N inputs
mean_on_time = 10 * ones(1, N); %10 for 1 second
mean_off_time = 3 * ones(1, N);  %10/40-20%; 10/15-40% 
sigma_on_time = 3 * ones(1, N); %set as 1/3 of mean_on_time
sigma_off_time = 1 * ones(1, N); %set as 1/3 of mean_off_time

% Set up the matrix of pump states initial conditions
% This matrix holds the values for ri ght now and the near future
% A value of 9 means the pump hasnt chosen its on or off time yet

A12 = 9 + zeros(600, N); % 256 columns (for the pumps) and 600 rows for 60 seconds of buffersize
% Initial state of each pump is off

A12(1, :) = zeros(1, N);
for i = 1:3000 % 30000 Allows randomization to sort it all out
    A12 = pump_update_independent(A12, mean_on_time, mean_off_time, sigma_on_time, sigma_off_time);
end

save ic12 A12
load ic12

% Store actual matrix of on/off states for all time
time = 9000; %6000 should correspond to 10 minutes
A_master = zeros(time,N);
for kk = 1:time %180,000 = 5 hrs at 0.1 s update %24000 = 40 minutes at 0.1 s update %48000 = 80 minutes @ 0.1 s update
    A12 = pump_update_independent(A12, mean_on_time, mean_off_time, sigma_on_time, sigma_off_time);
    a = logical(A12(1, :));
    A_master(kk, :) = a;
end

 save j20_1s_80p_15min A_master kk %CHANGE in limit jet replacement twice, send_matrix to arduino once,  



