% Original code developed by Evan Variano with Professor Edwin Cowen
% This function will be used as a part of the RunJets program.
% This code will generate random matrices that control the state of each
% pump. Each row in the matrix represents a tenth of a second so the on
% time for each pump will be represented by the number of rows that contain a
% ’1’ in each column.

function out = pump_update_independent (A, mean_on_time, mean_off_time, sigma_on_time, sigma_off_time)
if sigma_on_time > mean_on_time;
    beep;
    'error - sigma_on_time should be less than 3x mean_on_time';
    return;
end

% Update procedure: for each pump, if the next step is going to be a 9, then
% choose an on or off time and fill in the corresponding number of next
% steps with zero or 1

% Update the matrix
for i = 1:32 % 64 % go through each pump
    if (A(2, i) == 9) % If next step is undefined
        if (A(1, i ) == 1) % If this pump is on right now
            offlength = round(normrnd(mean_off_time(i), sigma_off_time(i))); % Choose an off time
            if offlength < 1;
                offlength = 1;
            end
            A(2:offlength + 1, i) = 0;
            
% Fill in this off time as the number of rows. The +1 will be
% eliminated later and is only there for convenience

        elseif(A(1, i) == 0) % If pump is off right now
            onlength = round(normrnd(mean_on_time(i), sigma_on_time(i))); % Choose an on time
            if onlength < 1;
                onlength = 1;
            end
            A(2:onlength+1, i) = 1;
            
% Fill in this on time as the number of rows. The +1 will be
% eliminated later and is only there for convenience

        end
    end
end

% Correct for the +1 in order to normalize the run times
B = A;
A(1:599, :) = B(2:600, : ) ;
A(600, :) = 9;
% Send this new matrix
out = A; 
