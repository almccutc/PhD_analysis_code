clear;
close all;
clc;

%load j20_1s_80p_5min
%load j20_1s_40p_5min
load j20_1s_80p_15min
%load j20_0.6s_20p_5min


[nr nc] = size(A_master);
N = 32; % Number of jets
num_regs = (N/8)/2;

% Powers of 2 used to calculate binary to decimal are stored in vector2
vector2 = [1, 2, 4, 8, 16, 32, 64, 128];
Asum = zeros(nr, num_regs*2);
threshold = 16;% The maximum number of jets that can be on out of every set of 16
steps = nr;

% Calculate the percentage of jets that are on a t a time BEFORE changes and
% put in the matrix A_master_percentages column 1
A_master_percentages = zeros(steps, 2); % Initialize matrix
A_master_sums1 = squeeze(sum(A_master, 2)); % Sum of A_master across each row
A_master_percentages(:, 1) = A_master_sums1/nc; % Divide by the number of rows to get percentage

% Works for row 1
placement = 1:2:num_regs * 2;
for row=1:1 % For the first row only
    for shiftreg = 1:num_regs % Increment over every set of 16 jets
        
        Atemp = A_master(row, (shiftreg - 1) * 16 + 1 : shiftreg * 16); % Create a length 16 vector
        count = NumOn(Atemp); % Use NumOn function to find out the number of jets on
        p = find(Atemp==1) + ((shiftreg * 16) - 16); % Make a vector with indices of jets that are on
        % add something to it to make p
        % increment from 1:nr, instead of 1:16
        % y = randsample(p, length(p)); % Create a vector of randomly chosen, non-repeating, values from p
        
        for k = 1:(count-threshold)
            if count > threshold
                
                % Use RowOff function to turn a jet off for however long it
                % was going to be on for, returns updated
                % A_master and the last row needed to be turned off
                
                [r_off, col_off, A_master] = RowOff(row, p( k ), A_master);
                len = r_off - row + 1;
                [A_master, col_out] = Search(row, col_off, len, A_master);
                count = count - 1; % Make the count one less because a jet was turned off
            end
        end
        
        count = 0; % Reset the count for the next Atemp vector
        Atemp_half1 = Atemp(1:8);
        Atemp_half2 = Atemp(9:16);
        Atemp_powers1 = Atemp_half1 .* vector2; %mul t iply by the correct power of 2 for the individual 16
        Atemp_powers2 = Atemp_half2 .* vector2;
        summm1 = sum(Atemp_powers1); % Add up all those numbers to get the decimal equivalent
        summm2 = sum(Atemp_powers2);
        Asum(row, shiftreg * 2 - 1) = summm1; % Put in a matrix Asum the decimal number
        Asum(row, shiftreg * 2) = summm2;
    end
end

% Works for middle rows
for row = 2:nr-1 % For the middle row only (excludes the first and last rows)
    for shiftreg = 1:num_regs % Increment over every set of 16 jets
        Atemp = A_master(row, (shiftreg-1) * 16+1: shiftreg * 16);
        % Create a length 16 vector
        count = NumOn(Atemp) ; % Use NumOn function to find out the number of jets on
        p = find(Atemp==1) + ((shiftreg * 16)-16); % Make a vector with indices of jets that are on
        for k=1:length (p)
            % Use Jet_preference function to see if the jet was previously turned off
            if count > threshold && Jet_preference(row, p(k), A_master) == 1
                
                % Use RowOff funct ion to turn the preferenced jet off for
                % however long it was going to be on for, returns updated
                % A_master and the last row needed to be turned off
                [r_off, col_off, A_master] = RowOff(row, p(k), A_master);
                len = r_off - row + 1;
                [A_master, col_out] = Search2(row, col_off, len, A_master);
                count = count - 1;
            end
        end
        
        % Check to see you need to force more jet(s) off (in order)
        % because no more preferenced ones exist and the count is
        % still greater than the threshold
        
        if count > 5
            p_new = find(Atemp==1) + ((shiftreg * 16) - 16); % Make a new vector with indices of jets that are on
            for c = 1:length(p_new)
                if count > threshold
                    % Use RowOff function to turn a jet off for however long it
                    % was going to be on for, returns updated A_master
                    % and the last row needed to be turned off
                    [r_off, col_off, A_master] = RowOff(row, p_new(c), A_master);
                    len = r_off - row + 1;
                    [A_master, col_out] = Search2(row, col_off, len, A_master);
                    count = count - 1;
                end
            end
        end
        
        count = 0; % Reset the count for the next Atemp vec tor
        Atemp_half1 = Atemp(1:8);
        Atemp_half2 = Atemp(9:16);
        Atemp_powers1 = Atemp_half1 .* vector2; % Multiply by the correct power of 2 for the individual 16
        Atemp_powers2 = Atemp_half2 .* vector2;
        summm1 = sum(Atemp_powers1); % Addd up all those numbers to get the decimal equivalent
        summm2 = sum(Atemp_powers2);
        Asum(row, shiftreg * 2 - 1) = summm1; % Put in a matrix Asum the decimal number
        Asum(row, shiftreg * 2) = summm2;
    end
end

% Works for last row

for row = nr:nr % For the last row only
    for shiftreg = 1:num_regs % Increment over every set of jets
        Atemp = A_master (row, (shiftreg-1) * 16 + 1:shiftreg * 16); % Create a length 16 vector
        count = NumOn(Atemp); % Use NumOn function to find out the number of jets on
        p = find(Atemp==1) + ((shiftreg * 16 ) - 16); % Make a vector with indices of jets that are on
        
        for k = 1:length(p)
            if count > threshold && Jet_preference(row, p(k), A_master) == 1
                % Use LastRowOff function to turn off a jet, returns updated A_master
                [A_master] = LastRowOff(row, p(k), A_master);
                count = count - 1;
            end
        end
        
        % check to see if you need to force more jet(s) off (in order)
        % because no more preferenced ones exist and the count is
        % still greater than the threshold
        if count > threshold
            p_new = find(Atemp==1) + ((shiftreg * 16) - 16) ; %make a new vector with indices of jets that are on
            for c = 1:length(p_new)
                if count > threshold % Use LastRowOff function to turn off a jet, returns updated A_master
                    [A_master] = LastRowOff(row, p(k), A_master);
                    count = count - 1;
                end
            end
        end
        
        count = 0; % Reset the count for the next Atemp vector
        Atemp_half1 = Atemp(1:8);
        Atemp_half2 = Atemp(9:16);
        Atemp_powers1 = Atemp_half1 .* vector2; % Multiply by the correct power of 2 for the individual 16
        Atemp_powers2 = Atemp_half2 .* vector2;
        summm1 = sum(Atemp_powers1); % Add up all those numbers to get the decimal equivalent
        summm2 = sum(Atemp_powers2);
        Asum(row, shiftreg * 2 - 1) = summm1; % Put in a matrix Asum the decimal number
        Asum(row, shiftreg * 2) = summm2;
    end
end

[nrr,ncc] = size(Asum);
Asum(nrr,:) = zeros(1,ncc);

% calculate the percentage of jets that are on at a time AFTER changes and
% put in the matrix A_master_percentages column 2
A_master_sums2 = sum(A_master, 2); % Sum of A_master across each row
% A_master_percentages(:, 2) = A_master_sums2 ./ nr; % Divide by the number of rows to get percentage

%save j20_1s_20p_10min_limit10_lastrow0 Asum
%save j20_1s_40p_5min_limit10_lastrow0 Asum
save j20_1s_80p_15min_limit10_lastrow0 Asum


% save hannah_1s_10p_5min_test_limit5_lastrow0 
% save hannah_1s_10p_5hr_limit5_lastrow0
% save hannah_1s_5p_26min_limit5_lastrow0

