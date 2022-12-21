clear;
clc;
close all

%February 2, 2020j

%create identity matrix to individually check pumps

N = 32;
A_master = zeros(N+1,N);
A_master(1:N,1:N) = eye(N);

[nr nc] = size(A_master);
% N = 256; % Number of jets
num_regs = (N/8)/2;

% Powers of 2 used to calculate binary to decimal are stored in vector2
vector2 = [1, 2, 4, 8, 16, 32, 64, 128];
Asum = zeros(nr, num_regs*2);
steps = nr;

A_master_percentages = zeros(steps, 2); % Initialize matrix
A_master_sums1 = squeeze(sum(A_master, 2)); % Sum of A_master across each row
A_master_percentages(:, 1) = A_master_sums1/nc; % Divide by the number of rows to get percentage

for row=1:N % For the first row only
    for shiftreg = 1:num_regs % Increment over every set of 16 jets
        
        Atemp = A_master(row, (shiftreg - 1) * 16 + 1 : shiftreg * 16); % Create a length 16 vector

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

save identity_test32 Asum