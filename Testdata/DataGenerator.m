clear
clc
close all

N = 600;
A = zeros(N, N);

for i = 1:N
    for j = 1:N
        if (sqrt( (N/2 - i)^2 + (N/2 - j)^2) < 3*N/8)
            A(i, j) = 1;            
        end
    end
end

dlmwrite('test.txt',A,'delimiter',' ','precision', 1)

    