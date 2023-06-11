function [FREQ_MATRIX] = compute_freq2(U1,U2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
FREQ_MATRIX=zeros(size(U1,2),size(U2,2));

for i=1:size(U1,1)
    for j=1:size(U2,1)
        if i ~= j
%         FREQ_MATRIX(i,j)= abs((U1(i,1)-U2(j,1))/(max(U1(i,1),U2(j,1))));
        FREQ_MATRIX(i,j)= abs((U1(i,1)-U2(j,1)));
        end
    end
end

% Calcola la differenza tra due frequenze