function [MAC_MATRIX] = compute_mac(U1,U2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% mac

MAC_MATRIX=zeros(size(U1,2),size(U2,2));

for i=1:size(U1,2)
    for j=1:size(U2,2)
        MAC_MATRIX(i,j)=(dot((U1(:,i)),U2(:,j))^2)/dot(dot(U1(:,i),U1(:,i)),dot(U2(:,j),U2(:,j)));
    end
end

% Calcola semplicemente il MAC fra i modi FEM calcolati
% (eg U1 è la generica colonna del modo calcolato) e i modi ID identificati
% (eg U2 è la generica colonna del modo identificato)


