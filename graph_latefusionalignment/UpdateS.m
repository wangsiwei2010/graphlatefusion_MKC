function [Z, beta] = updateS(E, k, lambda)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
num = size(E, 1);
Z = zeros(num);
[SE, ind] = sort(E, 1);
Ek = SE(2:k+2, :);
Ek1 = SE(k+2, :);
Emin = Ek1 - Ek;
SZ = max(Emin./(sum(Emin)+eps),0);
for i=1:num
    Z(ind(2:k+2,i),i) = SZ(:,i);
end
beta = lambda*(sum(sum(Emin)))/num/2;
end

