function [block] = DOWNSAMPLE(DU, factors)
    % �²���
if isinteger(factors)
    factors = [factors,factors];
end
[H, W,~] = size(DU);
block = zeros(H / factors(1), W / factors(2),class(DU));

for r = 1:size(block,1)
    for c = 1:size(block,2)
        block(r,c) = DU(factors(1)*r, factors(2)*c);
    end
end
end