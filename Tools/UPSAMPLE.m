function [block]=UPSAMPLE(DU,factors)
if isinteger(factors)
    factors = [factors,factors];
end
[H, W] = size(DU);
block = imresize(DU,[H*factors(1),W*factors(2)],'bicubic');
end