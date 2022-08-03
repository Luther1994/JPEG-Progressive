function [res]=inverse_zigzag(src)
%{
    Inverse zigzag,transform zigzaged coefficients to coefficient matrix
    src: 1d array of coefficients ordered as zigzag from coefficient matrix
%}

res = ones(8,8);
if ~ismatrix(src)
    C = size(src,2);
    for i =1:C
        res(:, :, i) = Transform(src(:, i), true);
    end
else
    res = Transform(src,true);
end
end