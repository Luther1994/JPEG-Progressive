function [res] =zigzag(src)
% 	Zig_Zag coding the frequency coefficient
% 
%   zigzag code ,transform coefficient matrix to 1d array with zigzag order
%   src£ºmatrix of coefficients 
if ~ismatrix(src)
    H, W, C = size(src);
    res = zeros(H * W, C);
    for i =1:C
        res(:, i) = Transform(src(:, :, i),false);
    end
else
    res = Transform(src,false);
end
end