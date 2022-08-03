function [Q]= GetQuantizer(Q, quality)
% 	get_quantizer(Q,quality)
% 		In JPEG compress,quantizer should be controlled by the parameter quality
%
% 		Parameter
% 			Q : Standard Quantize Matrices in which
% 			QUANTIZER_LUMINANCE =  [[16, 11, 10, 16, 24, 40, 51, 61],
%                                   [12, 12, 14, 19, 26, 58, 60, 55],
%                                   [14, 13, 16, 24, 40, 57, 69, 56],
%                                   [14, 17, 22, 29, 51, 80, 87, 62],
%                                   [18, 22, 37, 56, 68, 109, 103, 77],
%                                   [24, 35, 55, 64, 81, 104, 113, 92],
%                                   [49, 64, 78, 87, 103, 121, 120, 121],
%                                   [72, 92, 95, 98, 112, 100, 103, 99]]
%
% 			QUANTIZER_CHROMINANCE =[[17, 18, 24, 47, 99, 99, 99, 99],
%                                   [18, 21, 26, 66, 99, 99, 99, 99],
%                                   [24, 26, 56, 99, 99, 99, 99, 99],
%                                   [47, 66, 99, 99, 99, 99, 99, 99],
%                                   [99, 99, 99, 99, 99, 99, 99, 99],
%                                   [99, 99, 99, 99, 99, 99, 99, 99],
%                                   [99, 99, 99, 99, 99, 99, 99, 99],
%                                   [99, 99, 99, 99, 99, 99, 99, 99]]
%
% 			quality: Parameter to define quantizer
Q = uint16(Q);                      % 预留内存，防止质量因子太大的时候数据溢出
assert( (0 < quality) && (quality < 100),...,
        'Quality must be GT 0 and LT 100.')
if quality <= 50
    Q = ceil(50 * Q / quality);
else
    Q = ceil((100 - quality) * Q / 50);
end
end