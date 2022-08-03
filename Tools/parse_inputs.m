function X = parse_inputs(varargin)

narginchk(1,1);  %ȷ����������ĸ���
X = varargin{1};
if (ndims(X) == 3)
    validateattributes(X,{'uint8','uint16','single','double'},{'real'}, ...
        mfilename,'RGB',1);
    if (size(X,3) ~= 3)
        error(message('images:rgb2yuv:invalidTruecolorImage'));
    end
    X = int16(X);
else
    error(message('images:rgb2yuv:invalidInputSize'))
end
end
