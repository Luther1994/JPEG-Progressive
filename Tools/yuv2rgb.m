function RGB = yuv2rgb(varargin)

%RGB2YUV Convert RGB color values to YUV color space.

yuv = parse_inputs(varargin{:});

T = [1 -0.00093 1.401687;...
    1 -0.3437 -0.71417; ...
    1 1.77216 0.00099];
offset = [0; 0; 0];

% Initialize output
RGB = zeros(size(yuv), 'like', yuv);

for p = 1:3
    RGB(:,:,p) = imlincomb(T(p,1),yuv(:,:,1),T(p,2),yuv(:,:,2)-128, ...
        T(p,3),yuv(:,:,3)-128,offset(p));
end
RGB = uint8(RGB);
end

