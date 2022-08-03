function YUV = rgb2yuv(varargin)

%RGB2YUV Convert RGB color values to YUV color space.

rgb = parse_inputs(varargin{:});

T = [ ...
    0.299 0.587 0.114;...
    -0.169 -0.331 0.5; ...
    0.5 -0.419 -0.081];
offset = [0; 128; 128];

% Initialize output
YUV = zeros(size(rgb), 'like', rgb);

for p = 1:3
    YUV(:,:,p) = imlincomb(T(p,1),rgb(:,:,1),T(p,2),rgb(:,:,2), ...
        T(p,3),rgb(:,:,3),offset(p));
end
YUV = uint8(YUV);
end

%%%
%Parse Inputs
%%%
