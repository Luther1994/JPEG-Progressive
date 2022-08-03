function  EncodeRestartInterval()
% 		"""
%
% 		Returns
%
% 		"""
ResetEncoder()

if CHANNELS == 3
    src = rgb2yuv(src);
    Y = src(:,:,1);
    U = src(:,:,2);
    V = src(:,:,3);
    
    Y = CompletePartialMCU(Y,SampleFactors(:,1));
    
    factors = [Hmax/SampleFactors(1,2),Wmax/SampleFactors(2,2)];
    U = DOWNSAMPLE(CompletePartialMCU(U,SampleFactors(:,2)), factors);
    
    factors = [Hmax/SampleFactors(1,3),Wmax/SampleFactors(2,3)];
    V = DOWNSAMPLE(CompletePartialMCU(V,SampleFactors(:,3)), factors);
    
    MCUcountsH = size(U,1) / 8;
    MCUcountsW = size(U,2) / 8;
else
    CHANNELS = 1;
    Y = src(:, :, 1);
    MCUcountsH = size(Y,1) /8;
    MCUcountsW = size(Y,2) /8;
end
MCUcounts = MCUcountsW * MCUcountsH;
DUcountsY = Hmax * Wmax;

    function [temp] = getMCU(idx)
        temp = zeros(DUcountsY+2,64);
        col = mod(idx , MCUcountsW);
        row = fix(double(idx) / MCUcountsW);
        microblocksY = Y(row * 8 * Hmax + 1 : (row + 1) * 8 * Hmax,...,
            col * 8 * Wmax + 1:(col + 1) * 8 * Wmax) - 2 ^ (precision - 1);
        for index = 0:DUcountsY-1
            r = mod(index,  Wmax);
            c = fix(double(index) / Wmax);
            coef = round(dct2(microblocksY(r* 8+1:(r+1) * 8,...,
                c * 8+1:(c+1)* 8)) ./ qlu);
            coef = zigzag(coef);
            temp(index+1,:) = coef;
        end
        microblock = U(row * 8+1:(row + 1) * 8, ...,
            col * 8+1:(col + 1) * 8) - 2 ^ (precision - 1);
        coef = round(dct2(microblock) ./ qchr);
        coef = zigzag(coef);
        temp(5,:) = coef;
        microblock = V(row * 8+1:(row + 1) * 8, ...,
            col * 8+1:(col + 1) * 8) - 2 ^ (precision - 1);
        coef = round(dct2(microblock) ./ qchr);
        coef = zigzag(coef);
        temp(6,:) = coef;
    end
code = '';
MCUidx = 0;
while MCUidx <= MCUcounts
    MCU = getMCU(MCUidx);
    for i = 1:length(MCU)
        if i < DUcountsY
            code = [code,encodeDU(MCU(i), 1)];
        else
            code = [code,encodeDU(MCU(i), i - DUcountsY + 1)];
        end
        while length(code) >= 8
            byte = bin2dec(code);
            write(byte,1)
            if byte == 255
                write('00',1)
            end
            code = code(8:end);
        end
    end
    MCUidx = MCU_idx + 1;
end
if ~isempty(code)
    PrepareforMarker(code)
end
end