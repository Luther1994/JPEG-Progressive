function [result]  = ReadNBytes(file,n,type)
if nargin==2
    if n == 0.5
        type = 'ubit4';
    elseif n == 2
        type = 'uint16';
    else
        type = 'uint8';
    end
end
if n == 0.5
    result = fread(file,1,type);
elseif n == 2
    result =swapbytes(uint16(fread(file,1,type)));
else
    result = fread(file,n,type);
end
end

