function WriteNBytes(fid,Bits,precision)
if ischar(Bits)
    Bits = hex2dec(Bits);
end
if precision == 0.5
    fwrite(fid, Bits, 'ubit4', 0, 'b');  % 一次写入半个字节（4个bit）
elseif precision == 1
    fwrite(fid, Bits, 'uint8', 0, 'b');  % 一次写入一个字节
elseif precision == 2
    fwrite(fid, Bits, 'uint16', 0, 'b'); % 一次写入两个字节
end
end
