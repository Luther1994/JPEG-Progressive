function [code,size] = EncodeDC(DIFF, codelength,Htable)
%{
            编码DC系数
        
            In JPEG Encoder, DC coefficients are encoded as the DIFF between
            current DC DCT coefficient && PRED which is the quantized DC values
            from the most recently coded 8x8 block from the SAME component.
            The DIFF is obtained from

                                    DIFF = ZZ(0)-PRED

            At the beginning of the scan && at the beginning of each restart
            interval, the PRED is initialed to 0.The TWO'S COMPLEMENT difference
            magnitudes are grouped into 12 categories——SSSS for each of WHICH
            a Huffman code is created.

            Args
                DCcoef  Quantized current DC coefficient of 8x8 block from one
                src component.
                whichchannel  Specify which component the encoding DC coefficient
                               from, 1 for Y && 2/3 for U/V

            Notes
                Below is table of the categories SSSS && corresponding DIFF value
                    SSSS							  DIFF values
                      0                                    0
                      1                                  -1, 1
                      2                              -3, -2, 2, 3
                      3                              -7……-4, 4……7
                      4                             -15……-8, 8……15
                      5                            -31……-16, 16……31
                      6                            -63……-32, 32……63
                      7                           -127……-64, 64……127
                      8                          -255……-128, 128……255
                      9						     -511……-256, 256……511
                      10						-1023……-512, 512……1023
                      11					   -2047……-1024, 1024……2047
%}
[XHUFCO,XHUFSI] = SortHuffTbl(Htable);
code = XHUFCO(codelength+1);    % Categary从0开始的，Matlab里面要+1
Size = XHUFSI(codelength+1);
code = dec2bin(code,Size);          % 有一些码字最高位是0，十进制不能直接转为二进制，需要指定码长
if DIFF ~= 0
    if DIFF < 0
        DIFF = DIFF +  2 ^ codelength - 1;
    end
    % 负数用反码存储
    code = [code dec2bin(DIFF,codelength)];
    % 最终结果是Huffman码加上additional bits，即码长的Huffman码加上
    % DIFF的二进制补码
end
size = length(code);
code = bin2dec(code);
end