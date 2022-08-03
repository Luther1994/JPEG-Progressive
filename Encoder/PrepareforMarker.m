function [code] = PrepareforMarker(code)
% Terminate the entropy-data segment
%
% This procedure terminates the entropy-data segment by
% 		a) padding a Huffman-coded segment with 1-bits to complete the final byte(&& if needed stuffing a 0 byte);
% 		b) invoking the procedure 'Flush' to terminate an arithmetic entropy-coded segment.
% 	In DCT based Huffman coded JFIF , if last bits are less than 8 to be converted to a byte, several 0s could be
% 	add in LSB to build a byte to fill the file.
% 如果所有数据全部写入之后最后几个bit不能写满一个字节，则在LSB端填充0直
% 到凑够8个bits，然后将其写入文件中。

assert( 0 < length(code) && length(code)< 8)
while length(code) < 8
    code = [code  '0'];
end

end