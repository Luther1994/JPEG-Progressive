function [code] = PrepareforMarker(code)
% Terminate the entropy-data segment
%
% This procedure terminates the entropy-data segment by
% 		a) padding a Huffman-coded segment with 1-bits to complete the final byte(&& if needed stuffing a 0 byte);
% 		b) invoking the procedure 'Flush' to terminate an arithmetic entropy-coded segment.
% 	In DCT based Huffman coded JFIF , if last bits are less than 8 to be converted to a byte, several 0s could be
% 	add in LSB to build a byte to fill the file.
% �����������ȫ��д��֮����󼸸�bit����д��һ���ֽڣ�����LSB�����0ֱ
% ���չ�8��bits��Ȼ����д���ļ��С�

assert( 0 < length(code) && length(code)< 8)
while length(code) < 8
    code = [code  '0'];
end

end