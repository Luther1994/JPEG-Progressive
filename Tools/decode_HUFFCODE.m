function[HuffCodes]= decode_HUFFCODE(HuffSize)
% 	"""Generate the Huffman code table from HUFFSIZE table
%
% 	Args:
% 		HUFFSIZE: Table containing a code for each size in HUFFSIZE
%
% 	Returns:
% 		HUFFCODES：Huffman code table which contain the huffman code words with code size in HS
%
% 	Notes:
% 		The procedure for generating HUFFCODE from HUFFSIZE is listed below:
% 		1、Code start with 0 whatever the minimum code size of first code is.
% 		   eg. if HUFFSIZE[0] = 3, HUFFCODE[0] = 000;if HUFFSIZE[0] = 5,HUFFCODE[0]=00000;
% 		   and so on.
% 		2、If next code size is equal to current code size,CODE_next = CODE_current + 1,
% 		   00->01、0000->0001->0010-> …………,else CODE_next = (CODE_current+1) << 1 until
% 		   len(CODE_next) = next code size,01->100,0010->00110->001110-> …………
% 		3、Terminate when going through	HUFFSIZE of code size = 0.
%
% 	"""
HuffCodes = zeros(1,length(HuffSize));  % 码字的个数
idx=1;                                  % 码字索引
CODE = 0;                               % 码字
SIZE = HuffSize(1);                     % 码长
while 1
    HuffCodes(idx) = CODE;          % 按照第一步，码字从0开始
    if idx == length(HuffCodes) || HuffSize(idx) == 0
        % 遍历完成或遇到HuffSize中出现0时退出
        break
    end
    idx = idx + 1;
    if HuffSize(idx) == SIZE
        %{
                    第二步，如果码长和之前的相同，则只是码字+1，不做其他操作
        %}
        CODE = CODE + 1;
    else
        %{
                    如果码长和前一个不同，则码字+1后左移，差几位左移几次
        %}
        CODE = bitshift(CODE+1,HuffSize(idx)-SIZE);
    end
    SIZE = HuffSize(idx);
end
end