function [Re]=EH(Htable)
%{
    编码，将4*N的huffman表编码成存储所用的格式
	Args:
		data: Huffman table specifications
        Returns:
%}
    MinCode = Htable(1,:);                   % 指定长度的最小码字
    MaxCode = Htable(2,:);                   % 指定长度的最大码字
    HuffVal = Htable(4,:);                   % Huffman码的值
    BITS = zeros(1,16);
    for i = 1:length(MinCode)
        if i > 16
            break
        end
        if MinCode(i) == -1 && MaxCode(i)==-1
            BITS(i)=0;
        else
            BITS(i) = MaxCode(i)-MinCode(i)+1;
        end
    end
    [HUFFSIZE,HUFFCODE] = Decode_Hufftbl(BITS);
    Re = {BITS,HuffVal,HUFFSIZE,HUFFCODE};
    % 结果为分别表示对应码长的码字个数、编码系数所需要的码长、
    % 编码该码长需要的码长、编码结果
    % 例如系数DIFF为35，编码35需要6个bit（这部分可以直接算出来）；
    % 编码6需要4个bit，编码结果为1110(这部分是熵编码部分)。
end
