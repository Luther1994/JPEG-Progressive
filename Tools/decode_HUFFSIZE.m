function [HUFFSIZE]=decode_HUFFSIZE(Bits)
%{
            Generate the table of Huffman code size
            通过Bits获得Huffman码长表
        	Args:
        		Bits: Table containing the number of codes of each size
        			  eg.Bits[2]=5 means there exists 5 code words with length of 3(counting starts in 0)
                      
        
        	Returns:
        		HUFFSIZE: Table of Huffman code size whose length represent number of existed code words
        				  and value represent size of code words
        
        	Notes:
        		Length of Bits should always be 16,represent 16 possible code sizes,Bits[i] will be set to 0
        		if no code word with code size i.
                注意这里另外一个地方，码长为i的码字个数n需满足n<2^(i-1)-2^(i-2),
                即不能出现全为1的码字，这个码字要保留作为下一个长度的码字的前缀，
                因为按照码字的生成方法，长度变化时上一个码字+1然后左移，如果
                上一个码字全为1，+1之后长度直接就是下一个码长了。
%}
HUFFSIZE = zeros(1,sum(Bits));  % 一共有多少个码字
idx =1;                         % 码字的索引
CODESIZE = 1;                   % 码长，1~16
COUNTS = 0;                     % 码长为指定CODESIZE的码字个数
%{
		    COUNTS is the number of code words with code size CODESIZE
		    which_channel is the index of HUFFSIZE table
%}
while CODESIZE <= 16
    while Bits(CODESIZE) > COUNTS
        % 码长为CODESIZE的码字个数为Bits(CODESIZE)
        HUFFSIZE(idx) = CODESIZE; %第idx个码字的码长为CODESIZE
        idx =idx+ 1;
        COUNTS =COUNTS+ 1;
    end
    CODESIZE = CODESIZE + 1;        % 下一个码长，个数初始化为1
    COUNTS = 0;
end
end