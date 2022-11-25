function  [EHUFCO,EHUFSI] = SortHuffTbl(HuffmanTable)
%{
    EHUFCO和EHUFSI中，待编码值和码字、码长一一对应，
    可以直接用值来索引码长、码字
%}
[BITS,HuffVal] = HuffmanTable{:};

HUFFSIZE = decode_HUFFSIZE(BITS);
HUFFCODE = decode_HUFFCODE(HUFFSIZE);

EHUFCO = zeros(1,max(HuffVal)+1);
EHUFSI = zeros(1,max(HuffVal)+1);

for K = 1:length(HuffVal)
    value = HuffVal(K);
    EHUFCO(value+1) = HUFFCODE(K);
    EHUFSI(value+1) = HUFFSIZE(K);
end
end