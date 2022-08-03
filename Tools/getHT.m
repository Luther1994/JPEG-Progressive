function  [HuffTable_DC, HuffTable_AC] = getHT(HS,HV_DC,HV_AC)
% 		"""Get Huffman table
%
% 		Returns
% 			HUFFSIZE Huffman Size Table which states the size of huffman code
% 			HV Huffman Value which states the corresponding symbols of huffman tree
%
% 		"""
HuffTable_DC = zeros(4,12);
HuffTable_AC = zeros(4,251);
for identifier = 1:4      % ['00', '10', '01', '11']
    [HUFFSIZE, HUFFCODE] = Decode_Hufftbl(HS(identifier,:));
    if identifier == 1
        HUFFVAL = HV_DC(1,:);
    elseif identifier == 2
        HUFFVAL = HV_AC(1,:);
    elseif identifier == 3
        HUFFVAL = HV_DC(2,:);
    else
        HUFFVAL = HV_AC(2,:);
    end
     
    SIZE = HUFFVAL(end) + 1;
    EHUFCO = zeros(1,SIZE);
    EHUFSI = zeros(1,SIZE);
    for K = 1:length(HUFFVAL)
        Index = HUFFVAL(K);
        EHUFCO(Index+1) = HUFFCODE(K);
        EHUFSI(Index+1) = HUFFSIZE(K);
    end
    if identifier == 1
        HuffTable_DC(1,:) = EHUFCO;
        HuffTable_DC(3,:) = EHUFSI;
    elseif identifier == 2
        HuffTable_AC(1,:) = EHUFCO;
        HuffTable_AC(3,:) = EHUFSI;
    elseif identifier == 3
        HuffTable_DC(2,:) = EHUFCO;
        HuffTable_DC(4,:) = EHUFSI;
    else
        HuffTable_AC(2,:) = EHUFCO;
        HuffTable_AC(4,:) = EHUFSI;
    end
end
end