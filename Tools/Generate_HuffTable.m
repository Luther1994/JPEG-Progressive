function Tbl = Generate_HuffTable(Values)
%{
    动态Huffman编码
%}

ValidValues = sort(unique(Values));
FREQ = zeros(1,256);
OTHERS = zeros(1,256)-1;
CODESIZE = zeros(1,256);
for i = 1:length(ValidValues)
    Val = ValidValues(i);
    FREQ(i) = sum(Values==Val)/numel(Values); % 按照顺序每个条目代表一个值
end
while 1
    least_val = min(FREQ(FREQ>0));
    index1 = find(FREQ==least_val);
    if least_val == 1
        break
    end
    if length(index1)>1
        % 如果一个概率对应几个不同的值，那么V1为最大的V值
        Vals = sort(ValidValues(index1));
        V1 = Vals(end);
        V2 = Vals(end-1);
        index1 = find(ValidValues==V1);
        index2 = find(ValidValues==V2);
    else
        next_val = min(FREQ(FREQ>least_val));
        if next_val ~=1
            index2 = find(FREQ==next_val);
        else
            break
        end
    end
    FREQ(index1) = FREQ(index1)+FREQ(index2);
    FREQ(index2) = 0;
    CODESIZE(index1) = CODESIZE(index1)+1;
    while OTHERS(index1) ~=-1
        index1 = OTHERS(index1);
        CODESIZE(index1) = CODESIZE(index1)+1;
    end
    OTHERS(index1)=index2;
    CODESIZE(index2) = CODESIZE(index2)+1;
    while OTHERS(index2)~=-1
        index2 = OTHERS(index2);
        CODESIZE(index2) = CODESIZE(index2)+1;
    end
end   
BITS = CountBits(CODESIZE);
Tbl = {BITS,ValidValues};
    function BITS = CountBits(CODESIZE)
        BITS = zeros(1,256);
        for l = 1:256
            if CODESIZE(l)~=0
                BITS(CODESIZE(l)) = BITS(CODESIZE(l))+1;
            end
        end
        if sum(BITS(17:end))>0
            BITS = AdjustBits(BITS);
        end
        BITS = BITS(1:16);
        MaxSymNum = zeros(1,16);
        MaxSymNum(1) = 1;
        for i = 2:16
            MaxSymNum(i) = 2*(MaxSymNum(i-1)+1-BITS(i-1))-1;
            if BITS(i)>MaxSymNum(i)
                BITS(i)=BITS(i)-1;
                BITS(i+1)=BITS(i+1)+1;
            end
        end
    end
    function [BITS] = AdjustBits(BITS)
        I = 32;
        while 1
            if BITS(I)==0
                I = I-1;
                if I==16
                    while BITS(I)==0
                        I=I-1;
                    end
                    BITS(I)=BITS(I)-1; 
                    break
                end
            else
                J = I-1;
                while BITS(J)==0
                    J=J-1;
                end
                BITS(I) = BITS(I)-2;
                BITS(I-1)=BITS(I-1)+1;
                BITS(J+1)=BITS(J+1)+2;
                BITS(J) = BITS(J)-1;
            end
        end
    end
end