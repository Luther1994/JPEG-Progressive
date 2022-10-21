function Tbl = Generate_HuffTable(Values)
% ==========================================
% Huffman码表的生成函数，P144~p148。fig K.1~4
% ==========================================

% 所有存在的待编码symbol
ValidValues = sort(unique(Values));

% 每个symbol所在树的另外一枝
OTHERS = zeros(1,256)-1;

% 码长为对应值的码字的个数
CODESIZE = zeros(1,256);

% 所有symbol的概率
Freq = zeros(1,256);
Freq(256) = 1;
for k = 1:length(ValidValues)
    Val = ValidValues(k);
    Freq(k) = sum(Values==Val)/numel(Values);
end

% 将待编码值按照概率大小排序，从而出现概率最大的值用最短的码长编码
[~,index] = sort(Freq(1:length(ValidValues)),'descend');
ValidValues = ValidValues(index);

% =======================================================
% 这一步只需要获取对应码长的码字的个数，不需要得到具体的码字
% =======================================================
while 1
    least_val = min(Freq(Freq>0));

    % 找到最小和次小概率所在的位置
    index1 = find(Freq==least_val);
    if length(index1)>1
        index2 = index1(2);
        index1 = index1(1);
    else
        next_val = min(Freq(Freq>least_val));
        if isempty(next_val)
            break
        end
        index2 = find(Freq==next_val,1);
    end

    % 概率求和赋值到最小概率的位置，次小概率置0
    Freq(index1) = Freq(index1)+Freq(index2);
    Freq(index2) = 0;

    % 相应的位置码长+1
    CODESIZE(index1) = CODESIZE(index1)+1;
    CODESIZE(index2) = CODESIZE(index2)+1;

    % 循环左子树
    sign = 1;
    while sign
        if OTHERS(index1)==-1
            sign = 0;
        else
            index1 = OTHERS(index1);
            CODESIZE(index1) = CODESIZE(index1)+1;
        end
    end

    % 链接到右子树并循环
    OTHERS(index1)=index2;
    while ~sign
        if OTHERS(index2) ==-1
            sign = 1;
        else
            index2 = OTHERS(index2);
            CODESIZE(index2) = CODESIZE(index2)+1;
        end
    end
end
BITS = CountBits(CODESIZE);
Tbl = {BITS,ValidValues};
end

function BITS = CountBits(CODESIZE)

% 对应长度的码字个数。
% 假设所有的symbol都有足够大的概率，因此不会
% 有码长大于32bit的码字存在，所以BITS最长为32
BITS = zeros(1,32);

% 码字统计
for l = 1:256
    if CODESIZE(l)~=0
        BITS(CODESIZE(l)) = BITS(CODESIZE(l))+1;
    end
end

% 调整码长，将码长大于16的码字都压缩到16bit以内
BITS = AdjustBits(BITS);
end

function [BITS] = AdjustBits(BITS)
% ===========================================
% 码长调整函数
% ===========================================

% 从最长的码字向下循环
I = 32;
while 1

    % 如果存在码长为I的码字，则对码长分布进行调整
    if BITS(I) > 0
        J = I-1;
        while 1
            J = J-1;
            if BITS(J) > 0
                break
            end
        end
        BITS(I) = BITS(I)-2;
        BITS(I-1) = BITS(I-1)+1;
        BITS(J+1) = BITS(J+1)+2;
        BITS(J) = BITS(J)-1;
    else
        I = I-1;
        if I == 16
            while BITS(I) == 0
                I = I - 1;
            end
            BITS(I) = BITS(I) - 1;
            BITS = BITS(1:16);
            break
        end
    end
end
end