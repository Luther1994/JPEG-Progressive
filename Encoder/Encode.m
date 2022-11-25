function Encode(src,savefile,quality,method)
%{
    The main function of Encoding.
%}
if ~exist("method",'var')
    % 默认的编码模式，这里是累进编码
    method = 2;
end
fid1 = fopen(savefile, 'wb');

BLOCKSIZE = 8;          %  block size of one 
DCTSIZE = 64;           %  Num of coeficients in a block
if max(max(src)) > 255
    precision  = 16;
else
    precision = 8;
end
if ndims(src) == 3 && size(src,3) == 3
    src = int16(rgb2yuv(src))-2^(precision-1);
    [Height,Width,C] = size(src);
else
    src = int16(src)  - 2^(precision -1);
    [Height,Width] = size(src);
    C = 1;
end
if mod(Height,BLOCKSIZE) || mod(Width,BLOCKSIZE)
    addition_rows = Height-mod(Height,BLOCKSIZE);
    addition_cols = Width - mod(Width,BLOCKSIZE);
    Height = Height + addition_rows;
    Width = Width + addition_cols;
    src = padarray(src,[addition_rows,addition_cols],'symmetric','post');
end


% 标准参考量化表
qtbl ={
    [[16, 11, 10, 16, 24, 40, 51, 61]
    [12, 12, 14, 19, 26, 58, 60, 55]
    [14, 13, 16, 24, 40, 57, 69, 56]
    [14, 17, 22, 29, 51, 80, 87, 62]
    [18, 22, 37, 56, 68, 109, 103, 77]
    [24, 35, 55, 64, 81, 104, 113, 92]
    [49, 64, 78, 87, 103, 121, 120, 121]
    [72, 92, 95, 98, 112, 100, 103, 99]];

    [[17, 18, 24, 47, 99, 99, 99, 99]
    [18, 21, 26, 66, 99, 99, 99, 99]
    [24, 26, 56, 99, 99, 99, 99, 99]
    [47, 66, 99, 99, 99, 99, 99, 99]
    [99, 99, 99, 99, 99, 99, 99, 99]
    [99, 99, 99, 99, 99, 99, 99, 99]
    [99, 99, 99, 99, 99, 99, 99, 99]
    [99, 99, 99, 99, 99, 99, 99, 99]]};
[LQ,pre1] = GetQuantizer(qtbl{1},quality);
[CQ,pre2] = GetQuantizer(qtbl{2},quality);
if C == 1
    SF = [1;1];
else
    SF = [2 1 1;2 1 1];
end

% 编码参数，一行中的参数分别为通道编号（-1表示全通道）、Ss、Se、Ah、Al，这里按照
% 一般累进编码的管理将一张图分为10次进行编解码，每次编码对应通道DCT系数的Ss:Se频
% 带的Ah:Al位。

Scans = {
    -1     0     0     0     1
    1     1     5     0     2
    -1     0     0     1     0
    3     1    63     0     1
    2     1    63     0     1
    1     6    63     0     2
    1     1    63     2     1
    3     1    63     1     0
    2     1    63     1     0
    1     1    63     1     0
    };

Markers = struct( ...
    %{
    ============================================
    Start of Frame.
    non-differential,Huffman coding
    ============================================
    %}
'SOF0',65472,...  % Baseline DCT
'SOF1',65473,...  % Extended sequential DCT
'SOF2',65474,...  % Progressive DCT
'SOF3',65475,...  % Lossless (sequential)
%{
    ============================================
    differential,Huffman coding
    ============================================
%}
'SOF5',65477,...  % Differential sequential DCT
'SOF6',65478,...  % Differential progressive DCT
'SOF7',65479,...  % Differential lossless (sequential)
%{
    ============================================
    non-differential,arithmetic coding
    ============================================
%}
'JPG',65480,...   % Reserved for JPEG Extensions
'SOF9',65481,...  % Extended sequential DCT
'SOF10',65482,... % Progreddive DCT
'SOF11',65483,... % Lossless (sequential)
%{
    ============================================
    differential,arithemtic coding
    ============================================
%}
'SOF13',65485,... % Differential sequential DCT
'SOF14',65486,... % Differential Progressive DCT
'SOF15',65487,... % Differential lossless(sequential)
%{
        Huffman table specification
%}
'DHT',65476,...   % Define Huffman table
%{
    ============================================
    Arithmetic coding conditioning specification
    ============================================
%}
'DAC',65484,...    % Define arithmetic coding conditioning(s)
%{
    ============================================
    Restart interval termination
    ============================================
%}
'RSTm',65488:65495,...        % Restart with modulo 8 count 'm'
%{
    ============================================
    Other markers
    ============================================
%}
'SOI',65496,...   % Start of Marker
'EOI',65497,...   % End of Image
'SOS',65498,...   % Start of Scan
'DQT',65499,...   % Define Quantization Table
'DNL',65500,...   % Define Number of lines
'DRI',65501,...   % Define restart interval
'DHP',65502,...   % Define hierarchical progression
'EXP',65503,...   % Expand reference conponent(s)
'APPn',65504:65519,... % Reserved for application segments
'JPGn',65520:65533,... % Reserved for JPEG extensions
'COM',65534,...   % Comment
%{
    ============================================
    Reserved markers
    ============================================
%}
'TEM',655281,... & For temporary private use in arithmetic coding
'RES',655282:655471);% Reserved

Comp = struct(...
    'component_id',0,...             % identifier for this component (0..255)
    'w_samp_factor',0,...            % sampling factor in width (1..4)
    'h_samp_factor',0,...            % sampling factor in height (1..4)
    'qtbl_index',0,...               % quantization table selector (0..3)
    'dc_tbl_no',0,...                % DC entropy table selector (0..3)
    'ac_tbl_no',0,...                % AC entropy table selector (0..3)
    'blocks_per_row',0,...           % num of blocks for current component,horizontally
    'blocks_per_col',0,...           % num of blocks for current component,vertically
    'downsampled_width',0,...        % actual component width in samples
    'downsampled_height',0,...       % actual component height in samples
    'MCU_width',0,...                % number of blocks per MCU, horizontally
    'MCU_height',0,...               % number of blocks per MCU, vertically
    'MCU_blocks',0,...               % MCU_width * MCU_height
    'last_cols',0,...                % non-dummy blocks across in last MCU
    'last_rows',0);                  % non-dummy blocks down in last MCU

% 编码器结构体
Encoder = struct(...
    'dqt_precision',[pre1,pre2],...
    'sample_factor',SF,...
    'Ss',0,...
    'Se',0,...
    'Ah',0,...
    'Al',0,...
    'blocks',0,...
    'MaxBlocks',0,...
    'savefile',savefile, ...
    'Buffer',[]);

Encoder.quanti_tbl = [zigzag(LQ) zigzag(CQ)];

% Image infomation
ImgInfo = struct(...
    'Precision',precision,...
    'components',C,...
    'Height',Height,...
    'Width',Width);

SFw_max = max(Encoder.sample_factor(1,:));
SFh_max = max(Encoder.sample_factor(2,:));
Encoder.MCUs_per_row = ceil(Width/SFw_max/BLOCKSIZE);
Encoder.MCUs_per_col = ceil(Height/SFh_max/BLOCKSIZE);
Encoder.MCUs = Encoder.MCUs_per_col*Encoder.MCUs_per_row;
Encoder.blks_in_mcu = sum(prod(Encoder.sample_factor));
for c = 1:ImgInfo.components
    ImgInfo.com_id(c) = c;
    Comp.component_id = c;
    sfw = Encoder.sample_factor(1,c);
    sfh = Encoder.sample_factor(2,c);
    Comp.w_samp_factor = sfw;
    Comp.h_samp_factor = sfh;
    if c == 1
        Comp.qtbl_index = 0;
        Comp.dc_tbl_no = 0;
        Comp.ac_tbl_no = 0;
    else
        Comp.qtbl_index = 1;
        Comp.ac_tbl_no = 1;
        Comp.dc_tbl_no = 1;
    end
    factors = [SFw_max/sfw,SFh_max/sfh];
    comp = DOWNSAMPLE(src(:,:,c),factors);
    [Height,Width] = size(comp);

    % blocks in height/width & num of blocks in one MCU
    Comp.MCU_height = sfh;
    Comp.MCU_width = sfw;
    Comp.MCU_blocks = sfh*sfw;

    % blocks in row/col and blocks in whole
    Comp.blocks_per_row = ceil(Width/BLOCKSIZE);
    Comp.blocks_per_col = ceil(Height/BLOCKSIZE);
    Comp.blocks = Comp.blocks_per_col*Comp.blocks_per_row;

    Encoder.blocks = Encoder.blocks + Comp.blocks;

    % remainder of cols/rows devided by BLOCKSIZE(8)
    Comp.last_cols = Comp.blocks_per_row*BLOCKSIZE-Width;
    Comp.last_rows = Comp.blocks_per_col*BLOCKSIZE-Height;

    % height/width of component to encode after padding and downsampling
    Comp.downsampled_height = Comp.blocks_per_col * BLOCKSIZE;
    Comp.downsampled_width = Comp.blocks_per_row * BLOCKSIZE;
    Comp.component = padarray(comp,[Comp.last_rows,Comp.last_cols],'symmetric','post');

    % initialized mem space to save coefficient
    Comp.coes = zeros(DCTSIZE,Comp.blocks_per_row*Comp.blocks_per_col);
    Encoder.component(c) = Comp;
end
Encoder.coefs = zeros(ImgInfo.components,DCTSIZE,Encoder.MCUs*SFh_max*SFw_max,'int16');
% 这里保存dc系数解码时的block排列顺序，用于应对多通道编码的交错存储模式
Encoder.dc_blk_id = ones(1,Encoder.blocks);
GetMemeberShip()

blkn = 1;
for ROW = 0:Encoder.MCUs_per_col - 1
    for COL = 0:Encoder.MCUs_per_row - 1
        for ci = 1:ImgInfo.components
            McuWidth = Encoder.sample_factor(1,ci);
            McuHeight = Encoder.sample_factor(2,ci);
            comp = Encoder.component(ci);
            for r = 0:McuHeight - 1
                for c = 0:McuWidth - 1
                    row = ROW * McuHeight + r;
                    col = COL * McuWidth + c;
                    blk_id = row * comp.blocks_per_row + col + 1;
                    Encoder.dc_blk_id(blkn) = blk_id;
                    blkn = blkn + 1;
                end
            end
        end
    end
end
    function GetMemeberShip()
        Encoder.member_ship = zeros(1,Encoder.blks_in_mcu);
        index = 1;
        for cc = 1:ImgInfo.components
            mcublks = Encoder.component(cc).MCU_blocks;
            while mcublks > 0
                Encoder.member_ship(index) = cc;
                index = index + 1;
                mcublks = mcublks - 1;
            end
        end
    end
    function WriteOneByte(Bytes)
        for i = 1:length(Bytes)
            Byte = Bytes(i);
            Encoder.Buffer(end+1) = Byte;
        end
    end
    function WriteTwoBytes(Bytes)
        for i = 1:length(Bytes)
            Byte = Bytes(i);
            MSB = bitshift(Byte,-8);
            LSB = bitand(Byte,255);
            WriteOneByte(MSB);
            WriteOneByte(LSB);
        end
    end
WriteSOI = @() WriteTwoBytes(Markers.SOI);
WriteEOI = @() WriteTwoBytes(Markers.EOI);

EncodeImage();

    function  EncodeImage()
        WriteSOI()
        WriteAppn()
        WriteDQT(0)
        WriteDQT(1)
        WriteSOF()
        EncodeRestartInterval()
        WriteEOI()
        Close()
    end

    function  WriteDQT(identifier)

        % ============================================
        % write DQT marker into file
        % ============================================
        WriteTwoBytes(Markers.DQT);

        % Determin which quanti table and precision to be write
        qtb = Encoder.quanti_tbl(:,identifier+1);
        precision = Encoder.dqt_precision(identifier+1);
        ComponentId = bitor(bitshift(precision,4),identifier);
        datalength = length(qtb) + 2 + 1;
        WriteTwoBytes(datalength);
        WriteOneByte(ComponentId);

        if precision
            WriteTwoBytes(qtb);
        else
            WriteOneByte(qtb)
        end
    end

    function  WriteSOF()

        % ===============================================================
        %  Write Start of Frame based on different encode mode.This part
        %  contains the essential information of image. See B.2.2 in P.35
        % ===============================================================

        % Choose encode method(Baseline or Progressive)
        if method == 2
            SOF = Markers.SOF2;
        end
        WriteTwoBytes(SOF)
        datalength = 8 + ImgInfo.components * 3;
        WriteTwoBytes(datalength)
        WriteOneByte(ImgInfo.Precision)
        WriteTwoBytes(ImgInfo.Height)
        WriteTwoBytes(ImgInfo.Width)
        WriteOneByte(ImgInfo.components)
        if ImgInfo.components == 1
            WriteOneByte(1);
            WriteOneByte(17);
            WriteOneByte(0);
        else
            for channel_id = 1: ImgInfo.components
                WriteOneByte(channel_id);
                sw = Encoder.sample_factor(1, channel_id);
                sh = Encoder.sample_factor(2, channel_id);
                SamFct = bitshift(sw,4)+sh;
                WriteOneByte(SamFct);
                WriteOneByte(Encoder.component(channel_id).qtbl_index)
            end
        end
    end

    function  WriteDHT(id,table)
        % 写入Huffman码表，以Marker‘FFC4’开头
        WriteTwoBytes(Markers.DHT)
        [BITS,VALUES] = table{:};
        datalength = length(BITS) + length(VALUES) + 3;
        WriteTwoBytes(datalength)
        WriteOneByte(id)
        WriteOneByte(BITS)
        WriteOneByte(VALUES)
        %{
            Huffman表的存储格式是先码长表BITS后值表VALUES，码长表一共占16个字节，
            按顺序存储码长为1到16的码字的个数，后面的VALUE表的长度等于BITS中码字
            的总数，与BITS中码字一一对应，表示该码字编码的符号值。
            
            这个数据段的长度是固定的， 两个字节存储数据长度、一个字节存储Huffman
            表编号；16个字节存储码长表BITS，每个码长占一个字节；VALUES表的长度由
            BITS表内容决定，每个数据用一个字节存储。
        %}
    end

    function  WriteSOS(comps)
        %{
            SOS数据段规定了每次编解码的数据内容
        %}
        WriteTwoBytes(Markers.SOS)
        C = length(comps);
        datalength = 6 + 2 * C;
        WriteTwoBytes(datalength)
        WriteOneByte(C);
        for i = 1:C
            ComponentId = comps(i);
            WriteOneByte(ComponentId)
            dc_id = Encoder.component(ComponentId).dc_tbl_no;
            ac_id = Encoder.component(ComponentId).ac_tbl_no;
            huffid = bitor(bitshift(dc_id,4),ac_id);
            WriteOneByte(huffid);
        end
        WriteOneByte(Encoder.Ss);
        WriteOneByte(Encoder.Se);
        PtTra = bitshift(Encoder.Ah,4)+Encoder.Al;
        WriteOneByte(PtTra)
        %{
            这一段主要指定每个通道采用的Huffman表，存储格式为：两字节数据长度、
            一字节通道数、一字节通道序号、一字节Huffman表编号、一字节存储0，表
            示ZigZag编码时的开始序号、一字节存储63表示ZigZag编码的结束序号、一
            字节存储00，分别为Ah值和Al值，对基于DCT的编码过程，这两个值固定为0.
        %}
    end

    function  WriteAppn()
        % ================================================
        % 写入为应用环境保留的信息，以Marker'FFE0-FFEF'开头
        % ================================================
        WriteTwoBytes(Markers.APPn(1))
        WriteTwoBytes(16)
        WriteOneByte([74, 70, 73, 70,  0,  1,  1,  0,  0,  1,  0,  1,  0,  0])
    end

    function  EncodeRestartInterval()
        for chan = 1:ImgInfo.components
            blkn = 1;
            comp = Encoder.component(chan);
            qtbl = Encoder.quanti_tbl(:,comp.qtbl_index+1);
            for R  = 0:comp.blocks_per_col-1
                for C = 0:comp.blocks_per_row-1
                    block = comp.component(R*BLOCKSIZE+1:(R+1)*BLOCKSIZE,...
                        C*BLOCKSIZE+1:(C+1)*BLOCKSIZE);
                    coef = zigzag(dct2(block)) ./ qtbl;
                    Encoder.coefs(chan,:,blkn) = coef;
                    blkn = blkn + 1;
                end
            end
        end

%         COEF = Encoder.coefs;
%         save("coe.mat","COEF");

        Encoder.LastDCVal = zeros(ImgInfo.components,1); % DC系数差分

        % 缓存的码流
        Encoder.BufferedVal = 0;
        Encoder.BufferedBits = 0;

%         Encoder.coefs  = load('coe.mat').COEF;
        for s = 1:size(Scans,1)
            [Comp,Encoder.Ss,Encoder.Se,...
                Encoder.Ah,Encoder.Al] = Scans{s,:};
            EncodeScan(Comp);
        end
    end

    function EncodeScan(ComponentId)
        % 每次扫描————即一次SOS，编码指定通道、指定频带、指定精度的系数
        if Encoder.Ss == 0
            if Encoder.Ah == 0
                EncodeDCFirst();
            else
                EncodeDCRefine();
            end
        else
            if Encoder.Ah == 0
                EncodeACFirst(ComponentId);
            else
                EncodeACRefine(ComponentId);
            end
        end
    end

    function EncodeDCFirst()
        % ===================================
        % DC系数的第一次编码，Al为编码最低位数
        % ===================================
        C = ImgInfo.components;

        % DC系数编码，编码的是存储对应系数需要的二进制码字的长度
        CodeLength = zeros(C,Encoder.MaxBlocks);

        % DC系数编码当前block系数和前一个block系数的差
        DIFF = ones(C,Encoder.MaxBlocks);
        coefs = bitshift(Encoder.coefs,-Encoder.Al,'int16');
        blk_id = zeros(1,C);
        for mcu_cnt = 0:Encoder.MCUs-1
            for blkn = 1:Encoder.blks_in_mcu
                ci = Encoder.member_ship(blkn);
                blk_cnt = Encoder.dc_blk_id(sum(blk_id)+1);
                blk_id(ci) = blk_id(ci) + 1;
                % 计算DC系数差值
                DcDiff = coefs(ci,1,blk_cnt) - Encoder.LastDCVal(ci);

                % 确定码字长度
                DIFF(ci,blk_id(ci)) = DcDiff;
                CodeLength(ci,blk_id(ci)) = EnsureGategory(DcDiff);
                Encoder.LastDCVal(ci) = coefs(ci,1,blk_cnt);
                
            end
        end
        
        % 写入DC亮度表，如果是编码RGB图的话再写入AC亮度表
        LuTbl = Generate_HuffTable(CodeLength(1,:));
        WriteDHT(0,LuTbl);
        if C > 1
            codes = [CodeLength(2,1:Encoder.component(2).blocks)...
                CodeLength(3,1:Encoder.component(3).blocks)];
            ChrTbl = Generate_HuffTable(codes);
            WriteDHT(1,ChrTbl);
        end

        % 写入SOS扫描头，包含当前被编码的通道的信息
        WriteSOS(1:ImgInfo.components);
        blk_id = ones(1,C);
        for mcu_cnt = 0:Encoder.MCUs-1
            for blkn = 1:Encoder.blks_in_mcu
                ci = Encoder.member_ship(blkn);
                if ci == 1
                    hufftable = LuTbl;
                else
                    hufftable = ChrTbl;
                end
                DcDiff = DIFF(ci,blk_id(ci));
                codelength= CodeLength(ci,blk_id(ci));
                [code,len] = EncodeDC(DcDiff,codelength,hufftable);
                UpdateBuffer(code,len);
                WriteCodeStream();
                blk_id(ci) = blk_id(ci) + 1;
            end
        end
        FillByte();
    end

    function EncodeDCRefine()
        % =======================================
        % DC系数的精度逼近
        % =======================================
        C = ImgInfo.components;
        WriteSOS(1:ImgInfo.components);
        LSB = squeeze(bitand(Encoder.coefs(:,1,:),2^Encoder.Al,'int16'));
        for mcu_cnt = 0:Encoder.MCUs-1
            for blkn = 1:Encoder.blks_in_mcu
                ci = Encoder.member_ship(blkn);
                blk_cnt = Encoder.dc_blk_id(mcu_cnt*Encoder.blks_in_mcu+blkn);
                Encoder.BufferedVal = ...
                    bitshift(Encoder.BufferedVal,1);
                Encoder.BufferedBits = Encoder.BufferedBits + 1;
                if LSB(ci,blk_cnt)
                    Encoder.BufferedVal = Encoder.BufferedVal + 1;
                end
            end
            WriteCodeStream();
        end
        FillByte();
    end


    function EncodeACFirst(ComponentId)
        % ==========================================
        % AC系数的第一次编码，编码AC系数的高位，得到轮廓
        % ==========================================
        comp = Encoder.component(ComponentId);

        % 初始频带和终止频带 & 编码的bit位，确定带宽和待编码系数
        band_width = Encoder.Se-Encoder.Ss+1;
        ACcoefs = squeeze(Encoder.coefs(ComponentId,Encoder.Ss+1:Encoder.Se+1,:));
        ACcoefs = idivide(ACcoefs,2^Encoder.Al);

        % 全为0的block的数量，累进编码编码AC系数与普通编码最大区别
        run_of_length = 0;

        % symbol to be encoded，RS for AC coefficient
        HUFFVAL = [];

        % NUM of blocks in current component

        for blk_id =  1:comp.blocks
            ZZ = ACcoefs(:,blk_id);
            R = 0;      % num of continuous 0 in block
            for K = 1:band_width
                if ZZ(K)
                    if run_of_length

                        % 全 0 block的个数不能超过2^14-1个（RS=15/0被占用于
                        % 表示ZRL码），如果超过，则先用一个RS=14/0表示32767
                        if run_of_length > 32767
                            HUFFVAL(end+1) = bitshift(14,4);
                            run_of_length = run_of_length-32767;
                        end

                        % 取2对数并向下取整，得到编码run_of_length需要的bit数
                        EOBRUN = bitshift(floor(log2(run_of_length)),4);
                        HUFFVAL(end+1) = EOBRUN;
                    end

                    % 根据系数值确定编码该系数需要的码长S,组合成游程码RS存入待编码
                    S = EnsureGategory(ZZ(K));
                    while R >= 16
                        RRRRSSSS = bitshift(15,4);
                        HUFFVAL(end+1) = RRRRSSSS;
                        R = R - 16;
                    end
                    RRRRSSSS = bitshift(R,4) + S;
                    HUFFVAL(end+1) = RRRRSSSS;
                    % 出现非0系数后，全0系数的block计数重置为0,block内连续0
                    % 的个数归零
                    run_of_length = 0;
                    R = 0;

                    % block内连续0的个数+1
                else
                    R = R + 1;
                end
            end
            if R
                % run_of_length == 1，解码时会直接跳过当前block，即当前block
                % 后面的系数全为0
                run_of_length = run_of_length + 1;
            end
            if blk_id == comp.blocks
                if run_of_length
                    HUFFVAL(end+1) = bitshift(floor(log2(run_of_length)),4);
                end
            end
        end
        Htable = Generate_HuffTable(HUFFVAL);
        if ComponentId == 1
            WriteDHT(16,Htable);
        else
            WriteDHT(17,Htable);
        end
        WriteSOS(ComponentId);
        [EHUFCO,EHUFSI] = SortHuffTbl(Htable);
        Encoder.run_of_length = 0;
        for blk_id = 1:comp.blocks
            R = 0;
            ZZ = ACcoefs(:,blk_id);
            for K = 1:band_width
                if ZZ(K)
                    if Encoder.run_of_length
                        [Code,Size] = Encode_run_length(EHUFCO,EHUFSI);
                        UpdateBuffer(Code,Size);
                        Encoder.run_of_length = 0;
                    end
                    while R >=16
                        UpdateBuffer(EHUFCO(241),EHUFSI(241));
                        R = R - 16;
                    end
                    S = EnsureGategory(ZZ(K));
                    RS = bitshift(R,4) + S;
                    UpdateBuffer(EHUFCO(RS+1),EHUFSI(RS+1));
                    if ZZ(K) < 0
                        ZZ(K) = ZZ(K) +  2 ^ S - 1;
                    end
                    UpdateBuffer(ZZ(K),S);
                    R = 0;
                    if K == band_width
                        break;
                    end
                else
                    R = R + 1;
                    if K == band_width
                        Encoder.run_of_length = Encoder.run_of_length + 1;
                        if Encoder.run_of_length == 32767
                            [Code,Size] = Encode_run_length(EHUFCO,EHUFSI);
                            UpdateBuffer(Code,Size);
                            Encoder.run_of_length = Encoder.run_of_length - 32767;
                            break;
                        end
                    end
                end
            end
            if blk_id == comp.blocks
                if Encoder.run_of_length
                    [Code,Size]=Encode_run_length(EHUFCO,EHUFSI);
                    UpdateBuffer(Code,Size);
                end
            end
        end
        FillByte();
    end

    function EncodeACRefine(ComponentId)
        % ==========================================
        % AC系数的第二次编码，编码AC系数的低位，得到细节
        % ==========================================
        comp = Encoder.component(ComponentId);

        % 初始频带和终止频带 & 编码的bit位，确定带宽和待编码系数
        assert(Encoder.Ah - Encoder.Al == 1,'1 bit would be encoded once.')

        MSB = squeeze(idivide(Encoder.coefs(ComponentId,:,:),2^Encoder.Ah));
        LSB = idivide(squeeze(Encoder.coefs(ComponentId,:,:)) ...
            - MSB*2^Encoder.Ah,2^Encoder.Al);
        run_of_length = 0;
        % symbol to be encoded，RS for AC coefficient
        HUFFVAL = [];
        for blk_id = 1:comp.blocks
            ZZ = LSB(:,blk_id);
            History = MSB(:,blk_id);

            % Num of continuous 0 in block，ignore 0 with non-zero history
            % while counting.
            R = 0;

            for K = Encoder.Ss + 1:Encoder.Se + 1
                if History(K)==0
                    if ZZ(K)
                        if run_of_length

                            % 全 0 block的个数不能超过2^14-1个（RS=15/0被占用于
                            % 表示ZRL码），如果超过，则先用一个RS=14/0表示32767
                            if run_of_length > 32767
                                HUFFVAL(end+1) = bitshift(14,4);
                                run_of_length = run_of_length-32767;
                            end

                            % 取2对数并向下取整，得到编码run_of_length需要的bit数
                            EOBRUN = bitshift(floor(log2(run_of_length)),4);
                            HUFFVAL(end+1) = EOBRUN;
                        end

                        S = EnsureGategory(ZZ(K));
                        assert(S == 1,'Size of newly coef in refine AC must be 1.')

                        while R >= 16
                            RRRRSSSS = bitshift(15,4);
                            HUFFVAL(end+1) = RRRRSSSS;
                            R = R - 16;
                        end
                        RRRRSSSS = bitshift(R,4) + S;
                        HUFFVAL(end+1) = RRRRSSSS;
                        run_of_length = 0;
                        R = 0;
                    else
                        R = R + 1;
                    end
                end
            end
            if R
                % run_of_length == 1，解码时会直接跳过当前block，即当前block
                % 后面的系数全为0
                run_of_length = run_of_length + 1;
            end
            if blk_id == comp.blocks
                if run_of_length
                    HUFFVAL(end+1) = bitshift(floor(log2(run_of_length)),4);
                end
            end
        end
        Htable = Generate_HuffTable(HUFFVAL);
        if ComponentId == 1
            WriteDHT(16,Htable);
        else
            WriteDHT(17,Htable);
        end
        WriteSOS(ComponentId);
        [EHUFCO,EHUFSI] = SortHuffTbl(Htable);
        Encoder.run_of_length = 0;

        % =================================================================
        % 在ACrefine的编码中，由于统计间隔的连续0的个数的时候忽略了当前位置的系数
        % 前一精度不为0的情况，所以在得到run_of_length之后要回头统计被忽略的系数，
        % 并用特定方法对其进行编码。SkippedCode就是用来存储这些被忽略的系数的编码
        % 结果的。每次出现RS≠0时，都回过头处理一次，然后再初始化为0.
        % =================================================================
        [SkippedCode,SkippedSize] = deal(0);

        for blk_id = 1:comp.blocks
            R = 0;
            ZZ = LSB(:,blk_id);
            History = MSB(:,blk_id);
            K = Encoder.Ss + 1;

            % 同样由于忽略的原因，遇到RS后需要重新遍历当前block的系数，对前一
            % 精度不为0的系数进行编码，LastK是上一次RS≠0时的系数位置索引。每
            % 一个block，LastK都初始化为Ss+1.编码过程中History指当前位置的系
            % 数前一精度下的值，Current表示该系数在当前精度下的值。
            LastK = Encoder.Ss + 1;
            while K <= Encoder.Se + 1
                if History(K) == 0
                    if ZZ(K) == 0
                        % History == 0 && Current == 0
                        R = R + 1;
                    else
                        % History == 0 && Current ≠ 0
                        if Encoder.run_of_length
                            [Code,Size] = Encode_run_length(EHUFCO,EHUFSI);
                            UpdateBuffer(Code,Size);

                            % 先编码Eobrun，再把中间略过的部分加上去，然后重置。
                            UpdateBuffer(SkippedCode,SkippedSize);
                            [SkippedCode,SkippedSize] = deal(0);
                            Encoder.run_of_length = 0;
                        end
                        while R >= 16
                            % R >=16 时，用RS == 15/0表示连续的16个0，
                            % 同样需要统计16个0中被忽略的部分
                            UpdateBuffer(EHUFCO(241),EHUFSI(241));
                            k = LastK;
                            i = 1;
                            while i <= 16
                                if History(k) ~= 0
                                    if ZZ(k)
                                        UpdateBuffer(1,1);
                                    else
                                        UpdateBuffer(0,1);
                                    end
                                else
                                    i = i + 1;
                                end
                                k = k + 1;
                            end
                            R = R - 16;
                            LastK = k;
                        end
                        RS = bitshift(R,4) + 1;   % S shuld always be 1.
                        UpdateBuffer(EHUFCO(RS+1),EHUFSI(RS+1));
                        % 根据当前系数符号写入对应二进制码字
                        if ZZ(K) > 0
                            UpdateBuffer(1,1);
                        else
                            UpdateBuffer(0,1);
                        end
                        for k = LastK:K
                            if History(k) ~= 0
                                if ZZ(k)
                                    UpdateBuffer(1,1);
                                else
                                    UpdateBuffer(0,1);
                                end
                            else
                                if R == 0
                                    break
                                end
                                R = R - 1;
                            end
                        end
                        LastK = K + 1;
                    end
                end
                K = K + 1;
            end
            % 对跳过的部分编码
            for k = LastK:Encoder.Se + 1
                if MSB(k,blk_id)
                    SkippedCode = bitshift(SkippedCode,1);
                    SkippedSize = SkippedSize + 1;
                    if LSB(k,blk_id)
                        % With correction
                        SkippedCode = SkippedCode +1;
                    end
                end
            end
            if R
                Encoder.run_of_length = Encoder.run_of_length + 1;
                if Encoder.run_of_length == 32767
                    [Code,Size] = Encode_run_length(EHUFCO,EHUFSI);
                    UpdateBuffer(Code,Size);
                    Encoder.run_of_length = Encoder.run_of_length - 32767;
                    break;
                end
            end
        end
        % 如果最后几列全为0，也要编码之后写入文件中。
        if Encoder.run_of_length
            [Code,Size] = Encode_run_length(EHUFCO,EHUFSI);
            UpdateBuffer(Code,Size);
        end
        FillByte();
    end

    function UpdateBuffer(code,size)
        % 这一步是为了避免编码出来的码字过大，导致数据溢出
        while size > 8
            MSB = bitshift(code,8-size);
            Encoder.BufferedVal = bitshift(Encoder.BufferedVal,8)+MSB;
            Encoder.BufferedBits = Encoder.BufferedBits+8;
            WriteCodeStream();
            size = size-8;
            code = bitand(code,bitshift(1,size)-1);
        end
        Encoder.BufferedVal = bitshift(Encoder.BufferedVal,size)+code;
        Encoder.BufferedBits = Encoder.BufferedBits + size;
        WriteCodeStream();
    end

    function WriteCodeStream()
        % ==================
        % 将编码结果写入文件中
        % ==================

        while Encoder.BufferedBits >= 8
            CC = bitshift(Encoder.BufferedVal, ...
                8-Encoder.BufferedBits);
            WriteOneByte(CC);
            if CC == 255
                WriteOneByte(0);
            end
            Encoder.BufferedVal = bitand(Encoder.BufferedVal, ...
                bitshift(1,Encoder.BufferedBits-8)-1);
            Encoder.BufferedBits = Encoder.BufferedBits - 8;
        end
    end

    function FillByte()
        %  =============================
        %  对最后不足一个字节的Bits进行填充
        %  =============================
        if Encoder.BufferedBits
            Encoder.BufferedVal = bitshift(Encoder.BufferedVal, ...
                8-Encoder.BufferedBits);
            WriteOneByte(Encoder.BufferedVal)
            Encoder.BufferedBits = 0;
            Encoder.BufferedVal = 0;
        end
    end


    function [code,size] = Encode_run_length(EHUFCO,EHUFSI)
        eobrun = floor(log2(Encoder.run_of_length));
        EOBn = bitshift(eobrun,4);

        % 计算编码run_length需要的bit数，然后转为EOBn进行编码
        code = EHUFCO(EOBn+1);
        size = EHUFSI(EOBn+1);

        % Run of Length超过2^eobrun的部分，直接用码长为eobrun的二进制表示
        Extension = Encoder.run_of_length-bitshift(1,eobrun);

        if eobrun
            code = bitshift(code,eobrun) + Extension;
            size = size + eobrun;
        end
    end
    function Close()
        fwrite(fid1,Encoder.Buffer,"uint8");
        fclose(fid1);
    end
end