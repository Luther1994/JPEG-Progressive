function Encode(src,savefile,quality,method)
%{
    The main function of Encoding.
%}
if ~exist("method",'var')
    % Ĭ�ϵı���ģʽ���������۽�����
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


% ��׼�ο�������
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

% ���������һ���еĲ����ֱ�Ϊͨ����ţ�-1��ʾȫͨ������Ss��Se��Ah��Al�����ﰴ��
% һ���۽�����Ĺ���һ��ͼ��Ϊ10�ν��б���룬ÿ�α����Ӧͨ��DCTϵ����Ss:SeƵ
% ����Ah:Alλ��

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

% �������ṹ��
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
% ���ﱣ��dcϵ������ʱ��block����˳������Ӧ�Զ�ͨ������Ľ���洢ģʽ
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
        % д��Huffman�����Marker��FFC4����ͷ
        WriteTwoBytes(Markers.DHT)
        [BITS,VALUES] = table{:};
        datalength = length(BITS) + length(VALUES) + 3;
        WriteTwoBytes(datalength)
        WriteOneByte(id)
        WriteOneByte(BITS)
        WriteOneByte(VALUES)
        %{
            Huffman��Ĵ洢��ʽ�����볤��BITS��ֵ��VALUES���볤��һ��ռ16���ֽڣ�
            ��˳��洢�볤Ϊ1��16�����ֵĸ����������VALUE��ĳ��ȵ���BITS������
            ����������BITS������һһ��Ӧ����ʾ�����ֱ���ķ���ֵ��
            
            ������ݶεĳ����ǹ̶��ģ� �����ֽڴ洢���ݳ��ȡ�һ���ֽڴ洢Huffman
            ���ţ�16���ֽڴ洢�볤��BITS��ÿ���볤ռһ���ֽڣ�VALUES��ĳ�����
            BITS�����ݾ�����ÿ��������һ���ֽڴ洢��
        %}
    end

    function  WriteSOS(comps)
        %{
            SOS���ݶι涨��ÿ�α�������������
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
            ��һ����Ҫָ��ÿ��ͨ�����õ�Huffman���洢��ʽΪ�����ֽ����ݳ��ȡ�
            һ�ֽ�ͨ������һ�ֽ�ͨ����š�һ�ֽ�Huffman���š�һ�ֽڴ洢0����
            ʾZigZag����ʱ�Ŀ�ʼ��š�һ�ֽڴ洢63��ʾZigZag����Ľ�����š�һ
            �ֽڴ洢00���ֱ�ΪAhֵ��Alֵ���Ի���DCT�ı�����̣�������ֵ�̶�Ϊ0.
        %}
    end

    function  WriteAppn()
        % ================================================
        % д��ΪӦ�û�����������Ϣ����Marker'FFE0-FFEF'��ͷ
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

        Encoder.LastDCVal = zeros(ImgInfo.components,1); % DCϵ�����

        % ���������
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
        % ÿ��ɨ�衪��������һ��SOS������ָ��ͨ����ָ��Ƶ����ָ�����ȵ�ϵ��
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
        % DCϵ���ĵ�һ�α��룬AlΪ�������λ��
        % ===================================
        C = ImgInfo.components;

        % DCϵ�����룬������Ǵ洢��Ӧϵ����Ҫ�Ķ��������ֵĳ���
        CodeLength = zeros(C,Encoder.MaxBlocks);

        % DCϵ�����뵱ǰblockϵ����ǰһ��blockϵ���Ĳ�
        DIFF = ones(C,Encoder.MaxBlocks);
        coefs = bitshift(Encoder.coefs,-Encoder.Al,'int16');
        blk_id = zeros(1,C);
        for mcu_cnt = 0:Encoder.MCUs-1
            for blkn = 1:Encoder.blks_in_mcu
                ci = Encoder.member_ship(blkn);
                blk_cnt = Encoder.dc_blk_id(sum(blk_id)+1);
                blk_id(ci) = blk_id(ci) + 1;
                % ����DCϵ����ֵ
                DcDiff = coefs(ci,1,blk_cnt) - Encoder.LastDCVal(ci);

                % ȷ�����ֳ���
                DIFF(ci,blk_id(ci)) = DcDiff;
                CodeLength(ci,blk_id(ci)) = EnsureGategory(DcDiff);
                Encoder.LastDCVal(ci) = coefs(ci,1,blk_cnt);
                
            end
        end
        
        % д��DC���ȱ�����Ǳ���RGBͼ�Ļ���д��AC���ȱ�
        LuTbl = Generate_HuffTable(CodeLength(1,:));
        WriteDHT(0,LuTbl);
        if C > 1
            codes = [CodeLength(2,1:Encoder.component(2).blocks)...
                CodeLength(3,1:Encoder.component(3).blocks)];
            ChrTbl = Generate_HuffTable(codes);
            WriteDHT(1,ChrTbl);
        end

        % д��SOSɨ��ͷ��������ǰ�������ͨ������Ϣ
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
        % DCϵ���ľ��ȱƽ�
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
        % ACϵ���ĵ�һ�α��룬����ACϵ���ĸ�λ���õ�����
        % ==========================================
        comp = Encoder.component(ComponentId);

        % ��ʼƵ������ֹƵ�� & �����bitλ��ȷ������ʹ�����ϵ��
        band_width = Encoder.Se-Encoder.Ss+1;
        ACcoefs = squeeze(Encoder.coefs(ComponentId,Encoder.Ss+1:Encoder.Se+1,:));
        ACcoefs = idivide(ACcoefs,2^Encoder.Al);

        % ȫΪ0��block���������۽��������ACϵ������ͨ�����������
        run_of_length = 0;

        % symbol to be encoded��RS for AC coefficient
        HUFFVAL = [];

        % NUM of blocks in current component

        for blk_id =  1:comp.blocks
            ZZ = ACcoefs(:,blk_id);
            R = 0;      % num of continuous 0 in block
            for K = 1:band_width
                if ZZ(K)
                    if run_of_length

                        % ȫ 0 block�ĸ������ܳ���2^14-1����RS=15/0��ռ����
                        % ��ʾZRL�룩�����������������һ��RS=14/0��ʾ32767
                        if run_of_length > 32767
                            HUFFVAL(end+1) = bitshift(14,4);
                            run_of_length = run_of_length-32767;
                        end

                        % ȡ2����������ȡ�����õ�����run_of_length��Ҫ��bit��
                        EOBRUN = bitshift(floor(log2(run_of_length)),4);
                        HUFFVAL(end+1) = EOBRUN;
                    end

                    % ����ϵ��ֵȷ�������ϵ����Ҫ���볤S,��ϳ��γ���RS���������
                    S = EnsureGategory(ZZ(K));
                    while R >= 16
                        RRRRSSSS = bitshift(15,4);
                        HUFFVAL(end+1) = RRRRSSSS;
                        R = R - 16;
                    end
                    RRRRSSSS = bitshift(R,4) + S;
                    HUFFVAL(end+1) = RRRRSSSS;
                    % ���ַ�0ϵ����ȫ0ϵ����block��������Ϊ0,block������0
                    % �ĸ�������
                    run_of_length = 0;
                    R = 0;

                    % block������0�ĸ���+1
                else
                    R = R + 1;
                end
            end
            if R
                % run_of_length == 1������ʱ��ֱ��������ǰblock������ǰblock
                % �����ϵ��ȫΪ0
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
        % ACϵ���ĵڶ��α��룬����ACϵ���ĵ�λ���õ�ϸ��
        % ==========================================
        comp = Encoder.component(ComponentId);

        % ��ʼƵ������ֹƵ�� & �����bitλ��ȷ������ʹ�����ϵ��
        assert(Encoder.Ah - Encoder.Al == 1,'1 bit would be encoded once.')

        MSB = squeeze(idivide(Encoder.coefs(ComponentId,:,:),2^Encoder.Ah));
        LSB = idivide(squeeze(Encoder.coefs(ComponentId,:,:)) ...
            - MSB*2^Encoder.Ah,2^Encoder.Al);
        run_of_length = 0;
        % symbol to be encoded��RS for AC coefficient
        HUFFVAL = [];
        for blk_id = 1:comp.blocks
            ZZ = LSB(:,blk_id);
            History = MSB(:,blk_id);

            % Num of continuous 0 in block��ignore 0 with non-zero history
            % while counting.
            R = 0;

            for K = Encoder.Ss + 1:Encoder.Se + 1
                if History(K)==0
                    if ZZ(K)
                        if run_of_length

                            % ȫ 0 block�ĸ������ܳ���2^14-1����RS=15/0��ռ����
                            % ��ʾZRL�룩�����������������һ��RS=14/0��ʾ32767
                            if run_of_length > 32767
                                HUFFVAL(end+1) = bitshift(14,4);
                                run_of_length = run_of_length-32767;
                            end

                            % ȡ2����������ȡ�����õ�����run_of_length��Ҫ��bit��
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
                % run_of_length == 1������ʱ��ֱ��������ǰblock������ǰblock
                % �����ϵ��ȫΪ0
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
        % ��ACrefine�ı����У�����ͳ�Ƽ��������0�ĸ�����ʱ������˵�ǰλ�õ�ϵ��
        % ǰһ���Ȳ�Ϊ0������������ڵõ�run_of_length֮��Ҫ��ͷͳ�Ʊ����Ե�ϵ����
        % �����ض�����������б��롣SkippedCode���������洢��Щ�����Ե�ϵ���ı���
        % ����ġ�ÿ�γ���RS��0ʱ�����ع�ͷ����һ�Σ�Ȼ���ٳ�ʼ��Ϊ0.
        % =================================================================
        [SkippedCode,SkippedSize] = deal(0);

        for blk_id = 1:comp.blocks
            R = 0;
            ZZ = LSB(:,blk_id);
            History = MSB(:,blk_id);
            K = Encoder.Ss + 1;

            % ͬ�����ں��Ե�ԭ������RS����Ҫ���±�����ǰblock��ϵ������ǰһ
            % ���Ȳ�Ϊ0��ϵ�����б��룬LastK����һ��RS��0ʱ��ϵ��λ��������ÿ
            % һ��block��LastK����ʼ��ΪSs+1.���������Historyָ��ǰλ�õ�ϵ
            % ��ǰһ�����µ�ֵ��Current��ʾ��ϵ���ڵ�ǰ�����µ�ֵ��
            LastK = Encoder.Ss + 1;
            while K <= Encoder.Se + 1
                if History(K) == 0
                    if ZZ(K) == 0
                        % History == 0 && Current == 0
                        R = R + 1;
                    else
                        % History == 0 && Current �� 0
                        if Encoder.run_of_length
                            [Code,Size] = Encode_run_length(EHUFCO,EHUFSI);
                            UpdateBuffer(Code,Size);

                            % �ȱ���Eobrun���ٰ��м��Թ��Ĳ��ּ���ȥ��Ȼ�����á�
                            UpdateBuffer(SkippedCode,SkippedSize);
                            [SkippedCode,SkippedSize] = deal(0);
                            Encoder.run_of_length = 0;
                        end
                        while R >= 16
                            % R >=16 ʱ����RS == 15/0��ʾ������16��0��
                            % ͬ����Ҫͳ��16��0�б����ԵĲ���
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
                        % ���ݵ�ǰϵ������д���Ӧ����������
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
            % �������Ĳ��ֱ���
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
        % ��������ȫΪ0��ҲҪ����֮��д���ļ��С�
        if Encoder.run_of_length
            [Code,Size] = Encode_run_length(EHUFCO,EHUFSI);
            UpdateBuffer(Code,Size);
        end
        FillByte();
    end

    function UpdateBuffer(code,size)
        % ��һ����Ϊ�˱��������������ֹ��󣬵����������
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
        % ��������д���ļ���
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
        %  �������һ���ֽڵ�Bits�������
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

        % �������run_length��Ҫ��bit����Ȼ��תΪEOBn���б���
        code = EHUFCO(EOBn+1);
        size = EHUFSI(EOBn+1);

        % Run of Length����2^eobrun�Ĳ��֣�ֱ�����볤Ϊeobrun�Ķ����Ʊ�ʾ
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