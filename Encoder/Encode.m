function Encode(src,savefile,quality,method)
%{
    The main function of Encoding.
%}
if ~exist("method",'var')
    method = 2;
end
fid = fopen(savefile, 'wb');     
BLOCKSIZE = 8;
DCTSIZE = 64;
if max(max(src)) > 255
    precision  = 16;
else
    precision = 8;
end
if ndims(src) == 3 && size(src,3) == 3
    src = int16(rgb2yuv(src))-2^(precision-1);
    [Height,Width,Components] = size(src);
else
    src = int16(src)  - 2^(precision -1);
    [Height,Width] = size(src);
    Components = 1;
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
if Components == 1
    SF = [1;1];
else
    SF = [2 1 1;2 1 1];
end

% �������ṹ��
Encoder = struct(...
    'dqt_precision',[pre1,pre2],...
    'sample_factor',SF,...
    'Ss',0,...
    'Se',0,...
    'Ah',0,...
    'Al',0,...
    'run_length',0,...
    'bitsleft',0,...
    'TemporaryBuffer',0,...
    'SOIFound',0, ...
    'savefile',savefile, ...
    'EOB_RL',0, ...
    'Re',char());

Encoder.quanti_tbl{1} = LQ;
Encoder.quanti_tbl{2} = CQ;


% Image infomation
ImgInfo = struct(...
    'Precision',precision,...
    'components',Components,...
    'Height',Height,...
    'Width',Width);

Markers = struct( ...
    %{
        Start of Frame.
        non-differential,Huffman coding
    %}
'SOF0',65472,...  % Baseline DCT
'SOF1',65473,...  % Extended sequential DCT
'SOF2',65474,...  % Progressive DCT
'SOF3',65475,...  % Lossless (sequential)
%{
        differential,Huffman coding
%}
'SOF5',65477,...  % Differential sequential DCT
'SOF6',65478,...  % Differential progressive DCT
'SOF7',65479,...  % Differential lossless (sequential)
%{
        non-differential,arithmetic coding
%}
'JPG',65480,...   % Reserved for JPEG Extensions
'SOF9',65481,...  % Extended sequential DCT
'SOF10',65482,... % Progreddive DCT
'SOF11',65483,... % Lossless (sequential)
%{
        differential,arithemtic coding
%}
'SOF13',65485,... % Differential sequential DCT
'SOF14',65486,... % Differential Progressive DCT
'SOF15',65487,... % Differential lossless(sequential)
%{
        Huffman table specification
%}
'DHT',65476,...   % Define Huffman table
%{
        Arithmetic coding conditioning specification
%}
'DAC',65484,...    % Define arithmetic coding conditioning(s)
%{
        Restart interval termination
%}
'RSTm',65488:65495,...        % Restart with modulo 8 count 'm'
%{
     Other markers
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
    Reserved markers
%}
'TEM',655281,... & For temporary private use in arithmetic coding
'RES',655282:655471);% Reserved

c = struct(...
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

SFw_max = max(Encoder.sample_factor(1,:));
SFh_max = max(Encoder.sample_factor(2,:));
Encoder.MCUs_in_row = ceil(Width/SFw_max/BLOCKSIZE);
Encoder.MCUs_in_col = ceil(Height/SFh_max/BLOCKSIZE);

for i = 1:ImgInfo.components
    c.component_id = i;
    sfw = Encoder.sample_factor(1,i);
    sfh = Encoder.sample_factor(2,i);
    c.w_samp_factor = sfw;
    c.h_samp_factor = sfh;
    if i == 1
        c.qtbl_index = 0;
        c.dc_tbl_no = 0;
        c.ac_tbl_no = 0;
    else
        c.qtbl_index = 1;
        c.ac_tbl_no = 1;
        c.dc_tbl_no = 1;
    end
    factors = [SFw_max/sfw,SFh_max/sfh];
    comp = DOWNSAMPLE(src(:,:,i),factors);
    [Height,Width] = size(comp);
    
    % blocks in height/width & num of blocks in one MCU
    c.MCU_height = sfh;
    c.MCU_width = sfw;
    c.MCU_blocks = sfh*sfw;

    % blocks in row/col
    c.blocks_per_row = Encoder.MCUs_in_row * sfw;
    c.blocks_per_col = Encoder.MCUs_in_col * sfh;
    
    % remainder of cols/rows devided by BLOCKSIZE(8)
    c.last_cols = c.blocks_per_row*BLOCKSIZE-Width;
    c.last_rows = c.blocks_per_col*BLOCKSIZE-Height;

    % height/width of component to encode after padding and downsampling
    c.downsampled_height = c.blocks_per_col * BLOCKSIZE;
    c.downsampled_width = c.blocks_per_row * BLOCKSIZE;
    c.component = padarray(comp,[c.last_rows,c.last_cols],'symmetric','post');

    % initialized mem space to save coefficient
    c.coes = zeros(DCTSIZE,c.blocks_per_row*c.blocks_per_col);
    Encoder.component(i) = c;
end
WriteOneByte = @(Bits) WriteNBytes(fid,Bits,1);
WriteTwoBytes = @(Bits) WriteNBytes(fid,Bits,2);
WriteFourBits = @(Bits) WriteNBytes(fid,Bits,0.5);

EncodeImage(); % ��ʼ����
    function  appendSOI()
        % write SOI marker into file
        WriteTwoBytes(Markers.SOI);
    end

    function  appendEOI()
        % write EOI marker into file
        WriteTwoBytes(Markers.EOI);
    end

    function  appendDQT(identifier)
        % write DQT marker into file
        WriteTwoBytes(Markers.DQT);

        % Determin which quanti table and precision to be write
        qtb = zigzag(Encoder.quanti_tbl{identifier+1});
        precision = Encoder.dqt_precision(identifier+1);
        i = bitor(bitshift(precision,4),identifier);
        datalength = length(qtb) + 2 + 1;
        WriteTwoBytes(datalength);
        WriteOneByte(i);

        if precision
            WriteTwoBytes(qtb);
        else
            WriteOneByte(qtb)
        end
    end

    function  appendSOF()
        %{
           Write Start of Frame based on different encode mode.
           This part contains the essential information of image.
           See B.2.2 in P.35
        %}
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
            WriteFourBits(1);
            WriteFourBits(1);
            WriteOneByte(0);
        else
            for channel_id = 1: ImgInfo.components
                WriteOneByte(channel_id);
                WriteFourBits(Encoder.sample_factor(1, channel_id))
                WriteFourBits(Encoder.sample_factor(2, channel_id))
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

    function  WriteSOS(channels)
        %{
            SOS���ݶι涨��ÿ�α���봦���
        %}
        WriteTwoBytes(Markers.SOS)
        datalength = 6 + 2 * channels;
        WriteTwoBytes(datalength)
        WriteOneByte(channels)
        for i = 1:channels
            WriteOneByte(i)
            dc_id = Encoder.component(i).dc_tbl_no;
            ac_id = Encoder.component(i).ac_tbl_no;
            c = bitor(bitshift(dc_id,4),ac_id);
            WriteOneByte(c);
        end
        WriteOneByte(Encoder.Ss);
        WriteOneByte(Encoder.Se);
        WriteFourBits(Encoder.Ah);
        WriteFourBits(Encoder.Al);
        %{
            ��һ����Ҫָ��ÿ��ͨ�����õ�Huffman���洢��ʽΪ�����ֽ����ݳ��ȡ�
            һ�ֽ�ͨ������һ�ֽ�ͨ����š�һ�ֽ�Huffman���š�һ�ֽڴ洢0����
            ʾZigZag����ʱ�Ŀ�ʼ��š�һ�ֽڴ洢63��ʾZigZag����Ľ�����š�һ
            �ֽڴ洢00���ֱ�ΪAhֵ��Alֵ���Ի���DCT�ı�����̣�������ֵ�̶�Ϊ0.
        %}
    end

    function  appendAPPn()
        % д��ΪӦ�û�����������Ϣ��
        % ��Marker'FFE0-FFEF'��ͷ
        WriteTwoBytes(Markers.APPn(1))
        WriteTwoBytes(16)
        WriteOneByte([74, 70, 73, 70,  0,  1,  1,  0,  0,  1,  0,  1,  0,  0])
    end

    function  ResetEncoder()
        Encoder.LastDCVal = zeros(ImgInfo.components,1);
    end

    function  EncodeImage()
        appendSOI()
        appendAPPn()
        appendDQT(0)
        appendDQT(1)
        appendSOF()
        EncodeRestartInterval()
        appendEOI()
        close()
    end
    function  EncodeRestartInterval()
        ResetEncoder()
        % ����ϵ���������ǰ��շǽ���ģʽ���ɵģ������ͨ�������ʱ��ǵð���
        % ����ģʽȥȡ����
%         for chan = 1:ImgInfo.components
%             comp = Encoder.component(chan);
%             qtbl = Encoder.quanti_tbl{comp.qtbl_index+1};
%             for R  = 0:comp.blocks_per_col-1
%                 for C = 0:comp.blocks_per_row-1
%                     block_id = R*comp.blocks_per_row + C + 1;
%                     block = comp.component(R*BLOCKSIZE+1:(R+1)*BLOCKSIZE,...
%                         C*BLOCKSIZE+1:(C+1)*BLOCKSIZE);
%                     coef = round(dct2(block) ./ double(qtbl));
%                     coef = zigzag(coef);
%                     comp.coes(:,block_id) = coef;
%                 end
%             end
%             Encoder.component(chan) = comp;
%         end
%         COEF = cell(3,1);
%         for cc = 1:ImgInfo.components
%             COEF{cc} = Encoder.component(cc).coes;
%         end
%         save("coe.mat","COEF");
        coef  = load('coe.mat').COEF;
        for j = 1:ImgInfo.components
            Encoder.component(j).coes = coef{j};
        end
        EncodeDCFirst(0);
        EncodeACFirst(1);
%         EncodeDCRefine();
    end
    function EncodeDCFirst(Al)
        % ===================================
        % DCϵ���ĵ�һ�α��룬AlΪ�������λ��
        % ===================================
        ResetEncoder();
        Encoder.Ss = 0 ;
        Encoder.Se = 0 ;
        Encoder.Ah = 0 ;
        Encoder.Al = 0 ;
        channels = ImgInfo.components;
        CodeLength = cell(1,channels);     
        DIFF = cell(1,channels);    
        for c_id = 1:channels
            comp = Encoder.component(c_id);
            coefs = comp.coes;
            coefs = bitshift(coefs(1,:),-Al,'int16');
            blocks = comp.blocks_per_col*comp.blocks_per_row;
            huffval = zeros(1,blocks);
            diff = zeros(1,blocks);
            index = 1;
            for row = 0:Encoder.MCUs_in_col-1
                for col = 0:Encoder.MCUs_in_row-1
                    for r = 0:comp.MCU_height-1
                        for c = 0: comp.MCU_width-1
%                             bandid = comp.blocks_per_row*comp.MCU_height*row...
%                                 +col*comp.MCU_width+...
%                                 r*comp.blocks_per_row+c+1;
                            d = coefs(index) - Encoder.LastDCVal(c_id);
                            % ȷ�������ֵ��Ҫ�����ֳ��ȣ�Ҳ��ʵ�ʵı������
                            diff(index) = d;
                            huffval(index) = EnsureGategory(d);
                            Encoder.LastDCVal(c_id) = coefs(index);
                            index = index + 1;
                        end
                    end
                end
            end
            CodeLength{c_id} = huffval; % ����ͨ���Ĵ�����ֵ
            DIFF{c_id} = diff;
        end
        
        LuTbl = Generate_HuffTable(CodeLength{1});
        WriteDHT(0,LuTbl);

        if channels > 1
            ChrTbl = Generate_HuffTable([CodeLength{2:3}]);
            WriteDHT(1,ChrTbl);
        end
        WriteSOS(ImgInfo.components);
        % Interleaved mode should be adopted while encoding multicomponents 
        % ���ͨ��ͬʱ���б����ʱ����Ҫ���ս���ģʽ�洢
        
        for chan = 1:ImgInfo.components
            if chan == 1
                hufftable = LuTbl;
            else
                hufftable = ChrTbl;
            end
            comp = Encoder.component(chan);
            blocks = comp.blocks_per_col*comp.blocks_per_row;
            Re =[];
            for bandid = 1:blocks
                d = DIFF{chan}(bandid);
                codelength = CodeLength{chan}(bandid);
                [code,len] = EncodeDC(d,codelength,hufftable);
                Re = [Re dec2bin(code,len)];
                while length(Re) >= 8
                    byte = bin2dec(Re(1:8));
                    WriteOneByte(byte)
                    if byte == 255
                        WriteOneByte('00')
                    end
                    Re = Re(9:end);
                end
            end
        end
        if ~isempty(Re) % �������һ���ֽڵ�Bits�������
            Re = PrepareforMarker(Re);
            WriteOneByte(bin2dec(Re))
        end
    end

    function EncodeDCRefine()
        ResetEncoder();
        Encoder.Ss = 0;
        Encoder.Se = 0;
        Encoder.Ah = 1;
        Encoder.Al = 0;
        %         appendDHT('11');
        coefs = string(bitand(Encoder.COEFICIENTS{:}(1,:),1*2^Encoder.Al));
        tableid = [0,0,0];
        WriteSOS(tableid);
        % �������DCϵ��֮����Ҫ����Ƶ�ʴӵ͵����𲽿�ʼ����ACϵ����ÿ�α���ѡ
        % ����Ƶ�������ղ�ͬ��ͨ���ֱ���б���
        while length(coefs) >= 8
            byte = bin2dec(join(coefs(1:8)));
            WriteOneByte(byte)
            %{
               ÿ�δ��ȥһ���ֽڣ�������FF����Ҫ�ں������ '00'����ֹ
               ��Marker��ͻ
            %}
            if byte == 255
                WriteOneByte('00')
            end
            coefs = coefs(9:end);
        end
        if ~isempty(coefs)    % �������һ���ֽڵ�Bits�������
            coefs = PrepareforMarker(coefs);
            WriteOneByte(bin2dec(coefs))
        end
    end

    function EncodeACFirst(channelid)
        % ======================================
        % ACϵ���ĵ�һ�α��룬����ACϵ���ĸ�λ��
        % �õ�����
        % ======================================
        
        comp = Encoder.component(channelid);
        
        % ��ʼƵ������ֹƵ�� & �����bitλ��ȷ������ʹ�����ϵ��
        Encoder.Ss = 1;
        Encoder.Se = 63;
        Encoder.Ah = 0;
        Encoder.Al = 0;
        band_width = Encoder.Se-Encoder.Ss+1;
        ACcoefs = Encoder.component(channelid).coes;  
        ACcoefs = idivide(int16(ACcoefs(Encoder.Ss+1:Encoder.Se+1,:)),2^Encoder.Al);
%         fid = fopen('COE.txt','w');
%         for i = 1:32400
%             fprintf(fid,'%-6i',i);
%         end  
%         fprintf(fid,'\n\n');
%         for i = 1:63
%             for j = 1:32400
%                 fprintf(fid,'%-6i',ACcoefs(i,j));
%             end
%             fprintf(fid,'\n');
%         end
%         fclose(fid);
        % ȫΪ0��block���������۽��������ACϵ������ͨ�����������
        run_of_length = 0;

        % symbol to be encoded��RS for AC coefficient
        HUFFVAL = [];
        
        % NUM of blocks in current component 
        blocks = comp.blocks_per_col* comp.blocks_per_row;

        for bandid =  1:blocks
            ZZ = ACcoefs(:,bandid);
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
        end
        
        % ��ȡHuffman���
        Htable = Generate_HuffTable(HUFFVAL);
        fid = fopen('code.txt','w');
        if channelid == 1
            WriteDHT('10',Htable);
        else
            WriteDHT('11',Htable);
        end
        WriteSOS(1);
        
        [EHUFCO,EHUFSI] = SortHuffTbl(Htable);
        Re = '';
        Encoder.run_of_length = 0;
        for bandid = 1:blocks
            code = '';
            R = 0;
            ZZ = ACcoefs(:,bandid);
            for K = 1:band_width
                if ZZ(K)
                    if Encoder.run_of_length
                        code = [code Encode_run_length(EHUFCO,EHUFSI)];
                    end
                    while R >=16
                        code =[code, dec2bin(EHUFCO(241),EHUFSI(241))];
                        R = R - 16;
                    end
                    code = [code Encode_R_ZZ(R,ZZ(K))];
                    R = 0;
                    if K == band_width
                        Re = [Re code];
                        break;
                    end
                else
                    R = R + 1;
                    if K == band_width
                        Re = [Re code];
                        Encoder.run_of_length = Encoder.run_of_length + 1;
                        if Encoder.run_of_length == 32767
                            code=Encode_run_length(EHUFCO,EHUFSI);
                            Re =[Re code];
                            Encoder.run_of_length = Encoder.run_of_length-32767;
                            break;
                        end
                    end
                end
            end
            while length(Re) >= 8
                byte = bin2dec(Re(1:8));
                WriteOneByte(byte)
                if byte == 255
                    WriteOneByte('00')
                end
                Re = Re(9:end);
            end
        end
        if ~isempty(Re) % �������һ���ֽڵ�Bits�������
            Re = PrepareforMarker(Re);
            WriteOneByte(bin2dec(Re))
        end
        function code = Encode_R_ZZ(R,Coe)
            SSSS = EnsureGategory(Coe);
            RS = bitshift(R,4) + SSSS;
            huffmancode = EHUFCO(RS+1);
            huffmansize = EHUFSI(RS+1);
            huffmancode = dec2bin(huffmancode,huffmansize);
            if Coe<0
                Coe = Coe +  2 ^ SSSS - 1;
            end
            code = [huffmancode dec2bin(Coe,SSSS)];
        end
        function huffmancode = Encode_run_length(EHUFCO,EHUFSI)
            eobrun = floor(log2(Encoder.run_of_length));
            EOBn = bitshift(eobrun,4);
            % �������run_length��Ҫ��bit����Ȼ��תΪEOBn���б���
            huffmancode = EHUFCO(EOBn+1);
            huffmansize = EHUFSI(EOBn+1);
            huffmancode = dec2bin(huffmancode,huffmansize);
            if eobrun
                huffmancode = [huffmancode dec2bin(Encoder.run_of_length-bitshift(1,eobrun),eobrun)];
            end
            Encoder.run_of_length = 0;
        end
    end
    function EncodeACRefine()
        
    end
end