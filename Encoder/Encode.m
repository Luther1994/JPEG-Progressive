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
SF = [2 1 1;2 1 1];
% �������ṹ��
Encoder = struct(...
    'dqt_precision',[pre1,pre2],...
    'QtblIdx',[0,1,1],...
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

component_info = struct(...
    'component_id',0,...             % identifier for this component (0..255)
    'component_index',0,...          % its index in SOF or cinfo->comp_info[]
    'w_samp_factor',0,...            % sampling factor in width (1..4)
    'h_samp_factor',0,...            % sampling factor in height (1..4)
    'quant_tbl_no',0,...             % quantization table selector (0..3)
    'dc_tbl_no',0,...                % DC entropy table selector (0..3)
    'ac_tbl_no',0,...                % AC entropy table selector (0..3)
    'width_in_blocks',0,...          % num of blocks for current component,horizontally
    'height_in_blocks',0,...         % num of blocks for current component,vertically
    'downsampled_width',0,...        % actual component width in samples
    'downsampled_height',0,...       % actual component height in samples
    'MCU_width',0,...                % number of blocks per MCU, horizontally
    'MCU_height',0,...               % number of blocks per MCU, vertically
    'MCU_blocks',0,...               % MCU_width * MCU_height
    'MCU_sample_width',0,...         % MCU width in samples, MCU_width*DCT_[h_]scaled_size
    'last_col_width',0,...           % non-dummy blocks across in last MCU
    'last_row_height',0);            % non-dummy blocks down in last MCU
SFw_max = max(Encoder.sample_factor(1,:));
SFh_max = max(Encoder.sample_factor(2,:));
for c = 1:ImgInfo.components
    component_info.component_id = c;
    sfw = Encoder.sample_factor(1,c);
    sfh = Encoder.sample_factor(2,c);
    component_info.w_samp_factor = sfw;
    component_info.h_samp_factor = sfh;
    if c == 1
        component_info.quant_tbl_no = 0;
        component_info.dc_tbl_no = 0;
        component_info.ac_tbl_no = 0;
    else
        component_info.quant_tbl_no = 1;
        component_info.ac_tbl_no = 1;
        component_info.dc_tbl_no = 1;
    end
    
    factors = [SFw_max/sfw,SFh_max/sfh];
    C = DOWNSAMPLE(src(:,:,1),factors);
    [Height,Width] = size(C);
    Encoder.component_info(c).image = C;
    
end
% Encoder.ComponentSize = [ImgInfo.Height./SFs(1,:);...
%     ImgInfo.Width./SFs(2,:)];
% Encoder.blocks_per_comp = uint16(Encoder.ComponentSize./BLOCKSIZE);
% Encoder.Components = {DOWNSAMPLE(src(:,:,1),SFs(:,1)),...
%     DOWNSAMPLE(src(:,:,2),SFs(:,2)),...
%     DOWNSAMPLE(src(:,:,3),SFs(:,2))};
% Encoder.DUcnts = sum(prod(Encoder.blocks_per_comp,1));

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
        c = bitor(bitshift(precision,4),identifier);
        datalength = length(qtb) + 2 + 1;
        WriteTwoBytes(datalength);
        WriteOneByte(c);

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
        for c = 1: ImgInfo.components
            WriteOneByte(c);
            WriteFourBits(Encoder.sample_factor(1, c))
            WriteFourBits(Encoder.sample_factor(2, c))
            WriteOneByte(Encoder.QtblIdx(c))
        end
    end

    function  appendDHT(id,table)
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

    function  appendSOS(ids)
        Components = length(ids);
        WriteTwoBytes(Markers.SOS)
        datalength = 6 + 2 * Components;
        WriteTwoBytes(datalength)
        WriteOneByte(Components)
        for c = 1:length(ids)
            WriteOneByte(c)
            WriteOneByte(ids(c))
        end
        WriteOneByte(Encoder.Ss);
        WriteOneByte(Encoder.Se);
        WriteFourBits(Encoder.Ah);
        WriteFourBits(Encoder.Al); % Ss = 0, Se = 63
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
        % ѡ��һ�ַ�ʽȥ�洢ϵ�����������ͨ������������÷ǽ���ģʽ�洢
        for chan = 1:ImgInfo.components
            blocks_per_row = ImgInfo.Width/BLOCKSIZE/Encoder.sample_factor(1,chan);
            blocks_in_ver = Encoder.blocks_per_comp(1,chan);
            blocks_in_hor = Encoder.blocks_per_comp(2,chan);
            for row = 0:blocks_in_ver-1
                for col = 0:blocks_in_hor-1
                    block = Encoder.Components{chan}(row * BLOCKSIZE + 1 :(row + 1) * BLOCKSIZE,...
                        col * BLOCKSIZE + 1:(col + 1) * BLOCKSIZE);
                    coef = round(dct2(block) ./ double(Encoder.QLUTbl));
                    coef = zigzag(coef);
                    Encoder.component_coes{chan}(row+1,col+1,:) = coef;
                end
            end
        end
        coef = Encoder.component_coes;
        save("coe.mat","coef");
        EncodeDCFirst();
%         EncodeACFirst(1);
%         EncodeDCRefine();
    end
    function EncodeDCFirst()
        Encoder.Ss = 0;
        Encoder.Se = 0;
        Encoder.Ah = 0;
        Encoder.Al = 1;
        ResetEncoder();
        % ���ͨ��ͬʱ���б����ʱ����Ҫ���ս���ģʽ�洢
        for c = 1:ImgInfo.components
            DCcoefs = Encoder.COEFICIENTS{c}(:,:,1);
            
            DCcoefs = bitshift(DCcoefs(:),-Encoder.Al,'int16');
            blocks = prod(Encoder.blocks_per_comp(:,c)); 
            HUFFVAL = zeros(1,blocks);
            DIFFs = zeros(1,blocks);
            for bandid = 1:blocks
                DIFF = DCcoefs(bandid) - Encoder.LastDCVal(c);
                codelength = EnsureGategory(DIFF);  % ȷ������DIFF��Ҫ��bit��
                HUFFVAL(bandid) = codelength;
                Encoder.LastDCVal(c) = DCcoefs(bandid);
                DIFFs(bandid) = DIFF;
            end
            HuffVals{c} = HUFFVAL; % ����ͨ���Ĵ�����ֵ
            Actual_Val{c} = DIFFs;
        end
        LuTbl = Generate_HuffTable(HuffVals{1});
        ChrTbl = Generate_HuffTable([HuffVals{2:3}]);
        appendDHT('00',LuTbl);
        appendDHT('01',ChrTbl);
        tableid = [0,17,17];
        appendSOS(tableid);
        for c = 1:ImgInfo.components
            if c == 1
                table = LuTbl;
            else
                table = ChrTbl;
            end
            blocks = prod(Encoder.blocks_per_comp(:,c));
            HUFFCODE = zeros(1,blocks);
            HUFFSIZE = zeros(1,blocks);
            for bandid = 1:blocks
                BITS = table{1};
                HUFFVAL = table{2};
                DIFF = Actual_Val{c}(bandid);
                codelength = HuffVals{c}(bandid);
                [code,length] = EncodeDC(DIFF,codelength,{BITS,HUFFVAL});
                HUFFCODE(bandid) = code;
                HUFFSIZE(bandid) = length;
            end
            HUFFCODEs{c} = HUFFCODE;
            HUFFSIZEs{c} = HUFFSIZE;
        end
        Re =[];
        bandid = [1,1,1];
        for c = 1:sum(prod(Encoder.blocks_per_comp,1))
            temp = mod(c,Encoder.MCUs_per_iMCU);
            if (1 <= temp) && (temp<= 4)
                channelid = 1;
            elseif temp == 5
                channelid = 2;
            else
                channelid = 3;
            end
            code = HUFFCODEs{channelid}(bandid(channelid));
            length = HUFFSIZEs{channelid}(bandid(channelid));
            Re = [Re dec2bin(code,length)];
            bandid(channelid) = bandid(channelid) + 1;
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
        appendSOS(tableid);
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
        %{
            ����Ҫ����ACϵ������˼��ȷ�ϴ�����ֵ��Ȼ�����ɶ�Ӧ��huffman��
        %}
        Encoder.Ss = 1;
        Encoder.Se = 5;
        Encoder.Ah = 0;
        Encoder.Al = 2;
        ACcoefs = bitshift(Encoder.COEFICIENTS{channelid}(Encoder.Ss+1:Encoder.Se,:),-Encoder.Al,'int8');
        RSs = zeros(1,prod(Encoder.blocks_per_comp(:,channelid)));
        %{
            ����ɫ��ͨ����HUFFMAN��������������Ҫ��ɣ�Ȼ��ͨ�����շǽ���ģʽ�洢
        %}
        run_of_length = 0;
        HUFFVAL = [];
        band_width = Encoder.Se-Encoder.Ss;
        for bandid = 1:prod(Encoder.blocks_per_comp(:,channelid))
            R = 0;
            ZZ = ACcoefs(:,bandid);
            for K = 1:band_width
                if ZZ(K)
                    if run_of_length
                        if run_of_length > 32767
                            HUFFVAL(end+1) = bitshift(14,4);
                            run_of_length = run_of_length-32767;
                        end
                        EOBRUN = bitshift(floor(log2(run_of_length)),4);
                        HUFFVAL(end+1) = EOBRUN;
                    end
                    S = EnsureGategory(ZZ(K));
                    run_of_length = 0;
                    RS = bitshift(R,4) + S;
                    HUFFVAL(end+1) = RS;
                else
                    R = R + 1;
                end
            end
            if R
                run_of_length = run_of_length + 1;
            end
        end
        Htable = Generate_HuffTable(HUFFVAL);
        if channelid == 1
            appendDHT('10',Htable);
        else
            appendDHT('11',Htable);
        end
        appendSOS(1);
        [EHUFCO,EHUFSI] = SortHuffTbl(Htable);
        % �����е�ACϵ�����д���,��ϵ�����շǽ���ģʽ��������һ��
        band_width = Encoder.Se-Encoder.Ss;
        blocks_in_col = prod(Encoder.blocks_per_comp(1,channelid));
        blocks_in_row = prod(Encoder.blocks_per_comp(2,channelid));
        blocks_in_component = prod(Encoder.blocks_per_comp(:,channelid));
        COEs = bitshift(Encoder.COEFICIENTS{channelid}(Encoder.Ss+1:Encoder.Se,:),-Encoder.Al,'int8');
        non_interleaved_coes = zeros(size(COEs));
        row = 0; 
        rows_cnt = Encoder.sample_factor(2,channelid);
        cols_cnt = Encoder.sample_factor(1,channelid);
        if Encoder.sample_factor(2,channelid) > 1
           while 1
               rows = COEs(:,row*blocks_in_row+1:(row+rows_cnt)*blocks_in_row);
               for id = 1:rows_cnt*blocks_in_row
                   new_id =fix((id-1)/rows_cnt/cols_cnt) + mod(id-1,cols_cnt)+1;
                   if mod(id-1,rows_cnt*cols_cnt)==0 || mod(id-1,rows_cnt*cols_cnt)==1
                       non_interleaved_coes(:,row*blocks_in_row+new_id) = rows(:,id); 
                   else
                       non_interleaved_coes(:,(row+1)*blocks_in_row+new_id) = rows(:,id); 
                   end
               end
               row = row + rows_cnt;
               if row >= blocks_in_col
                   break
               end
           end         
        end     
        Re = [];
        Encoder.run_of_length = 0;
        for bandid = 1:blocks_in_component
            R = 0;
            ZZ = COEs(:,bandid);
            for K = 1:band_width
                if ZZ(K)
                    code= Encode_run_length(EHUFCO,EHUFSI);
                    Re =[Re code];
                    while R >=16
                        code = dec2bin(EHUFCO(240),EHUFSI(240));
                        Re =[Re code];
                        R = R - 16;
                    end
                    code = Encode_R_ZZ(R,ZZ(K));
                    Re =[Re code];
                    R = 0;
                    if K == band_width
                        break;
                    end
                else
                    R = R+1;
                    if K == band_width
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
        while length(Encoder.Re) >= 8
            byte = bin2dec(Encoder.Re(1:8));
            WriteOneByte(byte)
            if byte == 255
                WriteOneByte('00')
            end
            Encoder.Re = Encoder.Re(9:end);
        end
        if ~isempty(Encoder.Re) % �������һ���ֽڵ�Bits�������
            Encoder.Re = PrepareforMarker(Encoder.Re);
            WriteOneByte(bin2dec(Encoder.Re))
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
        function code = Encode_run_length(EHUFCO,EHUFSI)
            if Encoder.run_of_length == 0
                code = '';
                return
            else
                eobrun = floor(log2(Encoder.run_of_length));
                EOBn = bitshift(eobrun,4);
                % �������run_length��Ҫ��bit����Ȼ��תΪEOBn���б���
                huffmancode = EHUFCO(EOBn+1);
                huffmansize = EHUFSI(EOBn+1);
                huffmancode = dec2bin(huffmancode,huffmansize);
                code = [huffmancode dec2bin(Encoder.run_of_length-bitshift(1,eobrun),eobrun)];
                Encoder.run_of_length = 0;
            end
        end
    end
    function EncodeACRefine()

    end
end