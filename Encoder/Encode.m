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
SF = [2 1 1;2 1 1];
% 编码器结构体
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
SFw_max = max(Encoder.sample_factor(1,:));
SFh_max = max(Encoder.sample_factor(2,:));
Encoder.MCUs_in_row = ceil(Width/SFw_max/BLOCKSIZE);
Encoder.MCUs_in_col = ceil(Height/SFh_max/BLOCKSIZE);

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
    'qtbl_index',0,...             % quantization table selector (0..3)
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
    c.MCU_height = sfh;
    c.MCU_width = sfw;
    c.MCU_blocks = sfh*sfw;
    c.blocks_per_row = Encoder.MCUs_in_row * sfw;
    c.blocks_per_col = Encoder.MCUs_in_col * sfh;
    c.last_cols = c.blocks_per_row*BLOCKSIZE-Width;
    c.last_rows = c.blocks_per_col*BLOCKSIZE-Height;
    c.downsampled_height = c.blocks_per_col * BLOCKSIZE;
    c.downsampled_width = c.blocks_per_row * BLOCKSIZE;
    c.component = padarray(comp,[c.last_rows,c.last_cols],'symmetric','post');
    c.coes = zeros(DCTSIZE,c.blocks_per_row*c.blocks_per_col);
    Encoder.component(i) = c;
end
WriteOneByte = @(Bits) WriteNBytes(fid,Bits,1);
WriteTwoBytes = @(Bits) WriteNBytes(fid,Bits,2);
WriteFourBits = @(Bits) WriteNBytes(fid,Bits,0.5);

EncodeImage(); % 开始编码
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
        for i = 1: ImgInfo.components
            WriteOneByte(i);
            WriteFourBits(Encoder.sample_factor(1, i))
            WriteFourBits(Encoder.sample_factor(2, i))
            WriteOneByte(Encoder.component(i).qtbl_index)
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

    function  WriteSOS(channels)
        %{
            SOS数据段规定了每次编解码处理的
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
            这一段主要指定每个通道采用的Huffman表，存储格式为：两字节数据长度、
            一字节通道数、一字节通道序号、一字节Huffman表编号、一字节存储0，表
            示ZigZag编码时的开始序号、一字节存储63表示ZigZag编码的结束序号、一
            字节存储00，分别为Ah值和Al值，对基于DCT的编码过程，这两个值固定为0.
        %}
    end

    function  appendAPPn()
        % 写入为应用环境保留的信息，
        % 以Marker'FFE0-FFEF'开头
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
%         % 选择一种方式去存储系数，如果按照通道划分最好是用非交错模式存储
%         for chan = 1:ImgInfo.components
%             comp = Encoder.component(chan);
%             qtbl = Encoder.quanti_tbl{comp.qtbl_index};
%             for r  = 0:comp.blocks_per_col-1
%                 for c = 0:comp.blocks_per_row-1
%                     block_id = r*comp.blocks_per_col + c + 1;
%                     block = comp.component(r*BLOCKSIZE+1:(r+1)*BLOCKSIZE,...
%                         c*BLOCKSIZE+1:(c+1)*BLOCKSIZE);
%                     coef = round(dct2(block) ./ double(qtbl));
%                     coef = zigzag(coef);
%                     comp.coes(:,block_id) = coef;
%                 end
%             end
%             Encoder.component(chan) = comp;
%         end
        coef  = load('coe.mat').COEF;
        for i = 1:ImgInfo.components
            Encoder.component(i).coes = coef{i};
        end
%         save("coe.mat","COEF");
        EncodeDCFirst();
        EncodeACFirst(1);
        EncodeDCRefine();
    end
    function EncodeDCFirst()
        Encoder.Ss = 0;
        Encoder.Se = 0;
        Encoder.Ah = 0;
        Encoder.Al = 1;
        ResetEncoder();
        channels = ImgInfo.components;
        % 多个通道同时进行编解码时，需要按照交错模式存储
        HUFFVAL = zeros(1,blocks);
        DIFFs = zeros(1,blocks);
        for c_id = 1:channels
            comp = Encoder.component(c_id); 
            DCcoefs = comp.coes;
            DCcoefs = bitshift(DCcoefs(1,:),-Encoder.Al,'int16');
            blocks = comp.blocks_per_col*comp.blocks_per_row;
            for bandid = 1:blocks
                DIFF = DCcoefs(bandid) - Encoder.LastDCVal(c_id);
                codelength = EnsureGategory(DIFF);  % 确定编码DIFF需要的bit数
                HUFFVAL(bandid) = codelength;
                Encoder.LastDCVal(c_id) = DCcoefs(bandid);
                DIFFs(bandid) = DIFF;
            end
            if c_id == 1
                LuTbl = Generate_HuffTable(HUFFVAL);
            else

            end
            HuffVals{c_id} = HUFFVAL; % 所有通道的待编码值
            Actual_Val{c_id} = DIFFs;
        end
        
        ChrTbl = Generate_HuffTable([HuffVals{2:3}]);
        WriteDHT(0,LuTbl);
        WriteDHT(1,ChrTbl);
        WriteSOS(ImgInfo.components);
        for i = 1:ImgInfo.components
            if i == 1
                table = LuTbl;
            else
                table = ChrTbl;
            end
            comp = Encoder.component(i); 
            blocks = comp.blocks_per_col*comp.blocks_per_row;
            HUFFCODE = zeros(1,blocks);
            HUFFSIZE = zeros(1,blocks);
            for bandid = 1:blocks
                BITS = table{1};
                HUFFVAL = table{2};
                DIFF = Actual_Val{i}(bandid);
                codelength = HuffVals{i}(bandid);
                [code,length] = EncodeDC(DIFF,codelength,{BITS,HUFFVAL});
                HUFFCODE(bandid) = code;
                HUFFSIZE(bandid) = length;
            end
            HUFFCODEs{i} = HUFFCODE;
            HUFFSIZEs{i} = HUFFSIZE;
        end
        Re =[];
        bandid = [1,1,1];
        for ROW = 0:Encoder.MCUs_in_col-1
            for COL = 0:Encoder.MCUs_in_row-1
                for C = 1:ImgInfo.components
                    comp = Encoder.component(C);
                    htbl = Encoder.ta
                    for row = 0:comp.MCU_height-1
                        for col =0:comp.MCU_width-1
                            
                        end
                    end

                end
            end

        end
        for i = 1:sum(prod(Encoder.blocks_per_comp,1))
            temp = mod(i,Encoder.MCUs_per_iMCU);
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
        if ~isempty(Re) % 对最后不足一个字节的Bits进行填充
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
        % 编码完成DC系数之后，需要按照频率从低到高逐步开始编码AC系数，每次编码选
        % 定的频带，按照不同的通道分别进行编码
        while length(coefs) >= 8
            byte = bin2dec(join(coefs(1:8)));
            WriteOneByte(byte)
            %{
               每次存进去一个字节，遇到‘FF’需要在后面添加 '00'，防止
               和Marker冲突
            %}
            if byte == 255
                WriteOneByte('00')
            end
            coefs = coefs(9:end);
        end
        if ~isempty(coefs)    % 对最后不足一个字节的Bits进行填充
            coefs = PrepareforMarker(coefs);
            WriteOneByte(bin2dec(coefs))
        end
    end
    function EncodeACFirst(channelid)
        %{
            首先要根据AC系数编码思想确认待编码值，然后生成对应的huffman表
        %}
        Encoder.Ss = 1;
        Encoder.Se = 5;
        Encoder.Ah = 0;
        Encoder.Al = 2;
        ACcoefs = bitshift(Encoder.COEFICIENTS{channelid}(Encoder.Ss+1:Encoder.Se,:),-Encoder.Al,'int8');
        RSs = zeros(1,prod(Encoder.blocks_per_comp(:,channelid)));
        %{
            后面色度通道的HUFFMAN码表的生成问题需要完成，然后单通道按照非交错模式存储
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
            WriteDHT('10',Htable);
        else
            WriteDHT('11',Htable);
        end
        WriteSOS(1);
        [EHUFCO,EHUFSI] = SortHuffTbl(Htable);
        % 对所有的AC系数进行处理,将系数按照非交错模式重新排序一遍
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
        if ~isempty(Encoder.Re) % 对最后不足一个字节的Bits进行填充
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
                % 计算编码run_length需要的bit数，然后转为EOBn进行编码
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