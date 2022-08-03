function Encode(src,savefile,quality)
%{
    编码主函数
%}
fid = fopen(savefile, 'wb');
BLOCKSIZE = 8;
DCTSIZE = 64;
if max(max(src)) > 255
    precision  = 16;
else
    precision = 8;
end

% 编码器结构体
Encoder = struct('QLUTbl',QuantiTbl(quality).Q_LU,...
    'QCHRTbl',QuantiTbl(quality).Q_CHR,...
    'QtblIdx',ones(1,3),...
    'SampleFactors',[2 1 1;2 1 1],...
    'Ss',0,'Se',0,'Ah',0,'Al',0,...
    'run_length',0,...
    'RowIdx',0,'ColIdx',0,...
    'bitsleft',0,...
    'EOIReached',false,...
    'TemporaryBuffer',0 ,...
    'SOIFound',0, ...
    'savefile',savefile, ... 
    'EOB_RL',0,'Re',char());
[Hmax,Vmax] = deal(max(Encoder.SampleFactors(1,:)),...
    max(Encoder.SampleFactors(2,:)));
SFs = [Hmax./Encoder.SampleFactors(1,:);Vmax./Encoder.SampleFactors(2,:)];
Encoder.BlocksInMCU = [Encoder.SampleFactors(1,:);Encoder.SampleFactors(2,:)];
Encoder.MCUs_per_iMCU = sum(prod(Encoder.BlocksInMCU));
assert (ndims(src)>=2,'Dim of input image must ≥ 2.')
[Height,Width,CHANNELS] = size(src);
Mcu_height = max(BLOCKSIZE*SFs(1,:));
Mcu_width = max(BLOCKSIZE*SFs(2,:));
if mod(Height,Mcu_height)
    addition_row = Mcu_height-mod(Height,Mcu_height);
    Height = Height + addition_row;
    src = vertcat(src,flipud(src(end-addition_row+1:end,:,:)));
end
if  mod(Width,Mcu_width)
    addition_col = Mcu_width-mod(Width,Mcu_width);
    Width = Width + addition_col;
    src = horzcat(src,fliplr(src(:,end-addition_col+1:end,:)));
end
if CHANNELS == 3
    src = int16(rgb2yuv(src))-2^(precision-1);
else
    src = int16(src)  - 2^(precision -1);
end

% 图像信息
ImgInfo = struct('Precision',precision,...
    'Channels',CHANNELS,...
    'Height',Height,...
    'Width',Width);

jpeg_component_info = struct(...
    'component_id',0,...             % identifier for this component (0..255)
    'component_index',0,...          % its index in SOF or cinfo->comp_info[]
    'h_samp_factor',0,...            % horizontal sampling factor (1..4)
    'v_samp_factor',0,...            % vertical sampling factor (1..4)
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

% 所有可能出现的Markers
Markers = struct('SOI',65496,...
    'SOS',65498,...
    'DQT',65499,...
    'DHT',65476,...
    'COM',65534,...
    'EOI',65497,...
    'DNL',65500,...
    'SOFs',[65472:65475,65477:65483,65485:65587],...
    'APPn',65504:65519);


% 实际的采样与给定的系数相反，即[2 1 1]表示Y通道不采样，U/V通道2采1.

Encoder.ComponentSize = [ImgInfo.Height./SFs(1,:);...
    ImgInfo.Width./SFs(2,:)];
Encoder.blocks_per_comp = uint16(Encoder.ComponentSize./BLOCKSIZE);
Encoder.Components = {DOWNSAMPLE(src(:,:,1),SFs(:,1)),...
    DOWNSAMPLE(src(:,:,2),SFs(:,2)),...
    DOWNSAMPLE(src(:,:,3),SFs(:,2))};
Encoder.DUcnts = sum(prod(Encoder.blocks_per_comp,1));

WriteOneByte = @(Bits) WriteNBytes(fid,Bits,1);
WriteTwoBytes = @(Bits) WriteNBytes(fid,Bits,2);
WriteFourBits = @(Bits) WriteNBytes(fid,Bits,0.5);

% 量化表精度，0表示8位，1表示16位
dqt_precision = QuantiTbl(quality).precision;

%{
    通过码长和值，得到huffman码字表，这个表从按照码长排序改为按照值的大小排序。
%}

EncodeImage(); % 开始编码
    function  appendSOI()
        % 写入图像开始标记  FFD8
        WriteTwoBytes(Markers.SOI);
    end

    function  appendEOI()
        % 写入图像结束标记  FFD9
        WriteTwoBytes(Markers.EOI);
    end

    function  appendDQT(identifier)
        % 写入量化表
        WriteTwoBytes(Markers.DQT);
        if identifier
            qtb = zigzag(Encoder.QCHRTbl);
        else
            qtb = zigzag(Encoder.QLUTbl);
        end
        if dqt_precision == 8
            qtb = uint8(qtb);
        end

        % Defining type of qtb to uint16 is to cope with the situation
        % that quality too small so that quantification steps are too big

        datalength = length(qtb) + 2 + 1;
        WriteTwoBytes(datalength);
        if dqt_precision == 8
            WriteFourBits(0);
        else
            WriteFourBits(1);
        end
        WriteFourBits(identifier);
        WriteOneByte(qtb);
        %{
            量化表数据段中， Marker后面两个字节存储数据长度，
            一个字节存储量化表编号，量化表每一个元素用一个字节存储
        %}
    end

    function  appendSOF()
        %{
            根据不同的编码方式写入不同的帧头
        %}
        WriteTwoBytes(Markers.SOFs(3))
        datalength = 8 + ImgInfo.Channels * 3;
        WriteTwoBytes(datalength)
        WriteOneByte(precision)
        WriteTwoBytes(ImgInfo.Height)
        WriteTwoBytes(ImgInfo.Width)
        WriteOneByte(ImgInfo.Channels)
        quan_table_id = [0, 1, 1];
        for index = 1: ImgInfo.Channels
            WriteOneByte(index);
            WriteFourBits(Encoder.SampleFactors(1, index))
            WriteFourBits(Encoder.SampleFactors(2, index))
            WriteOneByte(quan_table_id(index))
        end
        %{
            这部分包含图像基本信息。两个字节存储数据长度、一个字节存储图像位数、
            两个字节存储图像高度、两个字节存储图像宽度、一个字节存储图像通道数，
            后面存储的每个通道的采样因子和指定的量化表，按照一字节序号、半字节水
            平方向采样因子、半字节竖直方向采样因子、一字节量化表编号的格式存储，
            直到存完所有的通道信息。
        
        %}
    end

    function  appendDHT(id,table)
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

    function  appendSOS(ids)
        CHANNELS = length(ids);
        WriteTwoBytes(Markers.SOS)
        datalength = 6 + 2 * CHANNELS;
        WriteTwoBytes(datalength)
        WriteOneByte(CHANNELS)
        for i = 1:length(ids)
            WriteOneByte(i)
            WriteOneByte(ids(i))
        end
        WriteOneByte(Encoder.Ss);
        WriteOneByte(Encoder.Se);
        WriteFourBits(Encoder.Ah);
        WriteFourBits(Encoder.Al); % Ss = 0, Se = 63
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
        % 初始化DC系数的预测值，也就是前一个block的DC系数
        Encoder.LastDCVal = zeros(ImgInfo.Channels,1);
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
%         for chan = 1:ImgInfo.Channels
%             blocks_in_ver = Encoder.blocks_per_comp(1,chan);
%             blocks_in_hor = Encoder.blocks_per_comp(2,chan);
%             for row = 0:blocks_in_ver-1
%                 for col = 0:blocks_in_hor-1
%                     block = Encoder.Components{chan}(row * BLOCKSIZE + 1 :(row + 1) * BLOCKSIZE,...
%                         col * BLOCKSIZE + 1:(col + 1) * BLOCKSIZE);
%                     coef = round(dct2(block) ./ double(Encoder.QLUTbl));
%                     coef = zigzag(coef);
%                     Encoder.component_coes{chan}(row+1,col+1,:) = coef;
%                 end
%             end
%         end
%         coef = Encoder.component_coes;
%         save("coe.mat","coef");
        Encoder.COEFICIENTS = load('coe.mat').coef;
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
        % 多个通道同时进行编解码时，需要按照交错模式存储
        for i = 1:ImgInfo.Channels
            DCcoefs = Encoder.COEFICIENTS{i}(:,:,1);
            
            DCcoefs = bitshift(DCcoefs(:),-Encoder.Al,'int16');
            blocks = prod(Encoder.blocks_per_comp(:,i)); 
            HUFFVAL = zeros(1,blocks);
            DIFFs = zeros(1,blocks);
            for bandid = 1:blocks
                DIFF = DCcoefs(bandid) - Encoder.LastDCVal(i);
                codelength = EnsureGategory(DIFF);  % 确定编码DIFF需要的bit数
                HUFFVAL(bandid) = codelength;
                Encoder.LastDCVal(i) = DCcoefs(bandid);
                DIFFs(bandid) = DIFF;
            end
            HuffVals{i} = HUFFVAL; % 所有通道的待编码值
            Actual_Val{i} = DIFFs;
        end
        LuTbl = Generate_HuffTable(HuffVals{1});
        ChrTbl = Generate_HuffTable([HuffVals{2:3}]);
        appendDHT('00',LuTbl);
        appendDHT('01',ChrTbl);
        tableid = [0,17,17];
        appendSOS(tableid);
        for i = 1:ImgInfo.Channels
            if i == 1
                table = LuTbl;
            else
                table = ChrTbl;
            end
            blocks = prod(Encoder.blocks_per_comp(:,i));
            HUFFCODE = zeros(1,blocks);
            HUFFSIZE = zeros(1,blocks);
            for bandid = 1:blocks
                BITS = table{1};
                HUFFVAL = table{2};
                DIFF = Actual_Val{i}(bandid);
                codelength = HuffVals{i}(bandid);
                [code,size] = EncodeDC(DIFF,codelength,{BITS,HUFFVAL});
                HUFFCODE(bandid) = code;
                HUFFSIZE(bandid) = size;
            end
            HUFFCODEs{i} = HUFFCODE;
            HUFFSIZEs{i} = HUFFSIZE;
        end
        Re =[];
        bandid = [1,1,1];
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
            size = HUFFSIZEs{channelid}(bandid(channelid));
            Re = [Re dec2bin(code,size)];
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
        appendSOS(tableid);
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
            appendDHT('10',Htable);
        else
            appendDHT('11',Htable);
        end
        appendSOS(1);
        [EHUFCO,EHUFSI] = SortHuffTbl(Htable);
        % 对所有的AC系数进行处理,将系数按照非交错模式重新排序一遍
        band_width = Encoder.Se-Encoder.Ss;
        blocks_in_col = prod(Encoder.blocks_per_comp(1,channelid));
        blocks_in_row = prod(Encoder.blocks_per_comp(2,channelid));
        blocks_in_component = prod(Encoder.blocks_per_comp(:,channelid));
        COEs = bitshift(Encoder.COEFICIENTS{channelid}(Encoder.Ss+1:Encoder.Se,:),-Encoder.Al,'int8');
        non_interleaved_coes = zeros(size(COEs));
        row = 0; 
        rows_cnt = Encoder.SampleFactors(2,channelid);
        cols_cnt = Encoder.SampleFactors(1,channelid);
        if Encoder.SampleFactors(2,channelid) > 1
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