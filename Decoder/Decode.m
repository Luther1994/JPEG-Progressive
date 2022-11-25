function Decode(imgfile,savefile)
jpegfile=fopen(imgfile,'rb');

% struct of decoder.
Decoder = struct( ...
    'LastDCVal',[0,0,0],...
    'quanti_tbl_idx',[0,0,0],...
    'sample_factor',[1;1],...
    'Ss',0, ...
    'Se',0, ...
    'Ah',0, ...
    'Al',0,...
    'EOBRUN',0,...
    'SOIFound',0, ...
    'bits_in_buffer',0,...
    'temp_buffer',0,...
    'next_marker',0, ...
    'JpegBuffer',fread(jpegfile), ...
    'FilePtr',1,...
    'savefile',savefile);

% struct of Image infomation for public use
ImgInfo = struct( ...
    'Precision',0,...
    'Channels',0,...
    'Height',0,...
    'Width',0);

% All markers that can be appeared
Markers = struct( ...
    'SOI',hex2dec('FFD8'),...
    'SOS',218,...
    'DQT',219,...
    'DHT',196,...
    'COM',254,...
    'EOI',217,...
    'DNL',220,...
    'SOFs',[192:195,197:203,205:207],...
    'APPn',224:249);


    function Byte = ReadOneByte()
        Byte = Decoder.JpegBuffer(Decoder.FilePtr);
        Decoder.FilePtr = Decoder.FilePtr + 1;
    end

    function Byte = ReadTwoBytes
        Byte = Decoder.JpegBuffer(Decoder.FilePtr:Decoder.FilePtr + 1);
        Byte = bitshift(Byte(1),8) + Byte(2);
        Decoder.FilePtr = Decoder.FilePtr + 2;
    end

    function Bytes = ReadNBytes(n,precision)
        if ~exist("precision",'var')
            precision = 1;
        end
        assert(precision == 1 | precision==2,'Data precision must be 1 or 2.')
        if precision==1
            Bytes = Decoder.JpegBuffer(Decoder.FilePtr:Decoder.FilePtr+n-1);
        else
            Bytes = Decoder.JpegBuffer(Decoder.FilePtr:Decoder.FilePtr+2*n-1);
            for i = 1:n
                Bytes(i) = bitshift(Bytes(2*i-1),8) + Bytes(2*i);
            end
            Bytes = Bytes(1:n);
        end
        Decoder.FilePtr = Decoder.FilePtr + precision * n;
    end
HUFF_LOOKAHEAD = 8;
JPEG_BLOCK_SIZE = 64;
while 1
    if ~Decoder.SOIFound
        New = ReadTwoBytes();
        if New == Markers.SOI
            Decoder.SOIFound = 1;
            fprintf('SOI Found,Decode Start!\n')
        elseif isempty(New)
            error('No SOI found,bad file format!');
        end
    else
        %{
            In this project,we try to read a marker 2 bytes once,
         so it must be guaranteed that the first 2 bytes in a 
         marker is 0xFF.
            While a marker appeared in image datasegment,assigning 
         it to Decoder.next_marker to execute corresponding opearation.
         IF no marker found,that is reading data with known data length,
         maker is next 4 bytes.
        %}
        Marker = Decoder.next_marker;
        if Marker == 0
            NewByte = ReadOneByte();
            if isempty(NewByte)
                % JPEG file terminate accidently.
                GenerateImage(Decoder,ImgInfo);
                fclose(jpegfile);
                error("Decoder.next_marker");
            end
            try
            assert(NewByte == 255,'Next Byte Must be 0xFF')
            catch
                return
            end
            Marker = ReadOneByte();
        end
        if ismember(Marker,Markers.SOFs)
            % 解码帧头信息
            if Marker == 192
                fprintf('BaseLine DCT Format.\n');
            elseif Marker == 193
                fprintf('Extend Sequential DCT Format.\n');
            elseif Marker == 194
                fprintf('Progressive DCT Format.\n');
            elseif Marker == 195
                fprintf('Lossless(sequential) Format.\n');
            end
            ParseSOF();
            Decoder.next_marker = 0;
        elseif ismember(Marker,Markers.APPn)
            ParseAPPn();
        else
            switch Marker
                case Markers.SOS
                    DecodeImage();
                case Markers.EOI
                    fprintf('EOI found,decode terminate.image saved as %s\n',savefile);
                    GenerateImage();
                    fclose(jpegfile);
                    break;
                case Markers.COM
                    % 添加的信息，单纯的图像解码可以不管
                    data_length = ReadTwoBytes();
                    NextNBytes(data_length - 2)
                    Decoder.next_marker = 0;
                case Markers.DHT
                    ParseDHT();
                case Markers.DQT
                    ParseDQT();
                case Markers.DNL
                    Decoder.next_marker = 0;
                    return
            end
        end
    end
end
    function ResetDecoder()
        Decoder.LastDCVal = [0,0,0];
        Decoder.bits_in_buffer = 0;
        Decoder.temp_buffer = 0;
        Decoder.comp_id = 1;
    end
    function DecodeImage()
        %{
            Main function of decode.
        %}
        ParseSOS();
        for i = 0:Decoder.MCUs-1
            Decoder.DecodeMCU(i);
        end
    end

    function DecodeDCFirst(mcu_cnt)
        %{
            First decode of DC coefficients
        %}
        for blkn = 1:Decoder.blks_in_MCU
            ci = Decoder.member_ship(blkn);
            Htable = Decoder.dc_huff_tbl(Decoder.dc_tbl_id(ci)+1);
            blk_cnt = Decoder.dc_blk_id(mcu_cnt*Decoder.blks_in_MCU+blkn);
            S = DECODE(Htable);
            DIFF = bitshift(EXTEND(get_bits(S), S),Decoder.Al,'int16');
            Decoder.LastDCVal(ci) = Decoder.LastDCVal(ci) + DIFF;
            Decoder.Coes(ci,1,blk_cnt) = Decoder.LastDCVal(ci);
        end
    end
    function DecodeACFirst(mcu_cnt)
        %   =================================
        %    First decode of AC coefficients
        %   =================================
        mcu_cnt = mcu_cnt + 1;
        HTable = Decoder.ac_huff_tbl(Decoder.ac_tbl_id(Decoder.comp_id)+1);

        % End of Band run length,if EOBRUN != 0,all left coefficient are 0
        if Decoder.EOBRUN > 0
            Decoder.EOBRUN = Decoder.EOBRUN - 1;
        else
            K = Decoder.Ss + 1;
            while K <= Decoder.Se+1
                RS = DECODE(HTable);
                % run_length & code size
                R = bitshift(RS,-4);
                S = bitand(RS,2^4-1);

                if S == 0
                    % code size == 0 means several continuous 0s.
                    if R == 15
                        % ZRL code
                        K = K + 16;
                    else
                        % note the method calucating EOURUN
                        Decoder.EOBRUN = bitshift(1,R) + get_bits(R)-1;
                        break
                    end
                else
                    % code size !=0 means non-zero coefficient
                    % appeared.
                    K = K + R;
                    Temp = EXTEND(get_bits(S),S);
                    Temp = bitshift(Temp,Decoder.Al,'int16');
                    Decoder.Coes(Decoder.comp_id,K,mcu_cnt) = Temp;
                    K = K + 1;
                end
            end
        end
    end

    
    function DecodeDCRefine(mcu_cnt)
        %   ================================
        %   Refine decode of  DC coefficient
        %   ================================
        P = 2^Decoder.Al;
        for blkn = 1:Decoder.blks_in_MCU
            ci = Decoder.member_ship(blkn);
            blk_cnt = Decoder.dc_blk_id(mcu_cnt*Decoder.blks_in_MCU+blkn);
            if get_bits(1)
                Decoder.Coes(ci,1,blk_cnt) = bitor( ...
                    Decoder.Coes(ci,1,blk_cnt),P,'int32');
            end
        end
    end

 function DecodeACRefine(mcu_cnt)
        %   ================================
        %   Refine decode of AC coefficients,
        %   One component decoded once.
        %   ================================
        mcu_cnt = mcu_cnt + 1;

        % End of Block
        EOBRUN = Decoder.EOBRUN;

        % index of component
        Id = Decoder.comp_id;

        HTable = Decoder.ac_huff_tbl(Decoder.ac_tbl_id(Id)+1);
        P = bitshift(1,Decoder.Al);

        K = 2;
        
        % 如果没有遇到end of band，则解码当前block
        if EOBRUN == 0
            while K <= 64
                RS = DECODE(HTable);
                R = fix(RS / 16);                  
                S = mod(RS , 16);
                if S
                    % 只要遇到ZZ（K）!= 0 & History(K)==0的情况就记为RS，同时
                    % 用一个bit来表示  
                    if(S ~= 1)
                        % size of new coef should always be 1
                        error('Corrupt JPEG data: bad Huffman code!');
                    end
                    if get_bits(1)
                        % ================================================
                        % Rule a. One additional bit codes the sign of 
                        % newly non-zero coef,with a 1-bits codes positive 
                        % sign and a 0-bit codes negative sign.
                        % =================================================
                        S = P;
                    else
                        S = -P;
                    end
                else
                    % S == 0 means sevral continuous 0s.  
                    if R ~= 15 
                        % End of Block,Calc run_length.
                        EOBRUN = 2 ^ R;
                        if R
                            EOBRUN = EOBRUN + get_bits(R);
                        end
                        break
                    end
                end
                while K <= 64
                    % the refine value is depended on sign of history val.
                    % ===================================================
                    % 注意这里R的计算没有考虑前一次解码结果不为0的系数，在统计
                    % R的时候跳过了前一精度下不为0的系数。
                    % ===================================================
                    ThisCoef = Decoder.Coes(Id,K,mcu_cnt);
                    if ThisCoef ~= 0
                        AssignCoef();
                    else
                        R = R - 1;
                        if R < 0
                            break
                        end
                    end
                    Decoder.Coes(Id,K,mcu_cnt) = ThisCoef;
                    K = K + 1;
                end
                if S
                    ThisCoef = S;
                    Decoder.Coes(Id,K,mcu_cnt) = ThisCoef;
                end
                K = K + 1;
            end
        end
        
        if EOBRUN > 0
            %  ============================================================
            %  Scan any remaining coefficient positions after the 
            %  end-of-band (the last newly nonzero coefficient, if any).  
            %  Append a correction bit to each already-nonzero coefficient. 
            %  A correction bit is 1 if the absolute value of the 
            %  coefficient must be increased.
            %  ============================================================
            while K <= 64
                ThisCoef = Decoder.Coes(Id,K,mcu_cnt);
                if ThisCoef ~= 0
                    AssignCoef();
                end
                Decoder.Coes(Id,K,mcu_cnt) = ThisCoef;
                K = K + 1;
            end
            EOBRUN = EOBRUN - 1;
        end
        Decoder.EOBRUN = EOBRUN;

        function AssignCoef()
            % Assign new coef to history existed coef.
            if get_bits(1)
                if bitand(abs(ThisCoef),P) == 0
                    if ThisCoef > 0
                        ThisCoef = ThisCoef + P;
                    else
                        ThisCoef = ThisCoef - P;
                    end
                end
            end
        end
    end

    function S = DECODE(Htable)
        %{
            Obtaining the value of symbol encoded by next bits from Huffman
            table

            Procedure which returns the 8-bit value associated with the next
            Huffman code in the compressed image data.
        %}

        if Decoder.bits_in_buffer < HUFF_LOOKAHEAD
            jpeg_fill_buffer();
        end
        buffer = Decoder.temp_buffer;
        bits_left = Decoder.bits_in_buffer;  % num of nused bits left in buffer 
        look = bitand(bitshift(buffer,-(bits_left-HUFF_LOOKAHEAD)),255);
        
        % Read one byte as the index of looktable,the high position of thre
        % result is the number of bits used to encode value,the 8 low
        % position is the encoding value.DIFF for DC coefficient and RS for
        % AC coefficient.
        temp = Htable.lookup(look+1);
        nb = bitshift(temp,-HUFF_LOOKAHEAD);
        if nb <= HUFF_LOOKAHEAD
            % if number of bits less than 8,lookup table. 
            S = bitand(temp,bitshift(1,HUFF_LOOKAHEAD) - 1);
            Decoder.bits_in_buffer = bits_left-nb;
            Decoder.temp_buffer = bitand(buffer,2^Decoder.bits_in_buffer-1);
        else
            % if number of bits bigger than 8,decode with tranditional
            % method.
            S = jpeg_huff_decode(Htable,nb);
        end
        
    end

    function ParseSOF()
        %{
            Start of Frame
        %}
        data_length = int32(ReadTwoBytes());
        data_length = data_length - 2;
        ImgInfo.Precision = ReadOneByte();
        ImgInfo.Height = double(ReadTwoBytes());
        ImgInfo.Width = double(ReadTwoBytes());
        ImgInfo.Channels = ReadOneByte();
        if (ImgInfo.Height <= 0 || ImgInfo.Width <= 0 ||...
                ImgInfo.Channels <= 0)
            error('"Empty JPEG image (DNL not supported)')
        end
        switch ImgInfo.Precision
            case 8
                type = 'uint8';
            case 16
                type = 'uint16';
        end
        ImgInfo.Img = zeros(ImgInfo.Height,ImgInfo.Width,ImgInfo.Channels,type);
        data_length = data_length - 6;
        fprintf('Image precision is %d bit,with shape of %d x %d x %d.\n',...
            ImgInfo.Precision,ImgInfo.Height,ImgInfo.Width,ImgInfo.Channels);
        
        if data_length ~= (ImgInfo.Channels * 3)
            error("Bogus marker length");
        end
        Decoder.init_blk_id = ones(1,ImgInfo.Channels);
        Idx = 1;
        for i = 1:ImgInfo.Channels
            % 通道编号和采样系数
            c_id = ReadOneByte();
            sample_factor = ReadOneByte();
            sf_w = bitshift(sample_factor,-4);
            sf_h = bitand(sample_factor,2^4-1);
            Decoder.sample_factor(:,c_id)=[sf_w,sf_h];

            % 量化表编号
            Decoder.quanti_tbl_idx(c_id) = ReadOneByte() + 1;
            
            % 各通道信息
            c1 = sf_h*sf_w;

            for j = 1:c1
                Decoder.member_ship(Idx) = i;
                Idx = Idx + 1;
            end
            data_length = data_length - 3;
        end
        
        sf_h = max(Decoder.sample_factor(1,:));
        sf_w = max(Decoder.sample_factor(2,:));

        DownSam = [sf_h;sf_w]./Decoder.sample_factor;
        Decoder.blks_per_row = ceil(ImgInfo.Width./DownSam(1,:)/8);
        Decoder.blks_per_col = ceil(ImgInfo.Height./DownSam(2,:)/8);
        Decoder.blks_per_comp = Decoder.blks_per_col.*Decoder.blks_per_row;
        Decoder.MCUs_per_row = ceil(ImgInfo.Width/8/sf_w);
        Decoder.MCUs_per_col = ceil(ImgInfo.Height/8/sf_h);

        Decoder.blks_in_MCU = sum(prod(Decoder.sample_factor));
        Decoder.Coes = zeros(ImgInfo.Channels,JPEG_BLOCK_SIZE,max(Decoder.blks_per_comp));
        
        % 这里保存dc系数解码时的block排列顺序，用于应对多通道编码的交错存储模式
        Decoder.dc_blk_id = ones(1,sum(Decoder.blks_per_comp));

        blkn = 1;
        for ROW = 0:Decoder.MCUs_per_col - 1
            for COL = 0:Decoder.MCUs_per_row - 1
                for ci = 1:ImgInfo.Channels
                    McuWidth = Decoder.sample_factor(1,ci);
                    McuHeight = Decoder.sample_factor(2,ci);
                    for r = 0:McuHeight - 1
                        for c = 0:McuWidth - 1
                            row = ROW * McuHeight + r;
                            col = COL * McuWidth + c;
                            blk_id = row * Decoder.blks_per_row(ci) + col + 1;
                            Decoder.dc_blk_id(blkn) = blk_id;
                            blkn = blkn + 1;
                        end
                    end
                end
            end
        end
        
        assert (data_length == 0,'Bad Data Length of SOF.')
        Decoder.next_marker = 0;
    end

    function ParseAPPn()
        %{
            Parse APPn data segment,this is about some comment infomation
        %}
        data_length = int32(ReadTwoBytes());
        DataSegment = ReadNBytes(data_length - 2);
        file_format = char();

        % file save format,usually JFIF
        if Marker == hex2dec('E0')
            for Idx =1:data_length - 2
                if (32 <= DataSegment(Idx))...
                        && (DataSegment(Idx)<= 126)
                    format_identifier = char(DataSegment(Idx));
                    file_format = [file_format format_identifier]; %#ok<AGROW>
                end
            end
        end
        Decoder.next_marker = 0;
        if file_format~='.'
            fprintf('File is of format %s.\n',file_format);
        end
    end

    function ParseDQT()
        %{
            Function parsing quantization tables.Quantization tbales
            consists of 64 elements with zigzag order,each of them is 
            storaged by 2 bytes.
        %}
        data_length=int32(ReadTwoBytes());
        data_length = data_length-2;
        while data_length > 0
            c = ReadOneByte();

            % precision of element,8/16 bits
            dqt_precision = bitshift(c,-4); 

            % identifier of table,0/1
            identifier = bitand(c,2^4-1);
            data_length = data_length -1;
            if dqt_precision == 0
                Qtbl = ReadNBytes(64);
            else
                Qtbl = ReadNBytes(64,2);
            end
            data_length = data_length - 64;

            % quantization table should be saved as zigzag order
            Decoder.dqt_ids{identifier+1} = inverse_zigzag(Qtbl);
        end
        Decoder.next_marker = 0;
    end

    function ParseSOS()
        %{
            SOS data segment,which includes channel information to be
            decoded.Number of channels to be decoded、which huffman table
            used by which channel、band specified for progressive
            mode、precision specified for approximation.
        %}

        % Initialization of channel infomation
        ResetDecoder();

        data_length = ReadTwoBytes();
        channels = ReadOneByte();

        Decoder.channels = channels;
        data_length = data_length - 3;

        for i = 1:channels
            c = ReadOneByte();
            Decoder.comp_id(i) = c;

            id = ReadOneByte();
            Decoder.dc_tbl_id(c) = bitshift(id,-4);
            Decoder.ac_tbl_id(c) = bitand(id,2^4-1);
            data_length = data_length - 2;
        end
        
        if channels == 1
            Decoder.MCUs = Decoder.blks_per_comp(Decoder.comp_id);
        else
            Decoder.MCUs = Decoder.MCUs_per_col*Decoder.MCUs_per_row;
        end
        % start and end band in progressive mode
        Decoder.Ss = ReadOneByte();
        Decoder.Se = ReadOneByte();
        if Decoder.Ss > Decoder.Se || Decoder.Se >= JPEG_BLOCK_SIZE
            error('Bad Ss=%d / Se=%d',Decoder.Ss,Decoder.Se)
        end

        % the highest/lowest bit of approximation mode
        NewByte = ReadOneByte();
        Decoder.Ah = bitshift(NewByte,-4);
        Decoder.Al = bitand(NewByte,15);
        
        if (Decoder.Ah ~= 0)
    		% Successive approximation refinement scan: must have Al = Ah-1
            if Decoder.Al ~= Decoder.Ah - 1
                error('Bad Ah(%d) / Al(%d)',Decoder.Ah,Decoder.Al)
            end
        end
        
        data_length = data_length - 3;
        
        % specify docode function vias readed parameters
        if Decoder.Ss == 0
            if Decoder.Ah == 0
                Decoder.DecodeMCU = @DecodeDCFirst;
            else
                Decoder.DecodeMCU = @DecodeDCRefine;
            end
        else
            if Decoder.Ah == 0
                Decoder.DecodeMCU = @DecodeACFirst;
            else
                Decoder.DecodeMCU = @DecodeACRefine;
            end
        end
        assert(data_length==0,"Bad Data Length of SOS.")
        Decoder.next_marker = 0;
    end

    function [V] = EXTEND(V, T)
        %{
            Extending the sign bit of a decoded value in V

            Procedure which converts the partially decoded DIFF value of precision T to the full
            precision difference by correcting sign bit in MSB.

            Args
                V Value which is unsigned of precision T
                T Code length with which V is encoded

            Returns
                V Corrected signed value

        %}
        if V < power(2,T-1)
            V = V - power(2,T)+1;
        end
    end

    function Byte = GetNewByte()
        Byte = 0;
        if ~Decoder.next_marker
            Byte = ReadOneByte();
            if Byte == hex2dec('FF')
                NextByte = ReadOneByte();
                while NextByte == hex2dec('FF')
                    NextByte = ReadOneByte();
                end
                if NextByte == 0
                    return;
                elseif ismember(NextByte,[196,217,218])
                    Decoder.next_marker = NextByte;
                else
                    error('Illegal marker!')
                end
            end
        end
    end
    
    function ParseDHT()
        %{
            Parse huffman table.
        %}
        bits = zeros(1,17);
        huffval = zeros(1,256);
        data_length = int32(ReadTwoBytes());
        data_length = data_length - 2;
        while data_length > 16
            index = ReadOneByte();
            isAC = bitshift(index,-4);
            index = bitand(index,2^4-1);
            bits(1) = 0;
            count = 0;
            for i = 2:17
                bits(i) = ReadOneByte();
                count = count +  bits(i);
            end
            data_length = data_length - 17;
            if count > 256 || count > data_length
                error('Bogus Huffman table definition')
            end
            for i = 1:count
                huffval(i) = ReadOneByte();
            end
            data_length =data_length - count;
            table.bits = bits;
            table.huffval = huffval;
            if isAC         % AC table definition
                Decoder.ac_huff_tbl(index+1)=jpeg_make_derived_tbl(0,table);
            else            % DC table definition
                Decoder.dc_huff_tbl(index+1)=jpeg_make_derived_tbl(1,table);
            end
        end
        assert(data_length == 0,('Counts of bits left must be 0.'))
        Decoder.next_marker = 0;
    end

    function dtbl = jpeg_make_derived_tbl(isDC, htbl)
        %{
            Note that huffsize() and huffcode() are filled in code-length order,
            paralleling the order of the symbols themselves in htbl.huffval().
        %}
        dtbl.pub = htbl;
        if isempty(htbl)
            MException('MATLAB:badargs','Huffman table is empty.');
            throw(MException)
        end

        % Figure C.1: make table of Huffman code length for each symbol

        index = 0;
        huffsize = zeros(1,sum(htbl.bits)+1);
        huffcode = zeros(1,sum(htbl.bits));
        for l = 2:17
            cnt = htbl.bits(l);
            if cnt < 0 || index + cnt > 256   % protect against table overrun
                ERREXIT(cinfo, JERR_BAD_HUFF_TABLE);
            end
            while cnt
                cnt = cnt-1;
                huffsize(index+1) = l-1;
                index = index+1;
            end
        end
        huffsize(index+1) = 0;
        numsymbols = index;

        % Figure C.2: generate the codes themselves
        % We also validate that the counts represent a legal Huffman code
        % tree.

        code = 0;
        l = huffsize(1);
        index = 1;
        while huffsize(index)
            while huffsize(index) == l
                huffcode(index) = code;
                index = index+1;
                code= code + 1;
            end
            if code >= bitshift(1,l)
                ERREXIT(cinfo, JERR_BAD_HUFF_TABLE);
            end
            code =bitshift(code,1);
            l = l + 1;
        end
        %{
         code is now 1 more than the last code used for codelength si; but
         it must still fit in si bits, since no code is allowed to be all
         ones.
         Figure F.15: generate decoding tables for bit-sequential decoding
        %}
        index = 1;
        for l = 2:17
            % num of code with length == l-1(Matlab starts with 1)
            num =htbl.bits(l); 
            if num
                %{
                    valoffset(l) = huffval() index of 1st symbol of code length l,
                    minus the minimum code of length l
                %}
                min_code = huffcode(index); 
                dtbl.valoffset(l) = index - min_code;

                % index of 1st symbol of code with next length l
                index = index + num; 

                % maximum code of length l
                dtbl.maxcode(l) = huffcode(index - 1); 
            else
                % -1 if no codes of this length
                dtbl.maxcode(l) = -1;    
            end
        end
        % avoid the index overflow array length
        % ensures jpeg_huff_decode terminates 
        dtbl.valoffset(18) = 0; 
        dtbl.maxcode(18)= hex2dec('FFFFF'); 
        %{
            Compute lookahead tables to speed up decoding.
            First we set all the table entries to 0, indicating "too long";
            then we iterate through the Huffman codes that are short enough and
            fill in all the entries that correspond to bit sequences starting
            with that code.
        %}
        for cnt=1: bitshift(1 ,HUFF_LOOKAHEAD)
            dtbl.lookup(cnt) = bitshift((HUFF_LOOKAHEAD + 1),HUFF_LOOKAHEAD);
        end
        index = 1;
        for L = 1:HUFF_LOOKAHEAD
            for cnt = 1:htbl.bits(L+1)
                % l = current code's length, p = its index in huffcode() & huffval()
                % Generate left-justified code followed by all possible bit sequences
                lookbits = bitshift(huffcode(index) ,(HUFF_LOOKAHEAD - L));
                for ctr = bitshift(1 , (HUFF_LOOKAHEAD - L)):-1:1
                    temp = bitshift(L,HUFF_LOOKAHEAD);
                    dtbl.lookup(lookbits+1) = bitor(temp,htbl.huffval(index));
                    lookbits=lookbits+1;
                end
                index = index + 1;
            end
        end
        %{
             Validate symbols as being reasonable.
             For AC tables, we make no check, but accept all byte values 0..255.
             For DC tables, we require the symbols to be in range 0..15.
             (Tighter bounds could be applied depending on the data depth and mode,
             but this is sufficient to ensure safe decoding.)
        %}
        if isDC
            for cnt = 1: numsymbols
                sym = htbl.huffval(cnt);
                if sym < 0 || sym > 15
                    error('MATLAB:badargs','bad huffman code');
                end
            end
        end
    end

    function re = jpeg_huff_decode(htbl,l)
        % HUFF_DECODE has determined that the code is at least min_bits 
        % bits long, so fetch that many bits in one swoop. 
        code = get_bits(l);

        % Collect the rest of the Huffman code one bit at a time. 
        while code > htbl.maxcode(l+1)
            code = bitor(bitshift(code,1),get_bits(1));
            l = l + 1;
        end

        % With garbage input we may reach the sentinel value l = 17. 
        if (l > 16)
            error("Corrupt JPEG data: bad Huffman code");
        end
        re =  htbl.pub.huffval(code + htbl.valoffset(l+1));
    end

    function jpeg_fill_buffer()
        buffer = Decoder.temp_buffer;
        nbits = Decoder.bits_in_buffer;
        while nbits <= 24   % number of bits read-in each time less than 32
            NewByte = GetNewByte();
            buffer = bitor(bitshift(buffer,HUFF_LOOKAHEAD),NewByte);
            nbits = nbits + 8;
        end
        Decoder.temp_buffer=buffer;
        Decoder.bits_in_buffer=nbits;
    end

    function re = get_bits(n)
        % read-in n bits data
        if Decoder.bits_in_buffer < n
            jpeg_fill_buffer()
        end
        buffer = Decoder.temp_buffer;
        nbits = Decoder.bits_in_buffer;
        nbits = nbits-n;
        re = bitshift(buffer,-nbits);
        buffer = bitand(buffer,2^nbits-1);
        Decoder.temp_buffer = buffer;
        Decoder.bits_in_buffer = nbits;
    end

     function GenerateImage()
        %{
            enerate and Display Image from Decoded Coefficients.
        %}
        VMax = max(Decoder.sample_factor(1,:));
        HMax = max(Decoder.sample_factor(2,:));
        SF = [VMax./Decoder.sample_factor(1,:);HMax./Decoder.sample_factor(2,:)];
        
        blk_cnt = 1;
        for ROW = 0:Decoder.MCUs_per_col-1
            for COL = 0:Decoder.MCUs_per_row-1
                MCU = zeros(8*VMax,8*HMax,ImgInfo.Channels);
                for ci = 1:ImgInfo.Channels
                    McuInComp = zeros(8*Decoder.sample_factor(1,ci),8*Decoder.sample_factor(2,ci));
                    for r = 0:Decoder.sample_factor(1,ci) - 1
                        for c = 0:Decoder.sample_factor(2,ci) - 1
                            index = Decoder.dc_blk_id(blk_cnt);
                            DQTable = Decoder.dqt_ids{Decoder.quanti_tbl_idx(ci)};
                            block = Decoder.Coes(ci,:,index);
                            temp = inverse_zigzag(block);
                            temp = idct2(temp .* DQTable)+ 2 ^ (ImgInfo.Precision - 1);
                            McuInComp(r*8+1:(r+1)*8,c*8+1:(c+1)*8) = temp;
                            blk_cnt = blk_cnt + 1;
                        end
                    end
                    MCU(:,:,ci) = UPSAMPLE(McuInComp,SF(:,ci));
                end
                ImgInfo.Img(ROW*8*VMax+1:(ROW+1)*8*VMax,COL*8*HMax+1:(COL+1)*8*HMax,:) = MCU;
            end
        end
        if ImgInfo.Channels == 3
            ImgInfo.Img = yuv2rgb(ImgInfo.Img);
        end
        imshow(ImgInfo.Img);
     end
end


