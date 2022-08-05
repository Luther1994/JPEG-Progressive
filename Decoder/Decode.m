function Decode(imgfile,savefile)
jpegfile=fopen(imgfile,'rb');
% ������
Decoder = struct( ...
    'LastDCVal',[0,0,0],...
    'quanti_tbl_idx',[0,0,0],...
    'sample_factor',[1;1],...
    'Ss',0, ...
    'Se',0, ...
    'Ah',0, ...
    'Al',0,...
    'EOBRUN',0,...
    'RowIdx',0, ...
    'ColIdx',0,...
    'bitsleft',0,...
    'EOIReached',false,...
    'SOIFound',0, ...
    'savefile',savefile, ...
    'SOScnt',0, ...
    'DHTcnt',0,...
    'bits_in_buffer',0,...
    'temp_buffer',0,...
    'marker_found',0,...
    'next_marker',0);

% ͼ����Ϣ
ImgInfo = struct( ...
    'Precision',0,...
    'Channels',0,...
    'Height',0,...
    'Width',0);

% ���п��ܳ��ֵ�Markers
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

ReadOneByte = @() ReadNBytes(jpegfile,1);
ReadTwoBytes = @() ReadNBytes(jpegfile,2);
ReadFourBits = @() ReadNBytes(jpegfile,0.5);

% ��������ѭ������marker�����ݶ����markerִ����Ӧ����
HUFF_LOOKAHEAD = 8;
JPEG_BLOCK_SIZE = 64;
while 1
    if ~Decoder.SOIFound
        New = ReadTwoBytes();
        if New == Markers.SOI
            % ����ffd8,��ʼ����
            Decoder.SOIFound = 1;
            fprintf('SOI Found,Decode Start!\n')
        elseif isempty(New)
            % ����ļ��𻵣��Ҳ���ffd8�������ļ�����
            error('No SOI found,bad file format!');
        end
    else
        % ȷ��ÿ��ͼƬ��Ϣ���ֽ�����ɺ���뵽��marker������FF��ͷ
        Marker = Decoder.next_marker;
        if Marker == 0
            NewByte = ReadOneByte();
            if isempty(NewByte)
                GenerateImage(Decoder,ImgInfo);
                fclose(jpegfile);
                break;
            end
            assert(NewByte == 255,'Next Byte Must be 0xFF')
            Marker = ReadOneByte();
        end
        if ismember(Marker,Markers.SOFs)
            % ����֡ͷ��Ϣ
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
                    Decoder.SOScnt = Decoder.SOScnt+1;
                    DecodeImage();
                case Markers.EOI
                    fprintf('EOI found,decode terminate.image saved as %s\n',savefile);
                    GenerateImage(Decoder,ImgInfo);
                    fclose(jpegfile);
                    break;
                case Markers.COM
                    % ��ӵ���Ϣ��������ͼ�������Բ���
                    SegLength = ReadTwoBytes();
                    NextNBytes(SegLength - 2)
                    Decoder.next_marker = 0;
                case Markers.DHT
                    Decoder.DHTcnt = Decoder.DHTcnt+1;
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
        Decoder.Channels = [];
        Decoder.LastDCVal = zeros(1,3);
        Decoder.block_idx = 0;
        Decoder.RowIdx = 0;
        Decoder.ColIdx = 0;
        Decoder.bits_in_buffer = 0;
        Decoder.temp_buffer = 0;
    end
    function DecodeImage()
        %{
            ��Ҫ���뺯����ÿ������FFDA�ؼ���ʱ����Ҫ����ָ��ͨ��ָ��ϵ���Ľ��롣
        %}

        % ����ɨ��ͷStart of Scan
        ParseSOS();
        % ѭ��������DCϵ��һ�ν���һ��MCU��ѭ��������MCU������ͬ��
        % ACϵ����ͨ������������ͨ��ˮƽ����һ�ν������С���ֱ����һ�ν���һ�У�
        % ɫ��ͨ��ˮƽ�������ֵ����ֱ�������У��У���ϵ���Ľ���ʹ洢��ʽҲ
        % ���۽�����ͻ��߱����һ����Ҫ��ͬ�㣬���ǽ���洢�ͽ���洢��
        % ѭ��������ָ��ͨ����ָ��ϵ��
        c1 = Decoder.channel_info.MCUs_per_row;
        c2 = Decoder.channel_info.MCUs_per_col;
        for r = 0:c1-1
            for c = 0:c2-1
                mcu_id = r*c2 + c;
                Decoder.DecodeMCU(mcu_id);
            end
        end
        %���һ�ν����Ҫ�������л�δ�����bits������
        Decoder.bitsleft = 0;
    end
    function DecodeDCFirst(id)
        %{
            DCϵ����һ�ν��룬���λ
        %}
        Width = Decoder.channel_info.MCUs_per_row;
        for i = 1:Decoder.channel_info.channels
            c = Decoder.channel_info.id(i); 
            MCU_width = Decoder.sample_factor(1,c);
            MCU_height = Decoder.sample_factor(2,c);
            row =  MCU_height * fix(id/Width);
            col = MCU_width * mod(id,Width);
            Htable = Decoder.dc_huff_tbl(Decoder.channel_info.dc_tbl_id(c)+1);
            for j = 1:MCU_height
                for k = 1:MCU_width
                    blk_idx = row*Width*MCU_width + col + k;
                    % Huffman������̣������JPEG�ĵ�
                    DIFF = DECODE(Htable);
                    % DCϵ��ͨ����ִ洢
                    Decoder.LastDCVal(c) = Decoder.LastDCVal(c) + DIFF;
                    % ������д�뻺����
                    Decoder.Coes{c}(1,blk_idx)= Decoder.LastDCVal(c);
                end
                row = row + 1;
            end
        end
    end
    function DecodeACFirst(mcu_cnt)
        %{
            ACϵ����һ�ν��룬���λ
        %}
        % EOB run length,�����ǰƵ���Ժ��Ƶ����Ϊ0������EOBRUN��ʾ
        EOBRUN = Decoder.EOBRUN;
        c = Decoder.channel_info.id;
        tbl_id = Decoder.channel_info.ac_tbl_id+1;
        HTable = Decoder.ac_huff_tbl(tbl_id);
        blks = Decoder.channel_info.blks_in_MCU;
        if EOBRUN > 0
            EOBRUN = EOBRUN - 1;
        else
            c1 = Decoder.channel_info.MCU_width;
            c2 = Decoder.channel_info.MCU_height;
            for i = 1:c1
                for j = 1:c2
                    temp = (i-1)*c2 + j;
                    blk_cnt = mcu_cnt*blks + temp;
                    K = Decoder.Ss+1;
                    while K <= Decoder.Se + 1
                        RS = DECODE(HTable,0);
                        R = bitshift(RS,-4);             % run_length
                        S = bitand(RS,2^4-1);            % code size
                        if S == 0
                            if R == 15
                                K = K + 16;           % ZRL code
                            else
                                % EOB, End of Band
                                EOBRUN = bitshift(1,R) + get_bits(R)-1;
                                break                 % EOB��End of Block
                            end
                        else
                            K = K + R;
                            Temp = get_bits(S);
                            Decoder.Coes{c}(K,blk_cnt) = bitshift(EXTEND(Temp, S),Decoder.Al,'int64');
                            K = K + 1;
                        end
                    end
                end
            end
        end
        Decoder.EOBRUN = EOBRUN;
    end
    function DecodeDCRefine(id)
        %{
            DCϵ�����������룬���Ƚ�һ���ƽ���ʼ����
        %}
        Width = Decoder.channel_info.MCUs_per_row;
        P = bitshift(1,Decoder.Al);
        for i = 1:Decoder.channel_info.channels
            c = Decoder.channel_info.id(i);
            MCU_width = Decoder.sample_factor(1,c);
            MCU_height = Decoder.sample_factor(2,c);
            row =  MCU_height * fix(id/Width);
            col = MCU_width * mod(id,Width);
            for j = 1:MCU_height
                for k = 1:MCU_width
                    blk_idx = row*Width*MCU_width + col + k;
                    if get_bits(1)
                        temp = Decoder.Coes{c}(1,blk_idx);
                        if temp > 0
                            Decoder.Coes{c}(1,blk_idx) = ...
                                bitor(temp,P);
                        else
                            Decoder.Coes{c}(1,blk_idx) = ...
                                bitor(temp + 2^8,P)-2^8;
                        end
                    end
                end
                row = row + 1;
            end
        end
    end
    function DecodeACRefine(mcu_cnt)
        %{
            ACϵ�����������룬��λ����������һ����ͬ��ע����ź�λ��
        %}
        EOBRUN = Decoder.EOBRUN;
        c = Decoder.channel_info.id;
        tbl_id = Decoder.channel_info.ac_tbl_id+1;
        HTable = Decoder.ac_huff_tbl(tbl_id);
        blks = Decoder.channel_info.blks_in_MCU;
        c1 = Decoder.channel_info.MCU_width;
        c2 = Decoder.channel_info.MCU_height;
        P = bitshift(1,Decoder.Al);
        for i = 1:c1
            for j = 1:c2
                temp = (i-1)*c2 + j;
                blk_cnt = mcu_cnt*blks + temp;
                K = Decoder.Ss + 1;
                if EOBRUN == 0
                    while K <= Decoder.Se+1
                        Temp = Decoder.Coes{c}(K,blk_cnt);
                        RS = DECODE(HTable,0);
                        R = fix(RS / 16);                    % run_length
                        S = mod(RS , 16);                    % code size
                        if S
                            if(S~=1)
                                error('Bad Huffman Code!');
                            end
                            if get_bits(1)
                                S = P;
                            else
                                S = -P;
                            end
                        else
                            if R~=15
                                EOBRUN = 2 ^ R;
                                if R
                                    EOBRUN = EOBRUN + get_bits(R);
                                end
                                break
                            end
                        end
                        while K <= Decoder.Se+1
                            Temp = Decoder.Coes{c}(K,blk_cnt);
                            if Temp ~= 0
                                if get_bits(1)
                                    if bitand(abs(Temp),P)
                                        if Temp>0
                                            Temp = bitor(Temp, P);
                                        else
                                            Temp = bitor(temp + 2^8,P)-2^8;
                                        end
                                    end
                                end
                            else
                                R = R-1;
                                if R < 0
                                    break
                                end
                            end
                            Decoder.Coes{c}(K,blk_cnt)=Temp;
                            K = K + 1;
                        end
                        if S
                            Temp = S;
                        end
                        Decoder.Coes{c}(K,blk_cnt)=Temp;
                        K = K + 1;
                    end
                end
                if EOBRUN >0
                    while K <= Decoder.Se+1
                        Temp = Decoder.Coes{c}(K,blk_cnt);
                        if Temp ~=0
                            if get_bits(1)
                                if bitand(abs(Temp),P) ==0
                                    if Temp >0
                                        Temp = Temp + P;
                                    else
                                        Temp = Temp - P;
                                    end
                                end
                            end
                        end
                        Decoder.Coes{c}(K,blk_cnt)=Temp;
                        K = K+1;
                    end
                    EOBRUN = EOBRUN - 1;
                end
                Decoder.EOBRUN = EOBRUN;
            end
        end
    end
    function S = DECODE(Htable,isDC)
        %{
            ��ñ���������Ҫ��bits��
            Obtaining the value of symbol encoded by next bits from Huffman table

            Procedure which returns the 8-bit value associated with the next
            Huffman code in the compressed image data.

            Args
                identifier id with which Huffman table the value is encoded

                        00  Huffman table of Luminance DC DIFF
                        10  Huffman table of Luminance AC coefficient
                        01  Huffman table of Chrominance DC DIFF
                        11  Huffman table of Chrominance AC coefficient
            Returns
                CODE_LEN�� Number of bits with which the image compressed
                           data is encoded
        %}
        if ~exist("isDC",'var')
            isDC = true;
        end
        if Decoder.bits_in_buffer < HUFF_LOOKAHEAD
            jpeg_fill_buffer();
        end
        buffer = Decoder.temp_buffer;
        bits_left = Decoder.bits_in_buffer;  % ������ʣ��δ�õ���bit
        look = bitand(bitshift(buffer,-(bits_left-HUFF_LOOKAHEAD)),255);
        % ����һ���ֽ���Ϊ���ұ�����������ҽ���ĸ�λΪ������Ҫ���ֽ�����
        % �Ͱ�λΪ�������ֵ(DCϵ��Ϊ����DIFF��Ҫ���볤��ACϵ��Ϊ�г���RS)
        temp = Htable.lookup(look+1);
        nb = bitshift(temp,-HUFF_LOOKAHEAD);
        if nb <= HUFF_LOOKAHEAD
            S = bitand(temp,bitshift(1,HUFF_LOOKAHEAD) - 1);
            Decoder.bits_in_buffer = bits_left-nb;
            Decoder.temp_buffer = bitand(buffer,2^Decoder.bits_in_buffer-1);
        else
            S = jpeg_huff_decode(Htable,nb);
        end
        if isDC
            DIFF = get_bits(S);
            S = bitshift(EXTEND(DIFF, S),Decoder.Al,'int16');
        end
    end

    function ParseSOF()
        %{
            Start of Frame
        %}
        data_length = ReadTwoBytes();
        data_length = data_length - 2;
        ImgInfo.Precision = ReadOneByte();
        ImgInfo.Height = double(ReadTwoBytes());
        ImgInfo.Width = double(ReadTwoBytes());
        ImgInfo.Channels = ReadOneByte();
        data_length = data_length - 6;
        fprintf('Image precision is %d bit,with shape of %d x %d x %d.\n',...
            ImgInfo.Precision,ImgInfo.Height,ImgInfo.Width,ImgInfo.Channels);
        Idx = 1;
        for i = 1:ImgInfo.Channels
            % ͨ����źͲ���ϵ��
            c_id = ReadOneByte();
            sample_factor = ReadOneByte();
            sf_h = bitshift(sample_factor,-4);
            sf_v = bitand(sample_factor,2^4-1);
            Decoder.sample_factor(:,c_id)=[sf_h,sf_v];

            % ��������
            Decoder.quanti_tbl_idx(c_id) = ReadOneByte() + 1;

            % ��ͨ����Ϣ
            c1 = sf_v*sf_h;
            Decoder.blks_in_MCU(i) = c1;
            for j = 1:c1
                Decoder.member_ship(Idx) = i;
                Idx = Idx + 1;
            end
            data_length = data_length-3;
        end

        sf_h_max = max(Decoder.sample_factor(1,:));
        sf_v_max = max(Decoder.sample_factor(2,:));
        for i = 1:ImgInfo.Channels
            sf_h = sf_h_max/Decoder.sample_factor(1,i);
            sf_v = sf_v_max/Decoder.sample_factor(2,i);
            c1 = ceil(ImgInfo.Width/8/sf_h);
            c2 = ceil(ImgInfo.Height/8/sf_v);
            Decoder.blks_per_row(i) = c1;
            Decoder.blks_per_col(i) = c2;
            Decoder.Coes{i} = zeros(JPEG_BLOCK_SIZE,c1*c2);
        end
        Decoder.MCUs_per_row = min(Decoder.blks_per_row);
        Decoder.MCUs_per_col = min(Decoder.blks_per_col);
        assert (data_length == 0,'Bad Data Length of SOF.')
        Decoder.next_marker = 0;
    end
    function ParseAPPn()
        SegLength = ReadTwoBytes();
        DataSegment = ReadNBytes(jpegfile,SegLength - 2);
        file_format = char();

        % �����ļ��洢��ʽ��һ�����ͨ�ø�ʽJFIF
        if Marker == hex2dec('E0')
            for Idx =1:SegLength - 2
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
            �������������
        %}

        SegLength=ReadTwoBytes();
        SegLength = SegLength-2;
        while SegLength > 0
            temp = ReadOneByte();
            dqt_precision = bitshift(temp,-4);
            identifier = bitand(temp,2^4-1);

            SegLength = SegLength -1;
            % ע����С�˴洢����ȡ�ʹ洢˳���෴
            if dqt_precision == 0
                Qtbl = ReadNBytes(jpegfile,64);
            else
                Qtbl = ReadNBytes(jpegfile,64,'uint16');
            end
            SegLength = SegLength - 64;
            Qtbl = inverse_zigzag(Qtbl);
            Decoder.dqt_ids{identifier+1} = Qtbl;
        end
        Decoder.next_marker = 0;
    end
    function ParseSOS()
        %{
            ����ɨ��ͷ����������˴�����ͨ����Ϣ��
            һ���ж��ٸ�ͨ����ÿ��ͨ��ʹ���ĸ�Huffman��
            DCϵ�����ĸ���ACϵ�����ĸ����۽������Ƶ
            ��ָ������αƽ��㷨�ľ��ȡ�
        %}

        % ��ʼ��
        ResetDecoder();
        data_length = ReadTwoBytes();
        channels = ReadOneByte();
        Decoder.channel_info.channels = channels;
        data_length = data_length - 3;
        if channels == 1
            %{
                single component,non-interleaved save mode,
                Noninterleaved (single-component) scan 
            %}
            % ��ǰSOS�������ͨ���ź���Ҫʹ�õ�huffman��ID
            c = ReadOneByte();
            id = ReadOneByte();
            Decoder.channel_info.id = c;
            Decoder.channel_info.dc_tbl_id = bitshift(id,-4);
            Decoder.channel_info.ac_tbl_id = bitand(id,2^4-1);

            % ��ǰSOS�д�����ͨ����MCU������ÿ��MCU�е�blocks����
            c1 = Decoder.blks_per_row(c);
            c2 = Decoder.blks_per_col(c);
            Decoder.channel_info.MCUs_per_row = c1;
            Decoder.channel_info.MCUs_per_col = c2;
            Decoder.channel_info.MCUs = c1*c2;
            Decoder.channel_info.MCU_width = 1;
            Decoder.channel_info.MCU_height = 1;
            Decoder.channel_info.blks_in_MCU = 1;
            data_length = data_length-2;
        else
            for i = 1:channels
                c = ReadOneByte();
                id = ReadOneByte();
                Decoder.channel_info.id(i) = c;
                Decoder.channel_info.dc_tbl_id(c)= bitshift(id,-4);
                Decoder.channel_info.ac_tbl_id(c) = bitand(id,2^4-1);
                data_length = data_length-2;
            end
            Decoder.channel_info.MCUs_per_row = Decoder.MCUs_per_row;
            Decoder.channel_info.MCUs_per_col = Decoder.MCUs_per_col;
            Decoder.channel_info.blks_in_MCU = sum(Decoder.blks_in_MCU);
        end

        % ָ���۽�����ĳ�ʼƵ������ֹƵ��
        Decoder.Ss = ReadOneByte();
        Decoder.Se = ReadOneByte();

        % ָ����αƽ������λ�����λ
        Decoder.Al = ReadFourBits();
        Decoder.Ah = ReadFourBits();
        data_length = data_length-3;

        if Decoder.Ss == 0
            if Decoder.Ah == 0
                Decoder.DecodeMCU = @DecodeDCFirst;
            else
                Decoder.DecodeMCU = @DecodeDCRefine;
            end
            % ACϵ������
        else
            if Decoder.Ah == 0
                Decoder.DecodeMCU = @DecodeACFirst;
            else
                Decoder.DecodeMCU = @DecodeACRefine;
            end
        end

        % debug�۲��ã�release�汾ɾ��
        if Decoder.Se == 0
            Band = '��DCϵ��';
        else
            Band = '��ACϵ��';
        end
        if channels == 3
            STR = 'ȫ��ͨ��';
        else
            switch Decoder.channel_info.id
                case 1
                    STR ='����ͨ��';
                case 2
                    STR = 'ɫ��ͨ��1';
                case 3
                    STR = 'ɫ��ͨ��2';
            end
        end
        assert(data_length==0,"Bad Data Length of SOS.")
        Decoder.next_marker = 0;
        fprintf('��%d��SOS.����%s\n',Decoder.SOScnt,strcat(STR,Band))
    end
    function [VALUE] = RECEIVE(SSSS)
        %{
            �����볤��÷���ֵ�������ֵ��û�и����ģ�����Ǹ��������ص��Ǹ�����
            ���롣

            Read code_length bits additional value

            Procedure which places the next SSSS bits of the entropy-coded
            segment into the LSB of DIFF as additional bits.It calls the
            NEXTBIT and returns the value of DIFF to the calling procedure.

            Args
                SSSS Code size with which the value of DIFF is encoded

            Returns
                VALUE: DIFF encoded in two's complement
        %}
        VALUE=0;
        for i = 1:SSSS
            if Decoder.next_marker == 0
                VALUE = bitshift(VALUE,1) + NEXTBIT();
            else
                break;
            end
        end
    end
    function [V] = EXTEND(V, T)
        %{
            ��չ���ţ��������ʾ��ֵ����Ϊ��ԭ���ĸ���ֵ

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
    function Bit = NEXTBIT()
        if ~Decoder.next_marker
            C = Decoder.bitsleft;
            if C == 0
                Decoder.unusebits = ReadOneByte();
                if Decoder.unusebits == hex2dec('FF')
                    nextbyte = ReadOneByte();
                    while nextbyte == hex2dec('FF')
                        nextbyte = ReadOneByte();
                    end
                    if nextbyte == 0
                        C = 8;
                    elseif ismember(nextbyte,[196,217,218])
                        Decoder.next_marker = nextbyte;
                        Decoder.unusebits = -1;
                        Bit = 0;
                        return;
                    else
                        error('Illegal marker!')
                    end
                else
                    C = 8;
                end
            end
            C = C - 1;
            Bit = bitshift(Decoder.unusebits,-C);
            Decoder.unusebits = mod(Decoder.unusebits,2^C);
            Decoder.bitsleft = C;
        end
    end
    function ParseDHT()
        bits = zeros(1,17);
        huffval = zeros(1,256);
        data_length = ReadTwoBytes();
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
        for l = 2:17
            cnt = htbl.bits(l);
            if cnt < 0 || index + cnt > 256   % protect against table overrun %}
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

        % Figure C.2: generate the codes themselves %}
        % We also validate that the counts represent a legal Huffman code tree. %}

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
         it must still fit in si bits, since no code is allowed to be all ones.
        %}
        % Figure F.15: generate decoding tables for bit-sequential decoding

        index = 1;
        for l = 2:17
            num =htbl.bits(l); % num of code with length == l-1(Matlab starts with 1)
            if num
                %{
                    valoffset(l) = huffval() index of 1st symbol of code length l,
                    minus the minimum code of length l
                    valoffset(length)��huffval���볤Ϊlength�ĵ�һ�����ŵ�������ȥ
                    �볤Ϊlength����С������
                %}
                min_code = huffcode(index); % minimum code of length l-1
                dtbl.valoffset(l) = index - min_code;
                index = index + num; % index of 1st symbol of code with next length l
                dtbl.maxcode(l) = huffcode(index - 1); % maximum code of length l
            else
                dtbl.maxcode(l) = -1;    % -1 if no codes of this length
            end
        end
        % ������Ϊ�˱���������ʱ���������
        dtbl.valoffset(18) = 0; % value in valoffset is 1 more than C/C++/Python
        dtbl.maxcode(18)= hex2dec('FFFFF'); %{ ensures jpeg_huff_decode terminates %}
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
        % �������п��ܵ�huffman���ֵ��볤�����huffman code size������8��bit�������ֱ��
        % ������������8���ֽھ���Ҫ�����Ӧ��bit���õ���Ӧ��ֵ��
        index = 1;
        for L = 1:HUFF_LOOKAHEAD
            for cnt = 1:htbl.bits(L+1)
                % l = current code's length, p = its index in huffcode() & huffval()
                % Generate left-justified code followed by all possible bit sequences
                lookbits = bitshift(huffcode(index) ,(HUFF_LOOKAHEAD - L));%code with codesize == L
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
        % HUFF_DECODE has determined that the code is at least min_bits */
        % bits long, so fetch that many bits in one swoop. */
        code = get_bits(l);
        % Collect the rest of the Huffman code one bit at a time. */
        while code > htbl.maxcode(l+1)
            code = bitor(bitshift(code,1),get_bits(1));
            l = l + 1;
        end
        % With garbage input we may reach the sentinel value l = 17. */
        if (l > 16)
            warning("Corrupt JPEG data: bad Huffman code");
        end
        try
            re =  htbl.pub.huffval(code + htbl.valoffset(l+1));
        catch
            return
        end
    end
    function jpeg_fill_buffer()
        [buffer,nbits] = deal(Decoder.temp_buffer,Decoder.bits_in_buffer);
        while nbits <= 24   % ÿ�ζ����bit����������32��
            buffer = bitor(bitshift(buffer,HUFF_LOOKAHEAD),RECEIVE(8));
            nbits = nbits + 8;
        end
        [Decoder.temp_buffer,Decoder.bits_in_buffer] = deal(buffer,nbits);
    end
    function re = get_bits(n)
        % ����n��bit����
        if Decoder.bits_in_buffer < n
            jpeg_fill_buffer()
        end
        [buffer,nbits] = deal(Decoder.temp_buffer,Decoder.bits_in_buffer);
        nbits = nbits-n;
        re = bitshift(buffer,-nbits);
        buffer = bitand(buffer,2^nbits-1);
        [Decoder.temp_buffer,Decoder.bits_in_buffer] = deal(buffer,nbits);
    end
end


