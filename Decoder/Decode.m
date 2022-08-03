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
            Decoder.next_marker = 0;
        else
            switch Marker
                case Markers.SOS
                    Decoder.SOScnt = Decoder.SOScnt+1;
                    DecodeImage();
                case Markers.EOI
                    fprintf('EOI found,decode terminate.image saved as %s\n',savefile);
                    %                     Decoder.Coes(1,:) = 0;
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
                    get_dht();
                    Decoder.next_marker = 0;
                case Markers.DQT
                    get_dqt();
                    Decoder.next_marker = 0;
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
        Decoder.block_idx = [1,1,1];
        for i = 1:Decoder.MCUs
           DecodeMCU();
        end
        %���һ�ν����Ҫ�������л�δ�����bits������
        Decoder.bitsleft = 0;
    end
    function DecodeMCU()
        %{
            ���Ľ��뺯��������ָ����MCU
        %}

        % ����SOS��ָ����Ƶ���;��ȷֱ���н���

        % DCϵ������
        if Decoder.Ss == 0
            if Decoder.Ah == 0
                DecodeDCFirst();
            else
                DecodeDCRefine();
            end
            % ACϵ������
        else
            if Decoder.Ah == 0
                DecodeACFirst();
            else
                DecodeACRefine();
            end
        end
    end
    function DecodeDCFirst()
        %{
            DCϵ����һ�ν��룬���λ
        %}

        % DCϵ������ʱ��ȫ��ͨ��һ�ν���
        for k = 1:Decoder.blks_in_MCU
            c_idx = Decoder.member_ship(k);
            Htable = Decoder.dc_huff_tbl(Decoder.DCTbls(c_idx)+1);
            % Huffman������̣������JPEG�ĵ�
            DIFF = DECODE(Htable);
            % DCϵ��ͨ����ִ洢
            Decoder.LastDCVal(c_idx) = Decoder.LastDCVal(c_idx) + DIFF;
            % ������д�뻺����
            Decoder.Coes{c_idx}(1,Decoder.block_idx(c_idx))= ...
            Decoder.LastDCVal(c_idx);
            Decoder.block_idx(c_idx) = Decoder.block_idx(c_idx)+1;
        end
    end
    function DecodeDCRefine()
        %{
            DCϵ�����������룬���Ƚ�һ���ƽ���ʼ����
        %}
        P = 2^Decoder.Al;
        for k = 1:sum(Decoder.blks_in_MCU_per_c)
            if RECEIVE(1)
                temp = Decoder.Coes(1,Decoder.block_idx+1);
                % ����Ҫ�����ѽ����ϵ����ȷ������
                % λ������൱�ڼӷ�
                if temp > 0
                    Decoder.Coes(1,Decoder.block_idx+1) = ...
                        bitor(temp,P);
                else
                    Decoder.Coes(1,Decoder.block_idx+1)  = ...
                        bitor(temp + 2^8,P)-2^8;
                end
            end
            Decoder.block_idx = Decoder.block_idx + 1;
        end
    end
    function DecodeACFirst()
        %{
            ACϵ����һ�ν��룬���λ
        %}
        HMax = max(Decoder.sample_factor(1,:));
        blockwidth = Decoder.sample_factor(1, Decoder.Channels);
        blockheight = Decoder.sample_factor(2, Decoder.Channels);

        % ˮƽ�����DU����
        HB = max(Decoder.blks_in_hor*blockwidth/HMax);
        P = 2^Decoder.Al;

        % EOB run length,�����ǰƵ���Ժ��Ƶ����Ϊ0������EOBRUN��ʾ
        EOBRUN = Decoder.EOBRUN;

        % ��/�е�����
        [RowIdx,ColIdx] = deal(Decoder.RowIdx,Decoder.ColIdx);

        % ˮƽ/��ֱ�����ϵĵڼ���MCU
        HMcuIdx = fix(ColIdx/blockwidth);
        VMcuIdx = fix(RowIdx/blockheight);
        
        % MCU ����
        McuIdx = VMcuIdx * max(Decoder.blks_in_hor)/HMax + HMcuIdx;

        % ���ݲ�ͬ��ͨ����ÿ�ν����ϵ��λ�ò�ͬ������ 1234 5 6 ��˳��
        % ����Coes�����С�ͬʱָ��Huffman��
        switch  Decoder.Channels
            case 1
                tbl_id = 1;
                if mod(RowIdx,blockheight)
                    temp = 3;
                else
                    temp = 1;
                end
            case 2
                tbl_id = 2;
                temp = 5;

            case 3
                tbl_id = 2;
                temp = 6;
        end
        HTable = Decoder.ac_huff_tbl(Decoder.ACTbls(tbl_id)+1);
        % ��ʼ���룬��������ĵ����۽����밴�շǽ���ģʽ���д洢����ÿ��
        % ֻ����һ��ͨ����
        StartBlockIdx = McuIdx * sum(Decoder.blks_in_MCU_per_c) + temp;
        EndBlockIdx = StartBlockIdx + blockwidth -1;
        for i = StartBlockIdx : EndBlockIdx
            K = Decoder.Ss+1;
            if EOBRUN > 0
                EOBRUN = EOBRUN - 1;
            else
                while K <= Decoder.Se + 1
                    RS = DECODE(HTable,0);
                    Decoder.temp(end+1) = RS;
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
                        Decoder.Coes(K,i) = EXTEND(Temp, S) * P;
                        K = K + 1;
                    end
                end
            end
            Decoder.block_idx = Decoder.block_idx + 1;
            Decoder.EOBRUN = EOBRUN;
            RowIdx = fix(Decoder.block_idx/HB);
            ColIdx = mod(Decoder.block_idx,HB);
        end
        Decoder.RowIdx = RowIdx;
        Decoder.ColIdx = ColIdx;
    end
    function DecodeACRefine()
        %{
            ACϵ�����������룬��λ����������һ����ͬ��ע����ź�λ��
        %}
        HMax = max(Decoder.sample_factor(1,:));
        blockwidth = Decoder.sample_factor(1,Decoder.Channels);
        blockheight = Decoder.sample_factor(2,Decoder.Channels);
        HB = max(Decoder.blks_in_hor*blockwidth/HMax);
        P = 2^Decoder.Al;
        EOBRUN = Decoder.EOBRUN;
        [RowIdx,ColIdx] = deal(Decoder.RowIdx,Decoder.ColIdx);
        VMcuIdx = fix(RowIdx/blockheight); % ��ֱ�����ϵĵڼ���MCU
        HMcuIdx = fix(ColIdx/blockwidth);
        McuIdx = VMcuIdx * max(Decoder.blks_in_hor)/HMax + HMcuIdx;
        switch Decoder.Channels
            case 1
                HTable = Decoder.AC_Lu_Table;
                if mod(RowIdx,blockheight)
                    S = 3;
                else
                    S = 1;
                end
            case 2
                S = 5;
                HTable = Decoder.AC_Chr_Table;
            case 3
                S = 6;
                HTable = Decoder.AC_Chr_Table;
        end
        StartBlockIdx = McuIdx * sum(Decoder.blks_in_MCU_per_c) + S;
        for i = StartBlockIdx : StartBlockIdx + blockwidth - 1
            K = Decoder.Ss + 1;
            if EOBRUN == 0
                while K <= Decoder.Se+1
                    Temp = Decoder.Coes(K,i);
                    RS = DECODE(HTable);
                    R = fix(RS / 16);                    % run_length
                    S = mod(RS , 16);                    % code size
                    if S
                        if(S~=1)
                            error('Bad Huffman Code!');
                        end
                        if RECEIVE(1)
                            S = P;
                        else
                            S = -P;
                        end
                    else
                        if R~=15
                            EOBRUN = 2 ^ R;
                            if R
                                EOBRUN = EOBRUN + RECEIVE(R);
                            end
                            break
                        end
                    end
                    while K <= Decoder.Se+1
                        Temp = Decoder.Coes(K,i);
                        if Temp ~= 0
                            if RECEIVE(1)
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
                        Decoder.Coes(K,i)=Temp;
                        K = K + 1;
                    end
                    if S
                        Temp = S;
                    end
                    Decoder.Coes(K,i)=Temp;
                    K = K + 1;
                end
            end
            if EOBRUN >0
                while K <= Decoder.Se+1
                    Temp = Decoder.Coes(K,i);
                    if Temp ~=0
                        if RECEIVE(1)
                            if bitand(abs(Temp),P) ==0
                                if Temp >0
                                    Temp = Temp + P;
                                else
                                    Temp = Temp - P;
                                end
                            end
                        end
                    end
                    Decoder.Coes(K,i)=Temp;
                    K = K+1;
                end
                EOBRUN = EOBRUN -1;
            end
            Decoder.block_idx = Decoder.block_idx + 1;
            Decoder.EOBRUN = EOBRUN;
            RowIdx = fix(Decoder.block_idx/HB);
            ColIdx = mod(Decoder.block_idx,HB);
        end
        Decoder.RowIdx = RowIdx;
        Decoder.ColIdx = ColIdx;
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
            ���������ͼ��ľ��ȡ��ߴ硢ͨ�����ã�����ÿ��ͨ���Ĳ���ϵ����ָ��������
        %}
        data_length = ReadTwoBytes();
        data_length = data_length-2;
        % ͼ�񾫶ȡ��ߴ硢ͨ����
        ImgInfo.Precision =ReadOneByte();
        data_length = data_length - 1; 
        ImgInfo.Height = double(ReadTwoBytes());
        ImgInfo.Width = double(ReadTwoBytes());
        ImgInfo.Channels = ReadOneByte();
        data_length = data_length - 5;
        fprintf('Image precision is %d bit,with shape of %d x %d x %d.\n',...
            ImgInfo.Precision,ImgInfo.Height,ImgInfo.Width,ImgInfo.Channels);

        % ͨ�����ã���������ϵ���͸�ͨ��ʹ�õ�������
        for n = 1:ImgInfo.Channels
            ChannelIdx = ReadOneByte();
            Decoder.sample_factor(:,ChannelIdx)=...
                [ReadFourBits(),ReadFourBits()];
            Decoder.quanti_tbl_idx(ChannelIdx) = ReadOneByte() + 1;
            data_length = data_length-3;
        end
        
        % ˮƽ����ʹ�ֱ����Ĳ���ϵ��
        fh_max = max(Decoder.sample_factor(1,:));
        fv_max = max(Decoder.sample_factor(2,:));

        % ˮƽ�������ֱ�����MCU����
        Decoder.MCUs_in_ver = ceil(ImgInfo.Height/fh_max/8);
        Decoder.MCUs_in_hor = ceil(ImgInfo.Width/fv_max/8);
        Decoder.MCUs = Decoder.MCUs_in_ver.*Decoder.MCUs_in_hor;

        % ÿ��MCU�е�DU��blks������������ͼ���е�MCU����
        Decoder.blks_in_MCU_per_c = prod(Decoder.sample_factor,1);
        Decoder.blks_in_MCU = sum(Decoder.blks_in_MCU_per_c);
        Idx = 1;
        for i = 1:ImgInfo.Channels
            for j = 1:Decoder.blks_in_MCU_per_c(i)
                Decoder.member_ship(Idx) = i;
                Idx = Idx + 1;
            end
            % �۽�����ͻ��߱����һ����Ҫ�������ϵ���Ĵ洢����Ҫһ���ϴ�Ļ���
            % ���洢���е�ϵ��
            Decoder.Coes{i} = zeros(JPEG_BLOCK_SIZE,...
                Decoder.MCUs*Decoder.blks_in_MCU_per_c(i));
        end
        assert (data_length == 0,'Bad Data Length of SOF.')
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
        fprintf('File is of format %s.\n',file_format);
    end
    function get_dqt()
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
        [Decoder.ACTbls,Decoder.DCTbls] = deal([]);
        Channels = ReadOneByte();
        data_length = data_length - 3;
        % ָ���ر���ı����
        for i = 1:Channels
            ChannelIdx = ReadOneByte();
            identifier = ReadOneByte();
            Decoder.DCTbls(ChannelIdx) = bitshift(identifier,-4);
            Decoder.ACTbls(ChannelIdx) = bitand(identifier,15);
            Decoder.Channels = [Decoder.Channels ChannelIdx];
            data_length = data_length-2;
        end
        % ָ���۽�����ĳ�ʼƵ������ֹƵ��
        Decoder.Ss = ReadOneByte();
        Decoder.Se = ReadOneByte();
        
        if Decoder.Se == 0
            Band = '��DCϵ��';
        else
            Band = '��ACϵ��';
        end
        switch sum(Decoder.Channels)
            case 6
                STR = 'ȫ��ͨ��';
            case 1
                STR ='����ͨ��';
            case 2
                STR = 'ɫ��ͨ��1';
            case 3
                STR = 'ɫ��ͨ��2';
        end
        % ָ����αƽ������λ�����λ
        Decoder.Al = ReadFourBits();
        Decoder.Ah = ReadFourBits();
        data_length = data_length - 1;
        assert(data_length == 0,'Bad Data Length of SOS.')
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
    function get_dht()
        bits = zeros(1,17);
        huffval = zeros(1,256);
        length = ReadTwoBytes();
        length = length - 2;
        while length > 16
            index = ReadOneByte();
            isAC = bitshift(index,-4);
            index = bitand(index,2^4-1);
            bits(1) = 0;
            count = 0;
            for i = 2:17
                bits(i) = ReadOneByte();
                count = count +  bits(i);
            end
            length = length - 17;
            if count > 256 || count > length
                error('Bogus Huffman table definition')
            end
            for i = 1:count
                huffval(i) = ReadOneByte();
            end
            length =length - count;
            table.bits = bits;
            table.huffval = huffval;
            if isAC         % AC table definition
                Decoder.ac_huff_tbl(index+1)=jpeg_make_derived_tbl(0,table);
            else            % DC table definition
                Decoder.dc_huff_tbl(index+1)=jpeg_make_derived_tbl(1,table);
            end
        end
        if (length ~= 0)
            error('Counts of bits left must be 0.')
        end
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


