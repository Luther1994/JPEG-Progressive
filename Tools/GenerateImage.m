function GenerateImage(Decoder,ImgInfo)
%{
   Generate and Display Image from Decoded Coefficients.
%}
switch ImgInfo.Precision
    case 8
        type = 'uint8';
    case 16
        type = 'uint16';
end
result = zeros(ImgInfo.Height,ImgInfo.Width,ImgInfo.Channels,type);
sf_v_max = max(Decoder.sample_factor(1,:));
sf_h_max = max(Decoder.sample_factor(2,:));
for row = 0:Decoder.MCUs_in_ver-1
    for col = 0:Decoder.MCUs_in_hor-1
        mcu_cnt = row*Decoder.MCUs_in_hor+col;
        for c = 1:ImgInfo.Channels
            DQTable = Decoder.dqt_ids{Decoder.quanti_tbl_idx(c)};
            blks = Decoder.blks_in_MCU_per_c(c);
            sf_v = Decoder.sample_factor(1,c);
            sf_h = Decoder.sample_factor(2,c);
            COE = Decoder.Coes{c}(:,mcu_cnt*blks+1:(mcu_cnt+1)*blks);
            re = [];
            for i = 0:sf_v-1
                for j = 0:sf_h-1
                    blk_cnt = i*Decoder.sample_factor(1,c)+j;
                    block = COE(:,blk_cnt+1);
                    temp = inverse_zigzag(block);
                    temp = idct2(temp .* DQTable)+ 2 ^ (ImgInfo.Precision - 1);
                    re(i*8+1:(i+1)*8,j*8+1:(j+1)*8) = temp;
                end
            end
            mcu = UPSAMPLE(re,[sf_h_max/sf_h,sf_v_max/sf_v]);
            result(row*8*sf_v_max+1:(row+1)*8*sf_v_max,...
                col*8*sf_h_max+1:(col+1)*8*sf_h_max,c)=mcu;
        end
    end
end
if ndims(result)==3 
    result = yuv2rgb(result);
end
imshow(result);
imwrite(result,Decoder.savefile);
end