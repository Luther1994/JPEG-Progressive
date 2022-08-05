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
for c = 1:ImgInfo.Channels
    DQTable = Decoder.dqt_ids{Decoder.quanti_tbl_idx(c)};
    sf_v = Decoder.sample_factor(1,c);
    sf_h = Decoder.sample_factor(2,c);
    COE = Decoder.Coes{c};
    height = Decoder.blks_per_col(c)*8;
    width = Decoder.blks_per_row(c)*8;
    re = zeros(height,width);
    id = 1;
    for row = 0:Decoder.blks_per_col(c)-1
        for col = 0:Decoder.blks_per_row(c)-1
            block = COE(:,id);
            id = id+1;
            temp = inverse_zigzag(block);
            temp = idct2(temp .* DQTable)+ 2 ^ (ImgInfo.Precision - 1);
            re(row*8+1:(row+1)*8,col*8+1:(col+1)*8) = temp;
        end
    end
    temp = UPSAMPLE(re,[sf_h_max/sf_h,sf_v_max/sf_v]);
    result(:,:,c) = temp(1:ImgInfo.Height,1:ImgInfo.Width);
end
if ndims(result)==3
    result = yuv2rgb(result);
end
imshow(result);
imwrite(result,Decoder.savefile);
end