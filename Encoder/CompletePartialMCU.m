function  [img]=CompletePartialMCU(img,sample_factors)
% 填充函数，对长宽不能整除8以及不能整除采样因子的图像进行边缘填充

[H,W] = size(img);
blocks_H = ceil(H/8);
blocks_W = ceil(W/8);

if mod(H,8) ~= 0
    for row = 1:mod(H,8)
        img(end+1,:) = img(end,:);  
    end
end
if  mod(W,8) ~= 0
    for col = 1:mod(W,8)
        img(:,end+1) = img(:,end); %#ok<*AGROW> 
    end
end

if mod(blocks_H,sample_factors(1)) ~=0
    img(end+1:end+8,:) = img(end-7:end,:);
end
if mod(blocks_W,sample_factors(2)) ~=0
    img(:,end+1:end+8) = img(:,end-7:end);
end

end