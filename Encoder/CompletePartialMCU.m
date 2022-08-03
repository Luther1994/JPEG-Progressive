function  [img]=CompletePartialMCU(img,sample_factors)
% ��亯�����Գ���������8�Լ����������������ӵ�ͼ����б�Ե���

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