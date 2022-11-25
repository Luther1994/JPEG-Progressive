img_path = 'DVA.jpg';
save_path = 'decode_re.jpg';
Decode(img_path,save_path);
%%
img = imread('DVA.jpg');
quality = 50;
Encode(img,'encode_re.jpg',quality);