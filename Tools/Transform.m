function [result] = Transform(src,inverse)
pos = [1,1];
if inverse
    assert (ndims(src) == 1 || size(src,2) == 1 || size(src,1) == 1 ,...
        'Input zigzaged coefficients vector must be one dimensions!')
    SIZE =sqrt(length(src));
    result = zeros(SIZE, SIZE);
    [cols,rows] = size(result);
    result(1,1) = src(1);
else
    assert (ismatrix(src), 'Input coefficients matrix must be two dimensions!')
    [cols,rows] = size(src);
    result = zeros(rows*cols,1);
    result(1) = src(pos(1),pos(2));
end

index = 2;
transform();

function step(direction)
if direction == 0
    % move down
    pos(1) = pos(1)+ 1;
elseif direction == 1
    % move right
    pos(2) = pos(2)+ 1;
elseif direction == 3
    % move in diagonal down
    pos(1) = pos(1) - 1;
    pos(2) = pos(2)+ 1;
elseif direction == 2
    % move in diagonal up
    pos(1) = pos(1) + 1;
    pos(2) = pos(2) - 1;
end
if inverse
    result(pos(1),pos(2)) = src(index);
else
    result(index) = src(pos(1),pos(2));
end

index = index +1;
end

function transform()
% 0 down 1 right 2 diagonal_down 3 diagonal_up
while index < len(src)
    if pos(2) <= cols-1
        step(1)
    else
        step(0)
    end
    while pos(1) < rows  &&  1 < pos(2)
        step(2)
    end
    if pos(1) <= rows-1
        step(0)
    else
        step(1)
    end
    while 1 < pos(1) && pos(2) < cols
        step(3)
    end
end
end
end

function [LEN] = len(src)

    [rows,cols] = size(src);
    LEN = rows*cols;
end




