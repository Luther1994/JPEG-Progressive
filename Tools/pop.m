function [item] = pop(Array,idx)
if (~exist('idx','var'))
    idx = 1;
end
item = Array(idx);
Array(idx)=[];
end
