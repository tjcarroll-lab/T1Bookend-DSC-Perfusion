function [images, names] = IA_ImportImages(varargin)

if length(varargin) < 3
    printf('error')
else
    images = varargin{1};
    names  = varargin{2};
    for i = 3:2:length(varargin)
        Index_adding = 1;
        for j = 1:length(names)
            if strcmp(varargin{i},names{j})
                images{j} = varargin{i+1};
                Index_adding = 0;
            end
        end
        if Index_adding
            names{length(names)+1} = varargin{i};
            images{length(names)} = varargin{i+1};
        end
    end
end

clear Index_adding