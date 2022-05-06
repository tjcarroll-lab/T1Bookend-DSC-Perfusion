function images = ReadDICOM(path_source,n_file_start,n_file_end,increment_file)

% images = ReadDICOM(path_source,n_file_start,n_file_end,increment_file)
% authors: Wanyong
%
% description: Read DCM images using dicomread.
%
% status: stable
%
% versions
%   06-09-21 (WYS): initial version



% parse arguments
if nargin < 4
    increment_file = 1;
    
    if nargin < 3
        n_file_end = 0;
        
        if nargin < 2
            n_file_start = 1;
        end
    end
end

% read files
files = dir(path_source);

if n_file_end == 0
    n_file_end = length(files)-2;
end

path_source = AddTrailingSlash(path_source);
count = 1;
for i = n_file_start:increment_file:n_file_end
    images{count} = dicomread([path_source num2str(i) '.dcm']);
    count = count + 1;
end
