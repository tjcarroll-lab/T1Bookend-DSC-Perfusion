function images = IA_ReadImages(path_source,n_file_start,n_file_end,increment_file)

% images = IA_ReadImages(path_source,n_file_start,n_file_end,increment_file)
% authors: Ty Cashen, Wanyong
%
% description: Read images for ImageAnalyzer.
% path_source is a string of path where images will be read, n_file_start
% and n_file_end are starting and ending number of images with a number of
% increment_file. output, images are a cell array.
%
% status: stable
%
% versions
%   ??-??-?? (TAC): initial version
%   06-09-21 (WYS): supported in Matlab 7.x.x 

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

% first determine if .mat file exists
files = what(path_source);

if ~isempty(files.mat)
    %handle_progress =  ProgressText('Reading MATLAB images...');
    
    images = load([path_source '/' files.mat{1}]);
    names_fields = fieldnames(images);
    images = getfield(images,names_fields{1});
    
    if ~iscell(images)
        images = Matrix2Cells(images);
    end
    if n_file_end == 0
        n_file_end = length(images);
    end
    images = images(n_file_start:increment_file:n_file_end);
    
else % otherwise use "FastReadDICOM", if it is availabe
    %handle_progress =  ProgressText('Reading dcm images...');
    verinfo = version;
    if str2num(verinfo) < 7 % if matlab 6.5.xx
        images = FastReadDICOM(path_source,n_file_start,n_file_end,increment_file);
    else 
        images = ReadDICOM(path_source,n_file_start,n_file_end,increment_file);
    end
end

%close(handle_progress)