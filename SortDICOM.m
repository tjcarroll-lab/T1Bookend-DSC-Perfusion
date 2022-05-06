function SortDICOM(path_source,path_target,convention)

% SortDICOM(path_source,path_target,convention)
%
% Sort DICOM Images
% author: Ty Cashen
%
% description: Sorts DICOM images from scanner into directories using a
% standardized directory and file naming convention.  "path_source" is the
% source path for the images.  "path_target" is the target path.
% "convention" is an optional argument that specifies the file naming
% convention.
%
% status: stable
%
% versions
%   03-06-09 (TAC): initial version
%   03-06-10 (TAC): fixed bug when series descriptions contain characters not allowed in directory names
%   03-06-18 (TAC): replaced SequenceName with SeriesDescription
%   03-09-27 (TAC): only reads DICOM (.dcm or .ima) files
%   03-11-05 (TAC): replaced "dicominfo" with "FastReadDICOMHeader" to optimize speed
%   04-04-04 (TAC): replaced "copyfile" with "fastcopyfile" to optimize speed

% set file naming convention to "standard" if not specified
if nargin == 2
    convention = 'standard';    
end

% add trailing "/" to paths if not already present
path_source = AddTrailingSlash(path_source);
path_target = AddTrailingSlash(path_target);

files = dir(path_source);

% sort images
handle_progress = ProgressBar('Sorting DICOM images...');

count = 1;
for i = 1:length(files);
    if ~files(i).isdir && (~strcmp(files(i).name(end-3),'.') || strcmp(files(i).name(end-2:end),'dcm') || strcmp(files(i).name(end-2:end),'IMA')) % check if element is a directory or a non-DICOM file
         %header = FastReadDICOMHeader(strcat(path_source,files(i).name));
        header = dicominfo(strcat(path_source,files(i).name));
        
        % remove unallowed characters from series description, if present
        forwardslashes = strfind(header.SeriesDescription,'/');
        backslashes = strfind(header.SeriesDescription,'\');
        colons = strfind(header.SeriesDescription,':');
        asterisks = strfind(header.SeriesDescription,'*');
        questionmarks = strfind(header.SeriesDescription,'?');
        doublequotes = strfind(header.SeriesDescription,'"');
        lessthans = strfind(header.SeriesDescription,'<');
        greaterthans = strfind(header.SeriesDescription,'>');
        pipes = strfind(header.SeriesDescription,'|');
        bad_chars = [forwardslashes;backslashes;colons;asterisks;questionmarks;doublequotes;lessthans;greaterthans;pipes];
        for j = 1:length(bad_chars)
            header.SeriesDescription(bad_chars(j)) = '~';
        end
        
        if ~isdir(strcat(path_target,int2str(header.SeriesNumber),'_',header.SeriesDescription)) % check if series directory already exists
            mkdir([path_target,strcat(int2str(header.SeriesNumber),'_',header.SeriesDescription)])
        end
        
        if strcmp(convention,'eFilm')
            % fastcopyfile(strcat(path_source,files(i).name),strcat(path_target,int2str(header.SeriesNumber),'_',header.SeriesDescription,'/',int2str(header.SeriesNumber),'_',int2str(header.InstanceNumber-1),'.dcm'));
            copyfile(strcat(path_source,files(i).name),strcat(path_target,int2str(header.SeriesNumber),'_',header.SeriesDescription,'/',int2str(header.SeriesNumber),'_',int2str(header.InstanceNumber-1),'.dcm'));
        elseif strcmp(convention,'random')
            % fastcopyfile(strcat(path_source,files(i).name),strcat(path_target,int2str(header.SeriesNumber),'_',header.SeriesDescription,'/',int2str(count),'.dcm'));
            copyfile(strcat(path_source,files(i).name),strcat(path_target,int2str(header.SeriesNumber),'_',header.SeriesDescription,'/',int2str(count),'.dcm'));
            
            count = count+1;
        else % strcmp(convention,'standard')
            % fastcopyfile(strcat(path_source,files(i).name),strcat(path_target,int2str(header.SeriesNumber),'_',header.SeriesDescription,'/',int2str(header.InstanceNumber),'.dcm'));
            copyfile(strcat(path_source,files(i).name),strcat(path_target,int2str(header.SeriesNumber),'_',header.SeriesDescription,'/',int2str(header.InstanceNumber),'.dcm'),'f');
        end
    end
    
    ProgressBar(handle_progress,i/length(files))
end

close(handle_progress)


% ----------------------------------------------------------------
function path = AddTrailingSlash(path)

if (path(end) ~= '/') & (path(end) ~= '\')
    path = [path '/'];
end