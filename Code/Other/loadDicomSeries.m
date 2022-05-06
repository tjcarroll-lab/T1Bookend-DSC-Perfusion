function [data] = loadDicomSeries(folderpath)
%loadDicomSeries Load dicom series from folder
%   For COLLAT
%   Inputs:
%       - folderpath: path to folder with dicoms
%       - empty inputs to select folder via UI
%
%   Outputs:
%       - data: nx by ny by ntime by nslice matrix
%
% Author: Yong Ik Jeong
% Date: 2020-06-29
% Changelog:
%   - 20200629 YIJ: 

if nargin < 1
    folderpath = uigetdir();
end
currdatetime = datestr(datetime('now'),'yyyymmddHHMMSS');

% If zip, unzip
[filepath,name,ext] = fileparts(folderpath);
if strcmpi(ext,'.zip')
    % Copy contents to temporary folder first
    tmpfolder = fullfile(filepath,['tmp_' currdatetime]);
    unzipfiles = unzip(folderpath,tmpfolder);
    tmpdir = dir(tmpfolder);
    % Check if content is a folder or files
    if tmpdir(end).isdir
        movefile(fullfile(tmpdir(end).folder,tmpdir(end).name),filepath);
        name = tmpdir(end).name;
    else
        if ~exist(fullfile(filepath,name),'dir')
            mkdir(fullfile(filepath,name));
        end
        for ii = 1:length(unzipfiles)
            movefile(unzipfiles{ii},fullfile(filepath,name));
        end
    end
    rmdir(tmpfolder);
end

% Sort dicom file names
dcmdir = dir(fullfile(filepath,name,'\*.dcm'));
dcmnames = {dcmdir.name}';
[dcmnames,sortind] = sortFiles(dcmnames,'_.');
dcmdir = dcmdir(sortind);

% Get number of time points and slices
tmpimg = dicomread(fullfile(filepath,name,dcmdir(1).name));
tmphdr = dicominfo(fullfile(filepath,name,dcmdir(1).name));
nTP = tmphdr.NumberOfTemporalPositions;
if mod(length(dcmdir),nTP) == 0
    nSlc = length(dcmdir)/nTP;
else
    error('Number of slices and timepoints do not fall off even.');
end

% Load dicoms into matrix
data = zeros(size(tmpimg,1),size(tmpimg,2),length(dcmdir));
for ii = 1:length(dcmdir)
    data(:,:,ii) = dicomread(fullfile(filepath,name,dcmdir(ii).name));
end
data = reshape(data,size(data,1),size(data,2),nTP,nSlc);

end

