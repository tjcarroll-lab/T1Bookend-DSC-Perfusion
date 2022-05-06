function [T1map, M0map, InvFmap, RESNORMmap] = ucdog_T1map(T1dir, N_Slices, SF, mask)

% Adapted from ncp_T1wholeIm
%
% Author: Yong Ik Jeong
% Date: 2016-03-16


%adapted from fT1map_Bookend
%Author: NRC 2013.10.31



%ORIGNAL HELP TEXT:
% [T1map, M0map, InvFmap] = fT1map(path, SF, mask)
%
% T1 map from fitting
% author:   Wanyong Shin
%
% description: calculate T1 values pixel by pixel by fitting.  fitting
% algorighm is set to fitting up to null point in case of tfiTI, and
% fitting whole curve in case of LL-EPI.
%
% references
%
% status:   stable
%
% versions
%   [05-08-19] (WYS): initial version
%   [07-10-02] (JJM): include small adaptation to VB13A scans headers

% read images fast & filtering
images  = Cells2Matrix(IA_ReadImages(T1dir));

%masking
if nargin < 4 %Jessy (6/25/2008)
    %YIJ 20170614 smooth mask
    mask = automaskns(imfilter(images,fspecial('gaussian',[3 3],2)));
    mask = imfill(mask);
end

%spatial smoothing
if nargin < 3 %Jessy (6/25/2008)
    SF = [1 1];
end

%YIJ 20170614 smooth
%images = imfilter(images,fspecial('gaussian',[3 3],1));

% SF
S_images = fSpatialFilter(images,SF,mask);

% initialize variables
T1map       = zeros(size(S_images,1),size(S_images,2));
M0map       = zeros(size(S_images,1),size(S_images,2));
InvFmap     = zeros(size(S_images,1),size(S_images,2));
RESNORMmap     = zeros(size(S_images,1),size(S_images,2));
info = dicominfo(fullfile(T1dir,'1.dcm'));

% make time tables
scantype = getscantype(info);
switch scantype
    case 'tifti'
        T = TimeTableforIRtfiti(T1dir,size(S_images,3),'20linear');
    case 'LL'
        T = ucdog_TimeTableforIRLLEPI(T1dir,size(S_images,3));
    case 'ep2d_perf'
        LLEPIfiles = dir(T1dir);
        NTimePoints = size(LLEPIfiles,1)-2;
%        NTimePoints = size(nc_lsdcm(T1dir),1); %other method will give wrong answer if non dcm files are in the directory
        T = TimeTableforIRLLEPI_Bookend(T1dir,NTimePoints,N_Slices);
    otherwise
        disp(['Sequence name: ' info.SequenceName])
        disp(['Protocol name: ' info.ProtocolName])
        error('could not identify sequence type for T1 scan')
end



% calculate & show remaining time
%handle_progress = ProgressBar('Calculating T1 map alone');
Ncount = sum(sum(mask));
count  = 0;

%fitting
switch scantype
    case {'LL','ep2d_perf'}
        %YIJ 20170616 parallel computing
        jsize = size(S_images,2);
%         if ~matlabpool('size')
%             matlabpool
%         end
        %parfor i = 1:size(S_images,1)
        for i = 1:size(S_images,1)
            for j = 1:jsize
                if mask(i,j)
%                     if i == 143 && j == 100
%                         disp('hi');
%                     end
                    signal = double(squeeze(S_images(i,j,:))');
                    %YIJ 20170615 smooth
                    signal = smooth(signal,5,'sgolay',3)';
                    [T1,M0,InvF,RESNORM] = fIR_LLEPI_fitting(T,signal);
                    T1map(i,j)   = T1;
                    M0map(i,j)   = M0;
                    InvFmap(i,j) = InvF;
                    RESNORMmap(i,j) = RESNORM;
                end
            end
            %count = count + sum(sum(mask(i,:)));
            %ProgressBar(handle_progress,count/Ncount);
        end
    case 'tifti'
        for i = 1:size(S_images,1)
            for j = 1:size(S_images,2)
                if mask(i,j)
                    signal = double(squeeze(S_images(i,j,:))');
                    [T1,M0,RESNORM] = fIR_tfiti_fittinguptoNull(T,signal);
                    T1map(i,j)   = T1;
                    M0map(i,j)   = M0;
                    InvFmap(i,j)  = 2;
                    RESNORMmap(i,j) = RESNORM;
                end
            end
            count = count + sum(sum(mask(i,:)));
            %ProgressBar(handle_progress,count/Ncount);
        end
end

%close(handle_progress)

end


function scantype = getscantype(info)

%if strfind(info.SequenceName,'tfiti')
%    scantype='tifti';
%elseif strfind(info.SequenceName,'LL')
%    scantype = 'LL';
if strfind(info.ProtocolName,'LL') %Jessy (10/2/2007): added the 'ProtocolName' option so that it is compatible with VB13A scans headers
    scantype = 'LL';
elseif strfind(info.ProtocolName,'ep2d_gre_ScalePWI') %Jessy (10/2/2007): added the 'ProtocolName' option so that it is compatible with VB13A scans headers
    scantype = 'ep2d_perf';
elseif strfind(info.ProtocolName,'ep2d_ScalePWI_674') %Jessy (10/2/2007): added the 'ProtocolName' option so that it is compatible with VB13A scans headers
    scantype = 'ep2d_perf';
elseif strfind(info.ProtocolName,'Bookend') %Jessy (6/25/2008): new Protocol Name for Bookend (single sequence)
    scantype = 'ep2d_perf';
elseif strfind(info.ProtocolName,'ep2d_perf_PWI') %Neil (07/08/2013): to make this compatible with RAPID protocol stuff.  Copied from 'ep2d_ScalePWI_674'
    scantype = 'ep2d_perf';
elseif strfind(info.ProtocolName,'ep2d_perf') %Neil (07/08/2013): to make this compatible with RAPID protocol stuff.  Copied from 'ep2d_ScalePWI_674'
    scantype = 'ep2d_perf';
elseif strfind(info.ProtocolName,'ep2d_AX PERF ScalePWI_674') %Neil (07/08/2013): to make this compatible with RAPID protocol stuff.  Copied from 'ep2d_ScalePWI_674'
    scantype = 'ep2d_perf';
elseif strfind(info.ProtocolName,'ep2d_ScalePWI') %Neil (10/01/2013): yet another entry
    scantype = 'ep2d_perf';
end

end
