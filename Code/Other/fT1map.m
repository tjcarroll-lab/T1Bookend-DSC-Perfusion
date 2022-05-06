function [T1map, M0map, InvFmap, RESNORMmap,signalStorage, signalFitStorage] = fT1map(path, SF, mask);

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
images  = Cells2Matrix(IA_ReadImages(path));

if nargin < 3
    mask = automaskns(images);
    if nargin < 2
        SF = [1 1];
    end
end

% SF
S_images = fSpatialFilter(images,SF,mask);

% initial variables
T1map       = zeros(size(S_images,1),size(S_images,2));
M0map       = zeros(size(S_images,1),size(S_images,2));
InvFmap     = zeros(size(S_images,1),size(S_images,2));
RESNORMmap     = zeros(size(S_images,1),size(S_images,2));
% info = dicominfo([path '\1.dcm']);
info = dicominfo([path '\1.dcm'],'Dictionary','dicom-dict.txt'); %Upon the upgrade of the dicom dictionary

% make time tables
if strfind(info.SequenceName,'tfiti')
    T = TimeTableforIRtfiti(path,size(S_images,3),'20linear');
elseif strfind(info.SequenceName,'LL') 
    T = TimeTableforIRLLEPI(path,size(S_images,3));
elseif strfind(info.ProtocolName,'LL') %Jessy (10/2/2007): added the 'ProtocolName' option so that it is compatible with VB13A scans headers 
    T = TimeTableforIRLLEPI(path,size(S_images,3));
end

% calculate & show remaining time 
handle_progress = ProgressBar('Calculating T1 map alone');
Ncount = sum(sum(mask));
count  = 0;

%fitting
signalStorage = zeros(size(S_images));
signalFitStorage = zeros(size(S_images));
if strfind(info.SequenceName,'LL')  
    for i = 1:size(S_images,1)
        for j = 1:size(S_images,2)
            if mask(i,j)
                signal = double(squeeze(S_images(i,j,:))');
                signalStorage(i,j,:) = S_images(i,j,:);
                [T1,M0,InvF,RESNORM,signalFit] = fIR_LLEPI_fitting(T,signal);
                T1map(i,j)   = T1;
                M0map(i,j)   = M0;
                InvFmap(i,j) = InvF;
                RESNORMmap(i,j) = RESNORM;
                signalFitStorage(i,j,:) = signalFit(:);
            end
        end
        count = count + sum(sum(mask(i,:)));
        ProgressBar(handle_progress,count/Ncount);
    end
    
elseif strfind(info.ProtocolName,'LL') %Jessy (10/2/2007): added the 'ProtocolName' option
        for i = 1:size(S_images,1)
        for j = 1:size(S_images,2)
            if mask(i,j)
                if i == 113 && j == 63
                    1;
                end
                signal = double(squeeze(S_images(i,j,:))');
                signalStorage(i,j,:) = S_images(i,j,:);
                [T1,M0,InvF,RESNORM,signalFit] = fIR_LLEPI_fitting(T,signal);
                T1map(i,j)   = T1;
                M0map(i,j)   = M0;
                InvFmap(i,j) = InvF;
                RESNORMmap(i,j) = RESNORM;
                signalFitStorage(i,j,:) = signalFit();
            end
        end
        count = count + sum(sum(mask(i,:)));
        ProgressBar(handle_progress,count/Ncount);
        end
    
elseif strfind(info.SequenceName,'tfiti')
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
        ProgressBar(handle_progress,count/Ncount);
    end
end
close(handle_progress)

clear S_images images mask T signal

