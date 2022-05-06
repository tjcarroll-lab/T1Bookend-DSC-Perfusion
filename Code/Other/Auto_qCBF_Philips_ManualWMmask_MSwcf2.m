function Auto_qCBF_Philips_ManualWMmask_MSwcf2(varargin)
%-------------------------------------------------------------------------%
% Format: Auto_qCBF([ID1 ID2 ID3 ...])
%
% Inputs:
% [ID1 ID2 ID3 ...]: The ID numbers for each patient
%
% Outputs:
% T1maps, AutoSS, AIFdata, DSCanalysis, Result files, and qCBV, qCBF, and
%   CMTT maps converted to DICOM format, saved in corresponding
%   directories
%
% This is the full file for automatic qCBF calculation revised to make
% targetpath current directory (no more valid - JJM)

%-------------------------------------------------------------------------%
% Run from main directory of study which is specified by the "targetpath"
% and which  must contain the following directories: T1mapping, AutoSS,
% AIFdata, DSCanalysis, Result, DICOM (including the different map types
% inside); don't forget to specify the full "targetpath" and the
% "Seqtarget"
%-------------------------------------------------------------------------%

% Author: Maulin Shah
%   Combined WanYong's code into one main file
% Date: 01/30/2007
%
% Modified by Jessy Mouannes:
%   07/24/2007:
%       - added DICOM conversion code using "autoWriteDICOM" to convert final
%           quantitation maps to DICOM format
%       - "targetpath" can no more be specified as the current directory
%   07/25/2007:
%       - "Seqtarget" is now only specified once in the main code section
%   08/03/2007:
%       - Different series numbers are given to different map series
%   08/06/2007
%       - Different "factor" value is given to different map series
%-------------------------------------------------------------------------%

%% Main Code
% targetpath = '../Bookend(LLEPI)';   % Enter targetpath here
% targetpath = 'C:\Documents and Settings\CWloka\My Documents\MATLAB\MS Processing\BookEnd Processing\Main Directory';
%targetpath = 'C:\Users\Grady\Dropbox\Toronto_Northwestern\OldCode\Main Directory - modified 6 April 2010'
%targetpath = pwd
%Seqtarget   = varargin{2};%'GE';    %'GE', 'SE', or 'all'
%global N_meas N_slices

global injectionNum;
global glblTargetPath;
global seqType;
global sampdelayfix;
global DDfix;
global MANUAL;

if nargin
    sampdelayfix = 1;
    DDfix = 0;
    MANUAL = 0;
    
    targetpath = varargin{1};
    glblTargetPath = targetpath;
    %input = varargin{1};
    input = dir([targetpath '/P*']);
    for ii = 1:length(input)
        t1 = tic;
        %PN = makePN(input(i))
        PN = input(ii).name;
        slclist = [];
        % YIJ: 20180118 Do only first set (for 2nd day data)
        if sum(strcmpi(PN,{'P001'}))%sum(strcmpi(PN,{'P001','P002','P003'}))
            tmpnum = PN(find(~ismember(PN,'P0'),1):end);
            injectionNum = str2double(tmpnum);
            fepath = fullfile(targetpath,PN,'ep2d_perf');
            sepath = fullfile(targetpath,PN,'ep2d_se');
            if exist(fepath,'dir')
                Seqtarget = 'GE';
                perfpath = fepath;
                seqType = 'GE_M';
            elseif exist(sepath,'dir')
                Seqtarget = 'SE';
                perfpath = sepath;
                seqType = 'SE_M';
            else
                error('blah');
            end
            perfdir = dir([perfpath '\*.dcm']);
            header = dicominfo([perfpath '\1.dcm']);
            
            if isfield(header,'NumberOfTemporalPositions')
                N_meas = header.NumberOfTemporalPositions;
                N_slices = length(perfdir)/N_meas;
            else
                for jj = 1:length(perfdir)
                    header = dicominfo([perfpath '\' num2str(jj) '.dcm']);
                    slclist(jj) = header.SliceLocation;
                end
                N_slices = length(unique(slclist));
                N_meas = length(perfdir)/N_slices;
            end
                          
            aifpath = fullfile(targetpath,'AIFdata');
            sspath = fullfile(targetpath,'AutoSS');
            dscpath = fullfile(targetpath,'DSCanalysis');
            resultpath = fullfile(targetpath,'Result_MSwcf2');
            t1mappath = fullfile(targetpath,'T1mapping');
            if ~exist(aifpath,'dir')
                mkdir(aifpath);
            end
            if ~exist(sspath,'dir')
                mkdir(sspath);
            end
            if ~exist(dscpath,'dir')
                mkdir(dscpath);
            end
            if ~exist(resultpath,'dir')
                mkdir(resultpath);
            end
            if ~exist(t1mappath,'dir')
                mkdir(t1mappath);
            end
            
            t1mapping(targetpath,PN);
            autoss(targetpath,PN);
        autoaif(targetpath,PN,Seqtarget,N_slices,N_meas);
        deconvolution(targetpath,PN,Seqtarget,N_slices,N_meas);
        autoquantification(targetpath,PN,Seqtarget,N_slices,N_meas);
        end
        toc(t1);
        %autoWriteDICOM(targetpath,PN,Seqtarget)
    end
end

fprintf('%s\n','All is done');

return

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%% T1 Mapping
function t1mapping(varargin)

targetpath = varargin{1};

input = varargin{2};

if (strfind(input,'P') == 1) & (length(input)==4)
    PN = input;
    fprintf('---------------------------------------------------------\n');
    tic
    loopingT1maps(targetpath, PN)
    toc
    fprintf('---------------------------------------------------------\n');
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loopingT1maps(varargin)

targetpath = varargin{1};
opt = varargin{2};

files    = dir(targetpath);
T1Matfiles = dir([targetpath '/T1mapping']);

% selective patients to run
PN = opt;
targetT1file    = [PN '_T1map.mat'];
if namefindinstrcut(files,PN)
    fprintf(['Checking T1map for ' PN ' : \n']);
    if namefindinstrcut(T1Matfiles,targetT1file,'exact')
        fprintf('%s\n',[targetpath '/' PN ' T1 mapping has been written...']);
    else
        fprintf('%s\n',['T1 mapping for patient ' targetpath '/' PN '...']);
        runitforT1mapping(targetpath,PN);
    end
    fprintf('T1 mapping done\n');
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function runitforT1mapping(targetpath,PN)

pathPN = [targetpath '/' PN];
Seqfiles = dir(pathPN);
count = 1;index=0;

for j = 1:length(Seqfiles);
    temp = Seqfiles(j);
    if ~isempty(strfind(temp.name,'LL_EPI_PRE')) |  ~isempty(strfind(temp.name,'LL_EPI_pre'))
        if isempty(strfind(temp.name,'LL_EPI_PRE2'))
            IRpathpre = [pathPN '/' temp.name];
            [T1map_pre,M0map_pre,InvFmap_pre,RESNORMmap_pre] = ucdog_T1map(IRpathpre);
            image_names{count}   = 'T1map_pre'; image_names{count+1} = 'M0map_pre'; image_names{count+2} = 'InvFmap_pre';   image_names{count+3} = 'RESNORMmap_pre';
            images{count}=T1map_pre;            images{count+1}=M0map_pre;          images{count+2}=InvFmap_pre;            images{count+3}=RESNORMmap_pre;
            count = count+4; index=1;
        end
    end
    if  ~isempty(strfind(temp.name,'LL_EPI_POST')) |  ~isempty(strfind(temp.name,'LL_EPI_post'))
        if isempty(strfind(temp.name,'LL_EPI_POST2'))
            IRpathpost = [pathPN '/' temp.name];
            [T1map_post,M0map_post,InvFmap_post,RESNORMmap_post] = ucdog_T1map(IRpathpost);
            image_names{count}   = 'T1map_post'; image_names{count+1} = 'M0map_post'; image_names{count+2} = 'InvFmap_post';    image_names{count+3} = 'RESNORMmap_post';
            images{count}=T1map_post;            images{count+1}=M0map_post;          images{count+2}=InvFmap_post;             images{count+3}=RESNORMmap_post;
            count = count+4; index=1;
        end
    end
end
if index
    save([targetpath '/T1mapping/' PN '_T1map.mat'], 'images', 'image_names');
    %save([targetpath '/' PN '_fitting.mat'], 'signalStoragePre', 'signalStoragePost', 'signalFitPre', 'signalFitPost');
end

return

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%% Auto SS
function autoss(varargin)

targetpath = varargin{1};

input = varargin{2};

if (strfind(input,'P') == 1) & (length(input)==4)
    PN = input;
    fprintf('---------------------------------------------------------\n');
    tic
    loopingautoss(targetpath, PN);
    toc
    fprintf('---------------------------------------------------------\n');
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loopingautoss(varargin)

targetpath = varargin{1};
opt = varargin{2};

files      = dir(targetpath);
T1Matfiles = dir([targetpath '/T1mapping']);
SSfiles    = dir([targetpath '/AutoSS']);

% selective patients to run
PN = opt;
targetT1file = [PN '_T1map.mat'];
targetSSfile    = [PN '_AutoSS.mat'];
if namefindinstrcut(files,PN,'exact')
    fprintf(['Checking AutoSS for ' PN ' : \n']);
    if namefindinstrcut(T1Matfiles,targetT1file,'exact')
        if namefindinstrcut(SSfiles,targetSSfile,'exact')
            fprintf('%s\n',[targetpath '/' PN ' autoSS has been written...']);
        else
            fprintf('%s\n',['autoSS for patient ' targetpath '/' PN '...']);
            runitforautoss(targetpath,targetT1file);          
        end
        fprintf('AutoSS done\n');
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function runitforautoss (targetpath,targetT1file)

load([targetpath '/T1mapping/' targetT1file]);
if ~isempty(strmatch('T1map_pre',image_names)) & ~isempty(strmatch('M0map_pre',image_names)) & ~isempty(strmatch('T1map_post',image_names)) & ~isempty(strmatch('M0map_post',image_names))
    T1map_pre  = images{strmatch('T1map_pre',image_names)};
    T1map_post = images{strmatch('T1map_post',image_names)};
    M0map_pre  = images{strmatch('M0map_pre',image_names)};
    M0map_post = images{strmatch('M0map_post',image_names)};
    RESNORMmap_pre = images{strmatch('RESNORMmap_pre',image_names)};
    RESNORMmap_post = images{strmatch('RESNORMmap_post',image_names)};
    [position_SS,SSmask] = AutoSS(T1map_pre,T1map_post,M0map_pre,M0map_post,RESNORMmap_pre,RESNORMmap_post);
    ROIs.positions.blood_pre = position_SS;
    ROIs.positions.blood_post = position_SS;
    masks_ROIs.SS = SSmask;
    save([targetpath '/AutoSS/' targetT1file(1:5) 'AutoSS.mat'], 'ROIs', 'masks_ROIs','images','image_names')
end

return

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%% AIF
function autoaif(varargin)

targetpath  = varargin{1};
PNinput     = varargin{2};
Seqtarget   = varargin{3};
option      = 'Skip';
N_slices = varargin{4};
N_meas = varargin{5};

fprintf('---------------------------------------------------------\n');
tic
loopingautoaif(targetpath,Seqtarget,PNinput,option,N_slices,N_meas);
toc
fprintf('---------------------------------------------------------\n');

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loopingautoaif(targetpath,Seqtarget,PNinput,option,N_slices,N_meas)

files       = dir(targetpath);
AIFfiles    = dir([targetpath '/AIFdata']);

runitforautoaif(targetpath,PNinput,Seqtarget,option,N_slices,N_meas);
fprintf('AutoAIF done\n');

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function runitforautoaif(targetpath,PN,Seqtarget,option,N_slices,N_meas)

pathPN = [targetpath '/' PN];
Seqfiles = dir(pathPN);
AIFfiles    = dir([targetpath '/AIFdata']);

if strmatch(Seqtarget,'GE')
    targetAIFfile = [PN 'GE_AIF_A.mat'];
    if namefindinstrcut(Seqfiles,'ep2d_perf')
        targetpathDSC = [pathPN '/' Seqfiles(namefindinstrcut(Seqfiles,'ep2d_perf')).name];
        fprintf('%s\n',['Checking AutoAIF for ' PN ' : ']);
        if namefindinstrcut(AIFfiles,targetAIFfile,'exact')
            fprintf('%s\n',[targetpath '\' targetAIFfile ' has been created...']);
            if strmatch(option,'Overwrite')
                runitforautoaifcore(targetpath, targetpathDSC,targetAIFfile,N_slices,N_meas);
            end
        else
            fprintf('%s\n',[targetpath '\' targetAIFfile ' is being created now...']);
            runitforautoaifcore(targetpath, targetpathDSC,targetAIFfile,N_slices,N_meas);
        end
    end
elseif strmatch(Seqtarget,'SE')
    targetAIFfile = [PN 'SE_AIF_A.mat'];
    if namefindinstrcut(Seqfiles,'ep2d_se')
        targetpathDSC = [pathPN '/' Seqfiles(namefindinstrcut(Seqfiles,'ep2d_se')).name];
        fprintf('%s\n',['Checking AutoAIF for ' PN ' : ']);
        if namefindinstrcut(AIFfiles,targetAIFfile,'exact')
            fprintf('%s\n',[targetpath '\' targetAIFfile ' has been created...']);
            if strmatch(option,'Overwrite')
                runitforautoaifcore(targetpath, targetpathDSC,targetAIFfile,N_slices,N_meas);
            end
        else
            fprintf('%s\n',[targetpath '\' targetAIFfile ' is being created now...']);
            runitforautoaifcore(targetpath, targetpathDSC,targetAIFfile,N_slices,N_meas);
        end
    end
elseif strmatch(Seqtarget,'all')
    targetAIFfile = [PN 'GE_AIF_A.mat'];
    if namefindinstrcut(Seqfiles,'ep2d_perf')
        targetpathDSC = [pathPN '/' Seqfiles(namefindinstrcut(Seqfiles,'ep2d_perf')).name];
        fprintf('%s\n',['Checking AutoAIF for ' PN ' : ']);
        if namefindinstrcut(AIFfiles,targetAIFfile,'exact')
            fprintf('%s\n',[targetpath '\' targetAIFfile ' has been created...']);
            if strmatch(option,'Overwrite')
                runitforautoaifcore(targetpath, targetpathDSC,targetAIFfile,N_slices,N_meas);
            end
        else
            fprintf('%s\n',[targetpath '\' targetAIFfile ' is being created now...']);
            runitforautoaifcore(targetpath, targetpathDSC,targetAIFfile,N_slices,N_meas);
        end
    end
    targetAIFfile = [PN 'SE_AIF_A.mat'];
    if namefindinstrcut(Seqfiles,'ep2d_se')
        targetpathDSC = [pathPN '/' Seqfiles(namefindinstrcut(Seqfiles,'ep2d_se')).name];
        fprintf('%s\n',['Checking AutoAIF for ' PN ' : ']);
        if namefindinstrcut(AIFfiles,targetAIFfile,'exact')
            fprintf('%s\n',[targetpath '\' targetAIFfile ' has been created...']);
            if strmatch(option,'Overwrite')
                runitforautoaifcore(targetpath, targetpathDSC,targetAIFfile,N_slices,N_meas);
            end
        else
            fprintf('%s\n',[targetpath '\' targetAIFfile ' is being created now...']);
            runitforautoaifcore(targetpath, targetpathDSC,targetAIFfile,N_slices,N_meas);
        end
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function runitforautoaifcore(targetpath,targetpathDSC,targetAIFfile,N_slices,N_meas)
% YIJ 20180129: autoaif_wy_Philips_v2 used for data in paper
global injectionNum;
if injectionNum == 0 %|| injectionNum == 1
    [targetpath '\AIFdata\' 'P001GE_AIF_A.mat']
    tmpdata = load([targetpath '\AIFdata\' 'P001GE_AIF_A.mat']);
    %tmpdata = [];
    [signal_AIF,cutoffs_AIF,position_AIF,conc_AIF,dsc_stack,n_slice_AIF,AIF,position_slice_AIF,aifmask] = autoaif_wy_Philips_v2([targetpathDSC],N_slices,N_meas,tmpdata);
else
    tmpdata = [];
    [signal_AIF,cutoffs_AIF,position_AIF,conc_AIF,dsc_stack,n_slice_AIF,AIF,position_slice_AIF,aifmask] = autoaif_wy_Philips_v2([targetpathDSC],N_slices,N_meas,tmpdata);
end

%[signal_AIF,cutoffs_AIF,position_AIF,conc_AIF,dsc_stack,n_slice_AIF,AIF,position_slice_AIF] = autoaif_wy_Philips([targetpathDSC],N_slices,N_meas);
if  (size(signal_AIF,1) && size(cutoffs_AIF,1) && size(position_AIF,1))
    ROIs.dsc_stack = dsc_stack;
    ROIs.data.AIF = signal_AIF;
    ROIs.data.conc = conc_AIF;
    ROIs.positions.AIF = position_AIF;
    ROIs.positions.n_slice_AIF = n_slice_AIF;
    ROIs.data.AIFslice = AIF;
    ROIs.positions.AIFslice = position_slice_AIF;
    cutoffs_ROIs.AIF = cutoffs_AIF;
    masks_ROIs.aifmask = aifmask;
    if(cutoffs_ROIs.AIF(1) == 0)
        cutoffs_ROIs.AIF(1) = 1;
    end
    save([targetpath '/AIFdata/' targetAIFfile], 'ROIs', 'cutoffs_ROIs','masks_ROIs');
    %YIJ 20170614
    %figure;plot(ROIs.data.AIF);
else
    fprintf([targetpath '/' targetAIFfile ' cannot be created because no contrast agent was injected.\n']);
end

return

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%% Deconvolution
function deconvolution(varargin)

targetpath  = varargin{1};
PNinput     = varargin{2};
Seqtarget   = varargin{3};
option      = 'Skip';
N_slices = varargin{4};
N_meas = varargin{5};

fprintf('---------------------------------------------------------\n');
tic
loopingdeconv(targetpath,Seqtarget,PNinput,option,N_slices,N_meas)
toc
fprintf('---------------------------------------------------------\n');

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loopingdeconv(targetpath,Seqtarget,PNinput,option,N_slices,N_meas)

files       = dir(targetpath);
AIFfiles    = dir([targetpath '/AIFdata']);
DSCfiles    = dir([targetpath '/DSCanalysis']);

runitfordeconv(targetpath,PNinput,Seqtarget,option,N_slices,N_meas);
fprintf('Deconv done\n');

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function runitfordeconv(targetpath,PN,Seqtarget,option,N_slices,N_meas)

pathPN = [targetpath '/' PN];
Seqfiles = dir(pathPN);
AIFfiles    = dir([targetpath '/AIFdata']);
DSCfiles    = dir([targetpath '/DSCanalysis']);
global injectionNum;
global DDfix;
tmpind = 0;
tmpinj = 'P001';
if strmatch(Seqtarget,'GE')
    targetAIFfile = [PN 'GE_AIF_A.mat'];
    if injectionNum == tmpind
        targetAIFfile = [tmpinj 'GE_AIF_A.mat']
    end
    targetDSCfile = [PN 'GE_DSC_A.mat'];
    if namefindinstrcut(Seqfiles,'ep2d_perf')
        targetpathDSC = [pathPN '/' Seqfiles(namefindinstrcut(Seqfiles,'ep2d_perf')).name];
        fprintf('%s\n',['Checking Deconv for ' PN ' : ']);
        if namefindinstrcut(AIFfiles,targetAIFfile,'exact') & namefindinstrcut(DSCfiles,targetDSCfile,'exact')
            fprintf('%s\n',[targetpath '\' targetDSCfile ' has been created...']);
            if strmatch(option,'Overwrite')
                fprintf('%s\n',[targetpath '\' targetDSCfile ' is being overwritten...']);
                %%
                images=cell(0);image_names=cell(0);
                load ([targetpath '\AIFdata\' targetAIFfile]);
                %                [DSC,headers,N_images] = ReadDICOMSeries(targetpathDSC,1);
                %                N_images = 15;
                N_images = N_slices; %Jessy (10/31/2008) %25 is correct Grady
                [CBV_DSC,CBF_nSVD,CMTT_nSVD,TmaxMap_nSVD,dBATmap,Conc_map,CmtWM] = Calc_nSVD_Philips(targetpathDSC,N_images,ROIs.data.AIF,cutoffs_ROIs.AIF,N_meas,ROIs.data.AIFslice,ROIs.positions.n_slice_AIF);
                cutoffs_ROIs.CmtWM = CmtWM;
                [images, image_names] = IA_ImportImages(images, image_names,'CBV_DSC',CBV_DSC,'CBF_nSVD',CBF_nSVD,'CMTT_nSVD',CMTT_nSVD,'TmaxMap_nSVD',TmaxMap_nSVD,'dBATmap',dBATmap,'Conc_map',Conc_map);
                save([targetpath '/DSCanalysis/' targetDSCfile], 'ROIs','cutoffs_ROIs','images', 'image_names');
                %%
            end
        elseif namefindinstrcut(AIFfiles,targetAIFfile,'exact')
            fprintf('%s\n',[targetpath '\' targetDSCfile ' is being created now...']);
            %%
            images=cell(0);image_names=cell(0);
            load ([targetpath '\AIFdata\' targetAIFfile]);
            %            [DSC,headers,N_images] = ReadDICOMSeries(targetpathDSC,1);
            %            N_images = 15;
            N_images = N_slices; %Jessy (10/31/2008) %25 is correct Grady
            [CBV_DSC,CBF_nSVD,CMTT_nSVD,TmaxMap_nSVD,dBATmap,Conc_map,CmtWM,ATDmap,IBATmap,BATmap,BRTmap,sigmap,VEINslice,veinmask,positionVein,n_Veinslice] = Calc_nSVD_Philips(targetpathDSC,N_images,ROIs.data.AIF,cutoffs_ROIs.AIF,N_meas,ROIs.data.AIFslice,ROIs.positions.n_slice_AIF);
            cutoffs_ROIs.CmtWM = CmtWM;
            ROIs.sigmap = sigmap;
            ROIs.data.VEINslice = VEINslice;
            masks_ROIs.veinmask = veinmask;
            ROIs.values.DDfix = DDfix;
            ROIs.positions.VEINslice = positionVein;
            ROIs.positions.n_slice_VEIN = n_Veinslice;
            [images, image_names] = IA_ImportImages(images, image_names,'CBV_DSC',CBV_DSC,'CBF_nSVD',CBF_nSVD,'CMTT_nSVD',CMTT_nSVD,'TmaxMap_nSVD',TmaxMap_nSVD,'dBATmap',dBATmap,'Conc_map',Conc_map,'ATDmap',ATDmap,'IBATmap',IBATmap,'BATmap',BATmap,'BRTmap',BRTmap);
            save([targetpath '/DSCanalysis/' targetDSCfile], 'ROIs','cutoffs_ROIs','images', 'image_names','masks_ROIs');
            %%
        end
    end
elseif strmatch(Seqtarget,'SE')
    targetAIFfile = [PN 'SE_AIF_A.mat'];
    if injectionNum == tmpind
        targetAIFfile = [tmpinj 'SE_AIF_A.mat'];
    end
    targetDSCfile = [PN 'SE_DSC_A.mat'];
    if namefindinstrcut(Seqfiles,'ep2d_se')
        targetpathDSC = [pathPN '/' Seqfiles(namefindinstrcut(Seqfiles,'ep2d_se')).name];
        fprintf('%s\n',['Checking Deconv for ' PN ' : ']);
        if namefindinstrcut(AIFfiles,targetAIFfile,'exact') & namefindinstrcut(DSCfiles,targetDSCfile,'exact')
            fprintf('%s\n',[targetpath '\' targetDSCfile ' has been created...']);
            if strmatch(option,'Overwrite')
                fprintf('%s\n',[targetpath '\' targetDSCfile ' is being overwritten...']);
                %%
                images=cell(0);image_names=cell(0);
                load ([targetpath '\AIFdata\' targetAIFfile]);
                %[DSC,headers,N_images] = ReadDICOMSeries(targetpathDSC,1);
                N_images = N_slices;
                [CBV_DSC,CBF_nSVD,CMTT_nSVD,TmaxMap_nSVD,dBATmap,Conc_map,CmtWM] = Calc_nSVD_Philips(targetpathDSC,N_images,ROIs.data.AIF,cutoffs_ROIs.AIF,N_meas);
                cutoffs_ROIs.CmtWM = CmtWM;
                [images, image_names] = IA_ImportImages(images, image_names,'CBV_DSC',CBV_DSC,'CBF_nSVD',CBF_nSVD,'CMTT_nSVD',CMTT_nSVD,'TmaxMap_nSVD',TmaxMap_nSVD,'dBATmap',dBATmap,'Conc_map',Conc_map);
                save([targetpath '/DSCanalysis/' targetDSCfile], 'ROIs','cutoffs_ROIs','images', 'image_names');
                %%
            end
        elseif namefindinstrcut(AIFfiles,targetAIFfile,'exact')
            fprintf('%s\n',[targetpath '\' targetDSCfile ' is being created now...']);
            %%
            images=cell(0);image_names=cell(0);
            load ([targetpath '\AIFdata\' targetAIFfile]);
            %[DSC,headers,N_images] = ReadDICOMSeries(targetpathDSC,1);
            N_images = N_slices;
            [CBV_DSC,CBF_nSVD,CMTT_nSVD,TmaxMap_nSVD,dBATmap,Conc_map,CmtWM] = Calc_nSVD_Philips(targetpathDSC,N_images,ROIs.data.AIF,cutoffs_ROIs.AIF,N_meas);
            cutoffs_ROIs.CmtWM = CmtWM;
            [images, image_names] = IA_ImportImages(images, image_names,'CBV_DSC',CBV_DSC,'CBF_nSVD',CBF_nSVD,'CMTT_nSVD',CMTT_nSVD,'TmaxMap_nSVD',TmaxMap_nSVD,'dBATmap',dBATmap,'Conc_map',Conc_map);
            save([targetpath '/DSCanalysis/' targetDSCfile], 'ROIs','cutoffs_ROIs','images', 'image_names');
            %%
        end
    end
elseif strmatch(Seqtarget,'all')
    targetAIFfile = [PN 'GE_AIF_A.mat'];
    if injectionNum == tmpind
        targetAIFfile = [tmpinj 'GE_AIF_A.mat'];
    end
    targetDSCfile = [PN 'GE_DSC_A.mat'];
    if namefindinstrcut(Seqfiles,'ep2d_perf')
        targetpathDSC = [pathPN '/' Seqfiles(namefindinstrcut(Seqfiles,'ep2d_perf')).name];
        fprintf('%s\n',['Checking Deconv for ' PN ' : ']);
        if namefindinstrcut(AIFfiles,targetAIFfile,'exact') & namefindinstrcut(DSCfiles,targetDSCfile,'exact')
            fprintf('%s\n',[targetpath '\' targetDSCfile ' has been created...']);
            if strmatch(option,'Overwrite')
                fprintf('%s\n',[targetpath '\' targetDSCfile ' is being overwritten...']);
                %%
                images=cell(0);image_names=cell(0);
                load ([targetpath '\AIFdata\' targetAIFfile]);
                %[DSC,headers,N_images] = ReadDICOMSeries(targetpathDSC,1);
                N_images = N_slices;
                [CBV_DSC,CBF_nSVD,CMTT_nSVD,TmaxMap_nSVD,dBATmap,Conc_map,CmtWM] = Calc_nSVD_Philips(targetpathDSC,N_images,ROIs.data.AIF,cutoffs_ROIs.AIF,N_meas);
                cutoffs_ROIs.CmtWM = CmtWM;
                [images, image_names] = IA_ImportImages(images, image_names,'CBV_DSC',CBV_DSC,'CBF_nSVD',CBF_nSVD,'CMTT_nSVD',CMTT_nSVD,'TmaxMap_nSVD',TmaxMap_nSVD,'dBATmap',dBATmap,'Conc_map',Conc_map);
                save([targetpath '/DSCanalysis/' targetDSCfile], 'ROIs','cutoffs_ROIs','images', 'image_names');
                %%
            end
        elseif namefindinstrcut(AIFfiles,targetAIFfile,'exact')
            fprintf('%s\n',[targetpath '\' targetDSCfile ' is being created now...']);
            %%
            images=cell(0);image_names=cell(0);
            load ([targetpath '\AIFdata\' targetAIFfile]);
            %[DSC,headers,N_images] = ReadDICOMSeries(targetpathDSC,1);
            N_images = N_slices;
            [CBV_DSC,CBF_nSVD,CMTT_nSVD,TmaxMap_nSVD,dBATmap,Conc_map,CmtWM] = Calc_nSVD_Philips(targetpathDSC,N_images,ROIs.data.AIF,cutoffs_ROIs.AIF,N_meas);
            cutoffs_ROIs.CmtWM = CmtWM;
            [images, image_names] = IA_ImportImages(images, image_names,'CBV_DSC',CBV_DSC,'CBF_nSVD',CBF_nSVD,'CMTT_nSVD',CMTT_nSVD,'TmaxMap_nSVD',TmaxMap_nSVD,'dBATmap',dBATmap,'Conc_map',Conc_map);
            save([targetpath '/DSCanalysis/' targetDSCfile], 'ROIs','cutoffs_ROIs','images', 'image_names');
            %%
        end
    end
    %SE
    targetAIFfile = [PN 'SE_AIF_A.mat'];
    if injectionNum == tmpind
        targetAIFfile = [tmpinj 'SE_AIF_A.mat'];
    end
    targetDSCfile = [PN 'SE_DSC_A.mat'];
    if namefindinstrcut(Seqfiles,'ep2d_se')
        targetpathDSC = [pathPN '/' Seqfiles(namefindinstrcut(Seqfiles,'ep2d_se')).name];
        fprintf('%s\n',['Checking Deconv for ' PN ' : ']);
        if namefindinstrcut(AIFfiles,targetAIFfile,'exact') & namefindinstrcut(DSCfiles,targetDSCfile,'exact')
            fprintf('%s\n',[targetpath '\' targetDSCfile ' has been created...']);
            if strmatch(option,'Overwrite')
                fprintf('%s\n',[targetpath '\' targetDSCfile ' is being overwritten...']);
                %%
                images=cell(0);image_names=cell(0);
                load ([targetpath '\AIFdata\' targetAIFfile]);
                %[DSC,headers,N_images] = ReadDICOMSeries(targetpathDSC,1);
                N_images = N_slices;
                [CBV_DSC,CBF_nSVD,CMTT_nSVD,TmaxMap_nSVD,dBATmap,Conc_map,CmtWM] = Calc_nSVD_Philips(targetpathDSC,N_images,ROIs.data.AIF,cutoffs_ROIs.AIF,N_meas);
                cutoffs_ROIs.CmtWM = CmtWM;
                [images, image_names] = IA_ImportImages(images, image_names,'CBV_DSC',CBV_DSC,'CBF_nSVD',CBF_nSVD,'CMTT_nSVD',CMTT_nSVD,'TmaxMap_nSVD',TmaxMap_nSVD,'dBATmap',dBATmap,'Conc_map',Conc_map);
                save([targetpath '/DSCanalysis/' targetDSCfile], 'ROIs','cutoffs_ROIs','images', 'image_names');
                %%
            end
        elseif namefindinstrcut(AIFfiles,targetAIFfile,'exact')
            fprintf('%s\n',[targetpath '\' targetDSCfile ' is being created now...']);
            %%
            images=cell(0);image_names=cell(0);
            load ([targetpath '\AIFdata\' targetAIFfile]);
            %[DSC,headers,N_images] = ReadDICOMSeries(targetpathDSC,1);
            N_images = N_slices;
            [CBV_DSC,CBF_nSVD,CMTT_nSVD,TmaxMap_nSVD,dBATmap,Conc_map,CmtWM] = Calc_nSVD_Philips(targetpathDSC,N_images,ROIs.data.AIF,cutoffs_ROIs.AIF,N_meas);
            cutoffs_ROIs.CmtWM = CmtWM;
            [images, image_names] = IA_ImportImages(images, image_names,'CBV_DSC',CBV_DSC,'CBF_nSVD',CBF_nSVD,'CMTT_nSVD',CMTT_nSVD,'TmaxMap_nSVD',TmaxMap_nSVD,'dBATmap',dBATmap,'Conc_map',Conc_map);
            save([targetpath '/DSCanalysis/' targetDSCfile], 'ROIs','cutoffs_ROIs','images', 'image_names');
            %%
        end
    end
end

return

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%% Quantification
function autoquantification(varargin)

targetpath  = varargin{1};
PNinput     = varargin{2};
Seqtarget   = varargin{3};
option      = 'Skip';
N_slices = varargin{4};
N_meas = varargin{5};

fprintf('---------------------------------------------------------\n');
tic
loopingautoquantification(targetpath,Seqtarget,PNinput,option,N_slices,N_meas)
toc
fprintf('---------------------------------------------------------\n');

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loopingautoquantification(targetpath,Seqtarget,PNinput,option,N_slices,N_meas)

files       = dir(targetpath);

runitforquantification(targetpath,PNinput,Seqtarget,option,N_slices,N_meas)
fprintf('Quantification done\n');

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function runitforquantification(targetpath,PN,Seqtarget,option,N_slices,N_meas)

files       = dir(targetpath);
AIFfiles    = dir([targetpath '/AIFdata']);
DSCfiles    = dir([targetpath '/DSCanalysis']);
T1files     = dir([targetpath '/T1mapping']);
SSfiles     = dir([targetpath '/AutoSS']);
Resultfiles = dir([targetpath '/Result_MSwcf2']);

targetT1file    = [PN '_T1map.mat'];
targetSSfile    = [PN '_AutoSS.mat'];

if strmatch(Seqtarget,'GE')
    targetDSCfile{1}   = [PN 'GE_DSC_A.mat'];
    targetResultfile{1}= [PN 'GE_M.mat'];
elseif strmatch(Seqtarget,'SE')
    targetDSCfile{1}   = [PN 'SE_DSC_A.mat'];
    targetResultfile{1}= [PN 'SE_M.mat'];
elseif strmatch(Seqtarget,'all')
    targetDSCfile{1}   = [PN 'GE_DSC_A.mat'];
    targetResultfile{1}= [PN 'GE_M.mat'];
    targetDSCfile{2}   = [PN 'SE_DSC_A.mat'];
    targetResultfile{2}= [PN 'SE_M.mat'];
end

for i = 1:length(targetDSCfile)
    targetdscfile = targetDSCfile{i};
    targetAIFfile = [targetdscfile(1:6) '_AIF_A.mat'];
    targetresultfile = targetResultfile{i};
    if namefindinstrcut(files,PN,'exact') & namefindinstrcut(T1files,targetT1file,'exact') & namefindinstrcut(SSfiles,targetSSfile,'exact') & namefindinstrcut(AIFfiles,targetAIFfile,'exact') & namefindinstrcut(DSCfiles,targetdscfile,'exact')
        fprintf('%s\n',['Checking Quant for ' PN ' : ']);
        if namefindinstrcut(Resultfiles,targetresultfile,'exact')
            fprintf('%s\n',[targetpath '\' targetresultfile ' has been created...']);
            if strmatch(option,'Overwrite')
                runitforquantificationcore(targetpath,PN,targetdscfile,targetresultfile,N_slices,N_meas);
            end
        else
            fprintf('%s\n',[targetpath '\' targetresultfile ' is being created now...']);
            runitforquantificationcore(targetpath,PN,targetdscfile,targetresultfile,N_slices,N_meas);
        end
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function runitforquantificationcore(targetpath,PN,targetdscfile,targetresultfile,N_slices,N_meas)

targetSSfile =[PN '_AutoSS.mat'];

% find pathLLEPI, PWI
pathPN  = [targetpath '/' PN];
Seqfiles = dir(pathPN);
if namefindinstrcut(Seqfiles,'LL_EPI_PRE')
    targetpathLLEPI = [pathPN '/' Seqfiles(namefindinstrcut(Seqfiles,'LL_EPI_PRE')).name];
end
if strmatch(targetdscfile(5:6),'GE')
    if namefindinstrcut(Seqfiles,'ep2d_perf')
        targetpathDSC = [pathPN '/' Seqfiles(namefindinstrcut(Seqfiles,'ep2d_perf')).name];
    end
elseif strmatch(targetdscfile(5:6),'SE')
    if namefindinstrcut(Seqfiles,'ep2d_se')
        targetpathDSC = [pathPN '/' Seqfiles(namefindinstrcut(Seqfiles,'ep2d_se')).name];
    end
end

%% load data
load([targetpath '/AutoSS/' targetSSfile])  % 'ROIs', 'masks_ROIs','images','image_names'
tempROIs = ROIs;        tempimages = images;    tempimage_names = image_names;
SS = masks_ROIs.SS;
load([targetpath '/DSCanalysis/' targetdscfile])    % 'ROIs','cutoffs_ROIs','images', 'image_names'
masks_ROIs.SS = SS;

% overload
ROIs.positions  = OverloadStructdata(tempROIs.positions, ROIs.positions);
images          = OverloadCelldata(tempimages,images);
image_names     = OverloadCelldata(tempimage_names,image_names);


%% update magneticfield strength
info = dicominfo([targetpathDSC '/1.dcm'],'Dictionary','dicom-dict.txt');
if ~round((info.MagneticFieldStrength - 1.5)*10)
    ROIs.values.Tesla = '1.5T';
elseif ~round((info.MagneticFieldStrength - 3.0)*10)
    ROIs.values.Tesla = '3.0T';
end

%% make mask on LL-EPI

% YIJ 20160830: added manual or auto dog brain masking-------------------------
% disp('Draw mask for wm seg.');
% f1 = figure;imshow(images{strmatch('T1map_pre',image_names)},[]);title('Draw mask of brain for WM seg.');
% dogbrainmask = roipoly;
% close(f1);
%roi_h = impoly;
%dogbrainmask = createMask(roi_h);
%dogbrainmask = autoDogBrainMask(images{strmatch('T1map_pre',image_names)});
% -------------------------------------------------------------------------
[masks_ROIs.WM_SS,cutoffs_ROIs.WM_SS,ROIs.data.distr.T1s,ROIs.data.distr.XT1s] = fAutoMaskWM_SS(images{strmatch('T1map_pre',image_names)},ROIs.values.Tesla);

%fprintf('Please select a WM ROI:')

%Kevin Modification, enable below 3 lines to draw first time and save
if exist([targetpath '/WM_Mask_' targetresultfile], 'file') == 2
    %Kevin Modification, once drawn and use following 3 lines to load
    load([targetpath '/WM_Mask_' targetresultfile]);
    masks_ROIs.WM_SS = WM_SS;
else
    masks_ROIs.WM_SS = fManualMaskWM_SS(images{strmatch('T1map_pre',image_names)},images{strmatch('T1map_post',image_names)});
    WM_SS = masks_ROIs.WM_SS;
    save([targetpath '/WM_Mask_' targetresultfile],'WM_SS');
end


%grady 3/9/2015

% % reconsider of auto WM_SS
% if strmatch(ROIs.values.Tesla, '1.5T')
%    if cutoffs_ROIs.WM_SS(2) > 750
%        fprintf(['Warning : ' pathPN ' might not be performed properly' ]);
%    end
% elseif strmatch(ROIs.values.Tesla, '3.0T')
%    if cutoffs_ROIs.WM_SS(2) > 100
%        fprintf(['Warning : ' pathPN ' might not be performed properly' ]);
%    end
% end

%% make mask on PWI
% [DSC,headers,Nslices] = ReadDICOMSeries(targetpathDSC,1);
Nslices = N_slices; %Jessy (10/31/2008): give slice number manually (unable to read header info)
[mask, nslice, match_index] = fAutoMaskWM_DSC_Philips(targetpathLLEPI,targetpathDSC,masks_ROIs.WM_SS,Nslices,N_meas);

%save data
ROIs.positions.n_slice_WM_DSC = nslice;
masks_ROIs.WM_DSC = mask;

% find cutoffs in cbvdsc distribution
tempCBV = images{strmatch('CBV_DSC',image_names)}(:,:,nslice);

tempCBVwm = tempCBV(find(masks_ROIs.WM_DSC));
%tempCBVwm = tempCBVwm(tempCBVwm > 0);

maxcbv = max(tempCBVwm);
%[ROIs.data.distr.CBVdsc ROIs.data.distr.XCBVdsc] = hist(tempCBV(find(masks_ROIs.WM_DSC)),0:maxcbv/50:maxcbv);
%grady 03.09.15
mincbv = min(tempCBVwm);
%figure;
%[ROIs.data.distr.CBVdsc ROIs.data.distr.XCBVdsc] = hist(tempCBVwm,mincbv:(maxcbv-mincbv)/20:maxcbv);
%title('CBV_DSC hist');

% YIJ 20200213: only count CBVdsc > 1, because less than 1 (small numbers like in csf) can bias
[ROIs.data.distr.CBVdsc edges] = histcounts(tempCBVwm(tempCBVwm > 1),'binmethod','fd');
%[ROIs.data.distr.CBVdsc edges] = histcounts(tempCBVwm,'binmethod','fd');
width = edges(2) - edges(1);
ROIs.data.distr.XCBVdsc = edges(1:end-1) + width/2;

cutoffs_ROIs.CBVdsc = zeros(1,2);
temp = find(ROIs.data.distr.CBVdsc);
cutoffs_ROIs.CBVdsc(1) = ROIs.data.distr.XCBVdsc(min(temp));
cutoffs_ROIs.CBVdsc(2) = ROIs.data.distr.XCBVdsc(max(temp));

% fit CBVdsc distribution
ROIs.data.fitting.CBVdsc = CalfitCBVdsc(ROIs.data.distr.XCBVdsc,ROIs.data.distr.CBVdsc,cutoffs_ROIs.CBVdsc,PN);
CBVdsc_WM = ROIs.data.fitting.CBVdsc.b1;

% create CBVss map
global injectionNum
switch injectionNum
    case 1
        [CBVssFast_per,CBVssSlow_per,T1s] = CalcCBVSSmap(images,image_names,ROIs.positions,masks_ROIs.WM_SS);
    case 2
        [CBVssFast_per,CBVssSlow_per,T1s] = CalcCBVSSmap(images,image_names,ROIs.positions,masks_ROIs.WM_SS);
    otherwise
        [CBVssFast_per,CBVssSlow_per,T1s] = CalcCBVSSmap(images,image_names,ROIs.positions,masks_ROIs.WM_SS);
        warning(['Injection number is greater than 2: ' num2str(injectionNum)]);
end

ROIs.values.T1s = T1s;

CBVssFastwm = CBVssFast_per(find(masks_ROIs.WM_SS));
%CBVssFastwm = CBVssFastwm(CBVssFastwm > 0);

% create CBVss map and calculate cutoffs
maxcbvss = max(CBVssFastwm);
mincbvss = min(CBVssFastwm);

%[ROIs.data.distr.CBVssFast ROIs.data.distr.XCBVssFast] = hist(CBVssFast_per(find(masks_ROIs.WM_SS))',-10:0.2:10);
%figure;
%[ROIs.data.distr.CBVssFast ROIs.data.distr.XCBVssFast] = hist(CBVssFast_per(find(masks_ROIs.WM_SS))',-40:0.2:40);
%YIJ -40:0.2:40 to mincbvss:maxcbvss/50:maxcbvss
%[ROIs.data.distr.CBVssFast ROIs.data.distr.XCBVssFast] = hist(CBVssFastwm,mincbvss:(maxcbvss-mincbvss)/20:maxcbvss);
%title('CBVssFast hist');

% YIJ 20200213: only count CBVss > .1 to try remove negative voxels
% (*20200217: doesn't really work)
[ROIs.data.distr.CBVssFast edges] = histcounts(CBVssFastwm,'binmethod','fd');
width = edges(2) - edges(1);
ROIs.data.distr.XCBVssFast = edges(1:end-1) + width/2;

cutoffs_ROIs.CBVssFast = zeros(1,2);
temp = find(ROIs.data.distr.CBVssFast);
cutoffs_ROIs.CBVssFast(1) = ROIs.data.distr.XCBVssFast(min(temp));
cutoffs_ROIs.CBVssFast(2) = ROIs.data.distr.XCBVssFast(max(temp));

% fit CBVss distribution
ROIs.data.fitting.CBVssFast = CalfitCBVss(ROIs.data.distr.XCBVssFast,ROIs.data.distr.CBVssFast,cutoffs_ROIs.CBVssFast, PN);
CBVss_WM = ROIs.data.fitting.CBVssFast.b1;

%YIJ 20170726
% tempConc = images{strmatch('Conc_map',image_names)}(:,:,nslice);
% maxconc = max(tempConc(find(masks_ROIs.WM_DSC)));
% minconc = min(tempConc(find(masks_ROIs.WM_DSC)));
% [ROIs.data.distr.Conc ROIs.data.distr.XConc] = hist(tempConc(find(masks_ROIs.WM_DSC)),minconc:(maxconc-minconc)/20:maxconc);
% cutoffs_ROIs.Conc = zeros(1,2);
% temp = find(ROIs.data.distr.Conc);
% cutoffs_ROIs.Conc(1) = ROIs.data.distr.XConc(min(temp));
% cutoffs_ROIs.Conc(2) = ROIs.data.distr.XConc(max(temp));
% ROIs.data.fitting.Conc = CalfitConc(ROIs.data.distr.XConc,ROIs.data.distr.Conc,cutoffs_ROIs.Conc,PN);


% quantification
[ROIs.values.CF,qCBF]   = Calc_QCBV_CBF_MS2(images{strmatch('CBF_nSVD',image_names)},ROIs.values.T1s,ROIs.data.fitting,ROIs.values.Tesla);
[ROIs.values.CF,qCBV]   = Calc_QCBV_CBF_MS2(images{strmatch('CBV_DSC',image_names)},ROIs.values.T1s,ROIs.data.fitting,ROIs.values.Tesla);


[images, image_names] = IA_ImportImages(images, image_names,'qCBF_nSVD',qCBF,'qCBV_DSC',qCBV);

[images, image_names] = IA_ImportImages(images, image_names,'CBVssFast',CBVssFast_per);

% save the final result
if strmatch (match_index,'unmatch')
    file_length = length(targetresultfile)-5;
    targetresultfile = [targetresultfile(1:file_length) 'PA.mat'];
end
save([targetpath '/Result_MSwcf2/' targetresultfile],'ROIs','cutoffs_ROIs','images','image_names','masks_ROIs');

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = OverloadStructdata(input_adding, input_added)

output = input_added;
subfields = fieldnames(input_adding);
for i = 1:length(subfields)
    eval(['output.' subfields{i} ' = input_adding.' subfields{i} ';'])
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = OverloadCelldata(input_adding, input_added)

output = input_added;
count = length(input_added);
for i = 1:length(input_adding)
    output{count+i} = input_adding{i};
end

return

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%% DICOM conversion
function autoWriteDICOM(varargin)

targetpath  = varargin{1};
PNinput     = varargin{2};
Seqtarget   = varargin{3};
option      = 'Skip';

loopingautoWriteDICOM(targetpath,Seqtarget,PNinput,option)

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loopingautoWriteDICOM(targetpath,Seqtarget,PNinput,option)

files       = dir(targetpath);

runitforWriteDICOM(targetpath,PNinput,Seqtarget,option)
fprintf('done\n');

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function runitforWriteDICOM(targetpath,PN,Seqtarget,option)

pathPN = [targetpath '/' PN];
Seqfiles = dir(pathPN);
Resultfiles = dir([targetpath '/Result']);

metrics = {'qCBV','qCBF','MTT'};

for i = 1:length(metrics)
    metric = metrics{i};
    DICOMdirectories = dir([targetpath '/DICOM/' metric 'maps']);
    
    if strmatch(Seqtarget,'GE')
        targetresultfile{1} = [PN 'GE_A.mat'];
        targetdicomDir{1} = [PN 'GE_' metric];
        if namefindinstrcut(Seqfiles,'ep2d_perf')
            hdrfilepath{1} = [pathPN '/' Seqfiles(namefindinstrcut(Seqfiles,'ep2d_perf')).name '/1.dcm'];
        end
    elseif strmatch(Seqtarget,'SE')
        targetresultfile{1} = [PN 'SE_A.mat'];
        targetdicomDir{1} = [PN 'SE_' metric];
        if namefindinstrcut(Seqfiles,'ep2d_se')
            hdrfilepath{1} = [pathPN '/' Seqfiles(namefindinstrcut(Seqfiles,'ep2d_se')).name '/1.dcm'];
        end
    elseif strmatch(Seqtarget,'all')
        targetresultfile{1} = [PN 'GE_A.mat'];
        targetdicomDir{1} = [PN 'GE_' metric];
        if namefindinstrcut(Seqfiles,'ep2d_perf')
            hdrfilepath{1} = [pathPN '/' Seqfiles(namefindinstrcut(Seqfiles,'ep2d_perf')).name '/1.dcm'];
        end
        targetresultfile{2} = [PN 'SE_A.mat'];
        targetdicomDir{2} = [PN 'SE_' metric];
        if namefindinstrcut(Seqfiles,'ep2d_se')
            hdrfilepath{2} = [pathPN '/' Seqfiles(namefindinstrcut(Seqfiles,'ep2d_se')).name '/1.dcm'];
        end
    end
    
    for i = 1:length(targetresultfile)
        targetRESULTfile = targetresultfile{i};
        targetDICOMDir = targetdicomDir{i};
        headerfilepath = hdrfilepath{i};
        if strfind(targetRESULTfile, 'GE')
            Sequence = 'GE';
        elseif strfind(targetRESULTfile, 'SE')
            Sequence = 'SE';
        end
        if namefindinstrcut(Resultfiles,targetRESULTfile,'exact')
            fprintf('%s\n',['Checking ' PN ' : ']);
            if namefindinstrcut(DICOMdirectories,targetDICOMDir,'exact')
                fprintf('%s\n',[targetpath '/DICOM/' metric 'maps/' targetDICOMDir ' has been created...']);
                if strmatch(option,'Overwrite')
                    runitforWriteDICOMcore(targetpath,PN,targetRESULTfile,targetDICOMDir,metric,headerfilepath,Sequence);
                end
            else
                fprintf('%s\n', [targetpath '/DICOM/' metric 'maps/' targetDICOMDir ' is being created now...']);
                runitforWriteDICOMcore(targetpath,PN,targetRESULTfile,targetDICOMDir,metric,headerfilepath,Sequence);
            end
        end
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function runitforWriteDICOMcore(targetpath,PN,targetRESULTfile,targetDICOMDir,metric,headerfilepath,Sequence)

load([targetpath '/Result/' targetRESULTfile]);

switch Sequence
    case 'GE'
        seriesNumber = 90;
        factor = 1000;
    case 'SE'
        seriesNumber = 95;
        factor = 5000;
    otherwise
        seriesNumber = 99;
        factor = 9000;
end

switch metric
    case 'qCBV'
        image = images{15};
        SF    = 100; %scaling factor
    case 'qCBF'
        image = images{14};
        SF    = 10;
        seriesNumber = seriesNumber + 1;
        factor = factor + 1000;
    case 'MTT'
        image = images{3};
        SF    = 100;
        seriesNumber = seriesNumber + 2;
        factor = factor + 2000;
    otherwise
        fprintf('%s\n', ['Please specify the image set corresponding to: ' metric ]);
end

%move to right directory
cd([targetpath '/DICOM/' metric 'maps'])

writeDICOM_JJM(image,headerfilepath,seriesNumber,targetDICOMDir,factor,SF);

cd(targetpath)

return

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%% Other Functions
function PN = makePN(i)

if i < 10
    PN = ['P00' num2str(i)];
elseif i<100
    PN = ['P0' num2str(i)];
else
    PN = ['P' num2str(i)];
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function index = namefindinstrcut(inputstruct, string, opt)

index=0;
if nargin < 3
    for r = 1:length(inputstruct)
        if strfind(inputstruct(r).name,string)
            index=r;
            break
        end
    end
elseif strmatch(opt,'exact')
    for r = 1:length(inputstruct)
        if strmatch(inputstruct(r).name,string,'exact')
            index=r;
            break
        end
    end
end

return