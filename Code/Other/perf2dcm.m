function [] = perf2dcm(bp,PNin,savepath,savenamepattern,quant)
%perf2dcm Save perfusion data into dicoms
%   Saves qCBF, qCBV, and MTT
%   Inputs:
%       - bp: ex. '..\RESOURCES\DICOMS\COLLAT\COLLAT_31\d02_20200108'
%       - PNin: ex. {'P001','P002'}
%       - savepath: ex. '..\SAVES\COLLAT\COLLAT_31\d02'
%       - savenamepattern: ex. 'COLLAT_31_d02'
%   Outputs:
%       - dicoms will saved to [savepath '\P001\qCBF\' savenamepattern
%       '_P001_qCBF_01.dcm']
%
% Author: Yong Ik Jeong
% Date: 2020-07-13
% Changelog:
%   - 20200713 YIJ: Initial version
%   - 20200717 YIJ: Added compatibility to old and new results files and
%   some organization

targetpath = bp;
input = dir([targetpath '/P*']);
for ii = 1:length(input)
    PN = input(ii).name;
    slclist = [];
    
    % Run only specified PN (ex. P001)
    if sum(strcmpi(PN,PNin))
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
        
        % Get number of timepoints and slices
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
        
        % Load the Results mat file
        %resultdir = dir(fullfile(targetpath,'Result_MSwcf2',[PN '*.mat']));
        resultdir = dir(fullfile(targetpath,'DSCanalysis',[PN '*.mat']));
        if length(resultdir) == 1
            tmpmat = matfile(fullfile(resultdir(1).folder,resultdir(1).name));
            images = tmpmat.images;
            tmpvars = who(tmpmat);
            if sum(strcmpi(tmpvars,'image_names'))
                oldmethod = 1;
                image_names = tmpmat.image_names;
            else
                oldmethod = 0;
                ROIs = tmpmat.ROIs;
            end            
        elseif isempty(resultdir)
            continue;
        else
            error('More than one results file');
        end
        
        if quant
            mkdir(fullfile(savepath,PN,'qCBF'));
            mkdir(fullfile(savepath,PN,'qCBV'));
        end
        mkdir(fullfile(savepath,PN,'rCBF'));       
        mkdir(fullfile(savepath,PN,'rCBV'));
        mkdir(fullfile(savepath,PN,'MTT'));
        mkdir(fullfile(savepath,PN,'Tmax'));
        mkdir(fullfile(savepath,PN,'TTP'));
               
        % Get dicom header info from first image of every DSC slice and
        % save perfusion data as dicoms using that header
        count = 1;
        for ind = 1:N_meas:N_meas*N_slices
            tmphdr = dicominfo(fullfile(perfpath,[num2str(ind) '.dcm']));
            if oldmethod
                if quant
                    tmpqcbf = squeeze(images{strcmpi(image_names,'qCBF_nSVD')}(:,:,count));
                    tmpqcbv = squeeze(images{strcmpi(image_names,'qCBV_DSC')}(:,:,count));
                end
                tmprcbf = squeeze(images{strcmpi(image_names,'CBF_nSVD')}(:,:,count));                
                tmprcbv = squeeze(images{strcmpi(image_names,'CBV_DSC')}(:,:,count));
                tmpmtt = squeeze(images{strcmpi(image_names,'CMTT_nSVD')}(:,:,count));
                tmptmax = squeeze(images{strcmpi(image_names,'Tmax_map')}(:,:,count));
                tmpttp = squeeze(images{strcmpi(image_names,'TTP_map')}(:,:,count));
            else
                if ROIs.values.DDfix == 1
                    if quant
                        tmpqcbf = squeeze(images.DD.qCBF_SVD(:,:,count));
                        tmpqcbv = squeeze(images.DD.qCBV_DSC(:,:,count));
                    end
                    tmprcbf = squeeze(images.DD.CBF_SVD(:,:,count));                    
                    tmprcbv = squeeze(images.DD.CBV_DSC(:,:,count));
                    tmpmtt = squeeze(images.DD.CMTT_CVP(:,:,count));
                else
                    if quant
                        tmpqcbf = squeeze(images.noDD.qCBF_SVD(:,:,count));
                        tmpqcbv = squeeze(images.noDD.qCBV_DSC(:,:,count));
                    end
                    tmprcbf = squeeze(images.noDD.CBF_SVD(:,:,count));                   
                    tmprcbv = squeeze(images.noDD.CBV_DSC(:,:,count));
                    tmpmtt = squeeze(images.noDD.CMTT_CVP(:,:,count));
                end
                tmptmax = squeeze(images.noDD.Tmax_map(:,:,count));
                tmpttp = squeeze(images.noDD.TTP_map(:,:,count));
            end
            % Scale because integer
            if quant
                tmpqcbf = uint16(tmpqcbf*10);
                tmpqcbv = uint16(tmpqcbv*100);
            end
            tmprcbf = uint16(tmprcbf*10);
            tmprcbv = uint16(tmprcbv*100);
            tmpmtt(tmpmtt > 100) = 0;
            tmpmtt = uint16(tmpmtt*100);
            tmptmax = uint16(tmptmax*100);
            tmpttp = uint16(tmpttp*100);
            
            if quant
                filename = fullfile(savepath,PN,'qCBF',[savenamepattern '_' PN '_qCBF_' sprintf('%02d',count) '.dcm']);
                dicomwrite(real(tmpqcbf),filename,tmphdr);
                fprintf('Saved to %s\n',filename);
                
                filename = fullfile(savepath,PN,'qCBV',[savenamepattern '_' PN '_qCBV_' sprintf('%02d',count) '.dcm']);
                dicomwrite(real(tmpqcbv),filename,tmphdr);
                fprintf('Saved to %s\n',filename);
            end
            
            filename = fullfile(savepath,PN,'rCBF',[savenamepattern '_' PN '_rCBF_' sprintf('%02d',count) '.dcm']);           
            dicomwrite(real(tmprcbf),filename,tmphdr);
            fprintf('Saved to %s\n',filename);
                        
            filename = fullfile(savepath,PN,'rCBV',[savenamepattern '_' PN '_rCBV_' sprintf('%02d',count) '.dcm']);           
            dicomwrite(real(tmprcbv),filename,tmphdr);
            fprintf('Saved to %s\n',filename);
            
            filename = fullfile(savepath,PN,'MTT',[savenamepattern '_' PN '_MTT_' sprintf('%02d',count) '.dcm']);           
            dicomwrite(real(tmpmtt),filename,tmphdr);
            fprintf('Saved to %s\n',filename);
            
            filename = fullfile(savepath,PN,'Tmax',[savenamepattern '_' PN '_Tmax_' sprintf('%02d',count) '.dcm']);           
            dicomwrite(real(tmptmax),filename,tmphdr);
            fprintf('Saved to %s\n',filename);
            
            filename = fullfile(savepath,PN,'TTP',[savenamepattern '_' PN '_TTP_' sprintf('%02d',count) '.dcm']);           
            dicomwrite(real(tmpttp),filename,tmphdr);
            fprintf('Saved to %s\n',filename);
            count = count + 1;
        end
    end   
end
fprintf('Done\n');