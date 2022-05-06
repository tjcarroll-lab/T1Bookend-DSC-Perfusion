function [] = debugDDaif(images,AIF,Dt,Wij,echoTime,DDfit,kATDmap,kIBATmap,kBATmap,kBRTmap,k,Nslices)
%debugDDaif look at aif before and after DD, and regional conc. curve
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2018-03-06

global glblTargetPath;
global injectionNum;
global seqType;
tmpdir = dir([glblTargetPath '\ROIs\*ROI_20180125v2.mat']);
if ~isempty(tmpdir)
    roidata = load([glblTargetPath '\ROIs\' tmpdir(1).name]);
    roistack = roidata.roi_stack;
    
    if size(roistack) ~= Nslices
        error('roi mask and dsc Nslices do not match');
    end
    tmproi = roistack(:,:,k);
    roinums = unique(tmproi);
    roinums = roinums(roinums > 0);
    
    oAIF = AIF;
    signals = {};
    for roiind = 1:length(roinums)
        avgsignal = zeros(1,size(images,3));
        count = 0;
        avgibat = 0;
        avgbat = 0;
        avgbrt = 0;
        avgatd = 0;
        for i = 1:size(images,1)
            for j = 1:size(images,2)
                if tmproi(i,j) == roinums(roiind)
                    tmpsignal = squeeze(images(i,j,:))';
                    avgsignal = avgsignal + reshape(smooth(tmpsignal,5,'sgolay',3),[1 size(images,3)]);
                    count = count + 1;
                    avgibat = avgibat + kIBATmap(i,j);
                    avgbat = avgbat + kBATmap(i,j);
                    avgbrt = avgbrt + kBRTmap(i,j);
                    avgatd = avgatd + kATDmap(i,j);
                end
            end
        end
        avgsignal = avgsignal ./ count;
        avgibat = round(avgibat ./ count);
        avgbat = round(avgbat ./ count);
        avgbrt = round(avgbrt ./ count);
        avgatd = round(avgatd ./ count);
        
        signals{roiind} = avgsignal;        
        
        % get regional conc. curve
        Cmt = Sig2Conct_AIF(avgsignal, [avgibat avgbat avgbrt], echoTime);
        % apply DDfix to AIF
        AIF = applyDDfix(oAIF,Dt,avgatd,DDfit,Wij);
        
        figure;plot(oAIF.Ct);hold on;plot(AIF.Ct);plot(Cmt.Ct(avgbat:avgbrt));
        title(['k: ' num2str(k) ' ROI: ' num2str(roinums(roiind)) ' ATD: ' num2str(avgatd) ' ' num2str(avgibat) '-' num2str(avgbat) '-' num2str(avgbrt)]);
    end
end

end

