function [infarcttbl, normaltbl] = getDataTbl(output, cortexroi, infarctData, PCS, vols2hr)
%getDataTbl Returns infarct and normal side roi values in table format
%   Finds the common cases between perf, roi and infarct data. Then
%   calculates the average values for CBF, CBV, MTT, percent CBF diff, CBF
%   diff, Tmax, volume of increased CBF voxels.
%   Inputs:
%    - output: contains perf DD and noDD data
%    - cortexroi: contains cortical ROIs (infarct side and normal side)
%    - infarctData: contains infarct volume data
%    - PCS: pial collateral scores
%    - vols2hr: infarct volume at approx. 2 hrs
%
% Author: Yong Ik Jeong
% Date: 2020-11-02

allcbf = [];
allddcbf = [];
alltmax = [];
allddtmax = [];

% Get common case ids
[caseids,ia,ib] = intersect({output.id},{cortexroi.id});
caseids = cellstr(string(regexpi(caseids, '\d\d', 'match')));
if istable(PCS)
    [caseids,ia,ib] = intersect(caseids, PCS.Case);
else
    pcs = PCS;
    PCS = table(caseids', pcs', 'variablenames',{'Case','PCS'});
end
if ~isempty(infarctData)
    [caseids,ia,ib] = intersect(caseids, {infarctData.case});
end

% Organize data using finalized case ids
[c,ia,ib] = intersect(cellstr(string(regexpi({output.id}, '\d\d', 'match'))), caseids);
output = output(ia);
[c,ia,ib] = intersect(cellstr(string(regexpi({cortexroi.id}, '\d\d', 'match'))), caseids);
cortexroi = cortexroi(ia);
[c,ia,ib] = intersect(PCS.Case, caseids);
PCS = table2array(PCS(ia,'PCS'));

if ~isempty(infarctData)
    [c,ia,ib] = intersect({infarctData.case}, caseids);
    infarctData = infarctData(ia);
end

if isempty(vols2hr)
    vols2hr = zeros(length(output),1);
else
    vols2hr = reshape(vols2hr,[length(vols2hr) 1]);
end

infarcttbl = struct();
infarcttbl.infarctvol = vols2hr;
infarcttbl.pcs = PCS;
normaltbl = struct();
normaltbl.infarctvol = vols2hr;
normaltbl.pcs = PCS;
tmaxvols = [];
ddtmaxvols = [];
volerr = [];
TRthresh = [];

for ii = 1:length(output)
    t1_slc_ind = output(ii).result.ROIs.positions.n_slice_WM_DSC;
    TR = output(ii).header.RepetitionTime;
    tmax = output(ii).result.images.noDD.Tmax_map(:,:,t1_slc_ind);
    isbolus = output(ii).result.images.noDD.isbolus(:,:,t1_slc_ind);
    atdmap = output(ii).result.images.etc.ATDmap(:,:,t1_slc_ind);    
    ddcbf = output(ii).result.images.DD.qCBF_SVD(:,:,t1_slc_ind);
    cbf = output(ii).result.images.noDD.qCBF_SVD(:,:,t1_slc_ind);
    cbv = output(ii).result.images.DD.qCBV_DSC(:,:,t1_slc_ind);
    mtt = output(ii).result.images.DD.CMTT_CVP(:,:,t1_slc_ind);
    volunit = output(ii).header.PixelSpacing(1) * ...
        output(ii).header.PixelSpacing(2) * ...
        output(ii).header.SliceThickness;
    
    % Infarct--------------------------------------------------------------
    roi = logical(cortexroi(ii).infarct_roi);
    tmaxhole = sum(roi == 1 & isbolus == 0,'all') * volunit;
    %roi(isbolus == 0) = 0;
    roi(ddcbf < 1) = 0;
    roi(cbf < 1) = 0;
          
    avgcbf = mean(cbf(roi));
    avgddcbf = mean(ddcbf(roi));
    avgdiff = mean(ddcbf(roi) - cbf(roi));   
    avgdiffvol = sum((ddcbf(roi) - cbf(roi))./cbf(roi) > .10,'all') .* volunit;
    avgperdiff = mean((ddcbf(roi) - cbf(roi))./cbf(roi)) .* 100;
    avgatd = mean(atdmap(roi)).*(TR/1000);
    avgtmax = mean(tmax(roi & isbolus == 1));
    avgcbv = mean(cbv(roi));
    avgmtt = mean(mtt(roi));
    %mtt = 60*cbv./ddcbf;
    %roi(ddcbf < 1 | cbv < .1) = 0;
    %roi(ddcbf < .1) = 0;
    %avgmtt = mean(mtt(roi));
    
    perchange = (ddcbf-cbf)./cbf;
    %yv(perchange,[0 1],'ol',perchange > .1 & roi == 1)
    yv(ddcbf,[0 250],'ol',cortexroi(ii).infarct_roi == 1);
    
    infarcttbl.cbf(ii,1) = avgcbf;
    infarcttbl.ddcbf(ii,1) = avgddcbf;
    infarcttbl.tmax(ii,1) = avgtmax;
    infarcttbl.cbv(ii,1) = avgcbv;
    infarcttbl.mtt(ii,1) = avgmtt;
    infarcttbl.tmaxhole(ii,1) = tmaxhole;
    infarcttbl.diff(ii,1) = avgdiff;
    infarcttbl.diffvol(ii,1) = avgdiffvol;
    infarcttbl.perdiff(ii,1) = avgperdiff;
    infarcttbl.atd(ii,1) = avgatd;
    infarcttbl.caseid{ii,1} = output(ii).id;
     
    % Normal---------------------------------------------------------------
    roi = logical(cortexroi(ii).normal_roi);
    roi(ddcbf < 1) = 0;
    roi(cbf < 1) = 0;
    %roi(isbolus == 0) = 0; 
    
    avgcbf = mean(cbf(roi));
    avgddcbf = mean(ddcbf(roi));
    avgdiff = mean(ddcbf(roi) - cbf(roi));   
    avgdiffvol = sum((ddcbf(roi) - cbf(roi))./cbf(roi) > .10,'all') .* volunit;
    avgperdiff = mean((ddcbf(roi) - cbf(roi))./cbf(roi)) .* 100;
    avgatd = mean(atdmap(roi)).*(TR/1000);
    avgtmax = mean(tmax(roi & isbolus == 1));
    avgcbv = mean(cbv(roi));
    avgmtt = mean(mtt(roi));
    %mtt = 60*cbv./ddcbf;
    %roi(ddcbf < 1 | cbv < .1) = 0;
    %roi(ddcbf < .1) = 0;
    %avgmtt = mean(mtt(roi));
    
    normaltbl.cbf(ii,1) = avgcbf;
    normaltbl.ddcbf(ii,1) = avgddcbf;
    normaltbl.tmax(ii,1) = avgtmax;
    normaltbl.cbv(ii,1) = avgcbv;
    normaltbl.mtt(ii,1) = avgmtt;
    normaltbl.diff(ii,1) = avgdiff;
    normaltbl.diffvol(ii,1) = avgdiffvol;
    normaltbl.perdiff(ii,1) = avgperdiff;
    normaltbl.atd(ii,1) = avgatd;
    normaltbl.caseid{ii,1} = output(ii).id;
        
end

infarcttbl = struct2table(infarcttbl);
normaltbl = struct2table(normaltbl);
fprintf('Done\n');

end

