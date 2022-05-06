function [ isbolus,scndfilt ] = boluscheck( Signal,oSignal,oAIF,bats,echoTime )
%boluscheck check if bolus or too much noise
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2019-07-23
% ibat = bats(1);
% bat = bats(2);
% brt = bats(3);

%---------------------------------------------
% YIJ 20200113: signal is smoothed before this function
if isempty(Signal)
%tmpfitopt = fitoptions('method','smoothingspline','smoothingparam',.1);
%tmpfittype = fittype('smoothingspline');
%tmpfit = fit([1:length(oSignal)]',oSignal',tmpfittype,tmpfitopt);
%Signal = reshape(tmpfit(1:length(oSignal)),[1 length(oSignal)]);
%Signal(Signal < 0) = 0;
% 20210422 YIJ: remove smoothing because of arrival time shift
Signal = oSignal;
end
            
%Signal = smooth(oSignal,.05,'sgolay',1)';
%Signal = smooth(oSignal,.1,'sgolay',1)';
%Signal(Signal < 0) = 0;

if isempty(bats)
    [ibat,bat,brt] = findBAT_test(Signal);
else
    ibat = bats(1);
    bat = bats(2);
    brt = bats(3);
end

Cmt = Sig2Conct_AIF(Signal, [ibat bat brt], echoTime);
oCmt = Cmt;

% 20190729 YIJ:
% If concentration after brt remains very high compared to baseline,
% subtract by its average to set mean to zero (not sure what causes this)
if mean(Cmt.Ct(ibat:bat)) < mean(Cmt.Ct(end-20:end)) - 1.5*std(oCmt.Ct)
    Cmt.Ct(brt+1:end) = Cmt.Ct(brt+1:end) - mean(Cmt.Ct(brt+1:end));
end

% 20190729 YIJ:
% Omit initial few data points to remove noise, non-steadystate signal
% affecting peak analysis
Cmt.Ct = Cmt.Ct(ibat:end);
bat = bat-(ibat-1);
brt = brt-(ibat-1);

% 20191210 YIJ: SAVE ------------------------------------------------------
% % try
% [height,ind,width,prom] = findpeaks(Cmt.Ct,'minpeakdistance',round(.05*length(Signal)));
% % catch
% %     1;
% % end
% [~,maxheightind] = max(height);
% [~,minheightind] = min(height);
% maxprom = prom(maxheightind);
% if length(prom) > 2
% prom([maxheightind minheightind]) = [];
% end
% -------------------------------------------------------------------------

[height,ind,width,prom] = findpeaks(Cmt.Ct,'widthreference','halfheight');

if length(prom) > 1
lowprom = prom < 0.1*mean(prom);
height(lowprom) = [];
ind(lowprom) = [];
width(lowprom) = [];
prom(lowprom) = [];
end

% YIJ 20200123: bolus should be highest peak so check prominence of highest
% peak
[~,mhi] = max(height);
maxprom = prom(mhi);
%prom(prom < mean(prom)) = [];

%height = height(prom > mean(prom));
%ind = ind(prom > mean(prom));
[maxconc,maxind] = max(Cmt.Ct);
%if (maxconc > mean(height) + 2*std(height)) && maxind < brt && maxind > bat
scndfilt = 0;
[~,mpi] = max(prom);
[~,mwi] = max(width);
[~,mpwi] = max(prom.*width);
[~,mhwi] = max(height.*width);

if mhi ~= mpi
    1;
end

tmpprom = prom;
tmpind = ind;
tmpprom(mhi) = [];
tmpind(mhi) = [];

outsideBolusPeaks = sum( (height > 0.5*maxprom) & (ind > brt | ind < bat) );

if ~isempty(prom)
    
    criteria = [( (maxprom > mean(prom) + 1.96*std(oCmt.Ct)) || (maxprom > mean(prom) + 1.96*std(oCmt.Ct(ibat:bat))) ) ... %1
                (ind(mhi) < brt && ind(mhi) > bat) ... %2
                mpi == mpwi ... %3
                mpi == mhi ... %4
                mpi == mhwi ... %5
                sum(ind < brt & ind > bat) < 4 ... %6
                maxprom > sum(tmpprom( (tmpind < brt & tmpind > bat) )) ]; %7
        
    if sum(criteria) == length(criteria)
        isbolus = 1;
    elseif criteria(5) == 0 && sum(criteria) == length(criteria)-1
        if outsideBolusPeaks == 0
            isbolus = 1;
            scndfilt = 1;
        else
            isbolus = 0;
        end
    elseif length(prom) == 1 && ind(mhi) < brt && ind(mhi) > bat
        isbolus = 1;
    else % Secondary filter
        isbolus = 0;
        %prom(prom < mean(prom)) = [];
        test = sort(prom,'descend');       
        if length(prom) > 1
            if 0.7*maxprom > test(2) && ...
                    ind(mpi) < brt && ind(mpi) > bat && ...'
                    mpi == mpwi && ...
                    sum(ind < brt & ind > bat) < 3 
                isbolus = 0;
                %scndfilt = 1;
            else
                isbolus = 0;
            end
        end
    end
    
%     %if (maxprom > mean(prom) + 2*std(prom)) &&...
%     if (maxprom > mean(prom) + 1.96*std(oCmt.Ct)) &&...
%             ind(mhi) < brt && ind(mhi) > bat &&...
%             mpi == mpwi && mpi == mhi && mpi == mhwi &&...
%             sum(ind < brt & ind > bat) < 4 &&...
%             maxprom > sum(tmpprom( (tmpind < brt & tmpind > bat) ))
%         isbolus = 1;
%     elseif length(prom) == 1 && ind(mhi) < brt && ind(mhi) > bat
%         isbolus = 1;
%     elseif (maxprom > mean(prom) + 1.96*std(oCmt.Ct(ibat:bat))) && ...
%             ind(mhi) < brt && ind(mhi) > bat && ...
%             mpi == mpwi && mpi == mhi && mpi == mhwi &&...
%             sum(ind < brt & ind > bat) < 4 &&...
%             maxprom > sum(tmpprom( (tmpind < brt & tmpind > bat) ))
%         isbolus = 1;
%     else % Secondary filter
%         isbolus = 0;
%         %prom(prom < mean(prom)) = [];
%         test = sort(prom,'descend');       
%         if length(prom) > 1
%             if 0.7*maxprom > test(2) && ...
%                     ind(mpi) < brt && ind(mpi) > bat && ...'
%                     mpi == mpwi && ...
%                     sum(ind < brt & ind > bat) < 3 
%                 isbolus = 0;
%                 scndfilt = 1;
%             else
%                 isbolus = 0;
%             end
%         end
%     end   
    
else
    isbolus = 0;
end

