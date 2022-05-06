function [ vout ] = getAllPerf( perfmr, perfmcr, regions, co2data )
%getAllPerf Combine MRI and microsphere cbf values regional and hemispheric
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2017-07-30

casenum_mr = [perfmr.casenum];
casenum_mcr = [perfmcr.casenum];
casenum_regions = [regions.casenum];
casenum_co2 = [co2data.casenum];

regio_mr = [];
regio_mcr = [];
hemis_mr = [];
hemis_mcr = [];
allregionum = [];
allhemislbl = [];
allco2 = [];
allregiolbl = [];

count = 0;
count2 = 0;
case_sep = [];
case_sep2 = [];

for ii = 1:length(perfmr)
    casenum = perfmr(ii).casenum;
    mcrind = ismember(casenum_mcr,casenum);
    regioind = ismember(casenum_regions,casenum);
    co2ind = ismember(casenum_co2,casenum);
    if sum(mcrind) && sum(regioind) && sum(co2ind)
        regionum_mr = [perfmr(ii).cbfpre{:,1}];
        regionum_mcr = [perfmcr(mcrind).cbfpre{:,1}];
        regionum = [regions(regioind).regionlist{:,5}];
        
        baseco2 = co2data(co2ind).base;
        stressco2 = co2data(co2ind).stress;
        
        %allco2 = [allco2; [baseco2 stressco2]];
        
        [regionum_inter,tmpindmr,tmpindmcr] = intersect(regionum_mr,regionum_mcr); %returns sorted values, indices
        %tmpind = ismember(regionum_mr,regionum_mcr);
        tmppremr = [perfmr(ii).cbfpre{tmpindmr,2}]';
        tmppostmr = [perfmr(ii).cbfpost{tmpindmr,2}]';
        %tmpcvrmr = (tmppostmr-tmppremr)./tmppremr;
        %tmpcvrmr = (tmppostmr-tmppremr)./(stressco2-baseco2);
        tmpcvrmr = (tmppostmr-tmppremr)./tmppremr./(stressco2-baseco2) .* 100;
        tmpmat = [tmppremr tmppostmr tmpcvrmr];
                  
        regio_mr = [regio_mr; tmpmat]; 
        
        %tmpind = ismember(regionum_mcr,regionum_mr);
        tmppremcr = [perfmcr(mcrind).cbfpre{tmpindmcr,2}]'.*100;
        tmppostmcr = [perfmcr(mcrind).cbfpost{tmpindmcr,2}]'.*100;
        %tmpcvrmcr = (tmppostmcr-tmppremcr)./tmppremcr;
        %tmpcvrmcr = (tmppostmcr-tmppremcr)./(stressco2-baseco2);
        tmpcvrmcr = (tmppostmcr-tmppremcr)./tmppremcr./(stressco2-baseco2) .* 100;
        tmpmat = [tmppremcr tmppostmcr tmpcvrmcr];
        
        regio_mcr = [regio_mcr; tmpmat];
        
        count = count + 1;
        case_sep = [case_sep; count, count+length(regionum_inter)-1, casenum];
        count = count + length(regionum_inter)-1;
        
        allregionum = [allregionum; regionum_inter'];
        [~,tmpindmr,tmpindregio] = intersect(regionum_inter,regionum);
        allregiolbl = [allregiolbl; [regions(regioind).regionlist(tmpindregio,1)]];
        
        %regionum_inter = intersect(regionum_mr,regionum_mcr);
        
        % Semi-circular... (L and R per slice)
        %         tmpcount2 = count2 + 1;
        %         for jj = 1:size(regions(regioind).hemislist,1)
        %             tmpind = ismember(regionum_inter,[regions(regioind).hemislist{jj,3}]);
        %             if sum(tmpind)
        %             hemis_mr = [hemis_mr;...
        %                 mean(tmppremr(tmpind)) mean(tmppostmr(tmpind)) mean(tmpcvrmr(tmpind))];
        %             hemis_mcr = [hemis_mcr;...
        %                 mean(tmppremcr(tmpind)) mean(tmppostmcr(tmpind)) mean(tmpcvrmcr(tmpind))];
        %             count2 = count2 + 1;
        %             end
        %         end
        %         case_sep2 = [case_sep2; tmpcount2, count2, casenum];
        
        % Hemispheric
        tmpcount2 = count2 + 1;
        hemistypes = unique(regions(regioind).hemislist(:,2)); % L or R
        for jj = 1:length(hemistypes)
            hemisind = strcmpi(hemistypes{jj},regions(regioind).hemislist(:,2));
            tmpind = ismember(regionum_inter,[regions(regioind).hemislist{hemisind,3}]);
            if sum(tmpind)
                hemis_mr = [hemis_mr;...
                    mean(tmppremr(tmpind)) mean(tmppostmr(tmpind)) mean(tmpcvrmr(tmpind))];
                hemis_mcr = [hemis_mcr;...
                    mean(tmppremcr(tmpind)) mean(tmppostmcr(tmpind)) mean(tmpcvrmcr(tmpind))];
                count2 = count2 + 1;
            end
        end
        case_sep2 = [case_sep2; tmpcount2, count2, casenum];
        
        allhemislbl = [allhemislbl; hemistypes];
    end
end
            
vout = {regio_mr,regio_mcr, hemis_mr,hemis_mcr, case_sep,case_sep2, allregionum,allhemislbl,allregiolbl};

end

