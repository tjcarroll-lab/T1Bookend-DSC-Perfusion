function [ vout ] = getAllPerf( perfmr, perfmcr, regions )
%getAllPerf Combine MRI and microsphere cbf values regional and hemispheric
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2017-07-30

casenum_mr = [perfmr.casenum];
casenum_mcr = [perfmcr.casenum];
casenum_regions = [regions.casenum];

regio_mr = [];
regio_mcr = [];
hemis_mr = [];
hemis_mcr = [];
allregionum = [];
allhemislbl = [];
allregiolbl = [];

count = 0;
count2 = 0;
case_sep = [];
case_sep2 = [];

for ii = 1:length(perfmr)
    casenum = perfmr(ii).casenum;
    mcrind = ismember(casenum_mcr,casenum);
    regioind = ismember(casenum_regions,casenum);
    if sum(mcrind) && sum(regioind)
        regionum_mr = [perfmr(ii).cbfpre{:,1}];
        regionum_mcr = [perfmcr(mcrind).cbfpre{:,1}];
        regionum = [regions(regioind).regionlist{:,5}];
        
        [regionum_inter,tmpindmr,tmpindmcr] = intersect(regionum_mr,regionum_mcr); %returns sorted values, indices
        tmppremr = [perfmr(ii).cbfpre{tmpindmr,2}]';
        tmpmat = [tmppremr];
                  
        regio_mr = [regio_mr; tmpmat]; 
        
        tmppremcr = [perfmcr(mcrind).cbfpre{tmpindmcr,2}]'.*100;
        tmpmat = [tmppremcr];
        
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
                    mean(tmppremr(tmpind))];
                hemis_mcr = [hemis_mcr;...
                    mean(tmppremcr(tmpind))];
                count2 = count2 + 1;
            end
        end
        case_sep2 = [case_sep2; tmpcount2, count2, casenum];
        
        allhemislbl = [allhemislbl; hemistypes];
    end
end
            
vout = {regio_mr,regio_mcr, hemis_mr,hemis_mcr, case_sep,case_sep2, allregionum,allhemislbl,allregiolbl};

end

