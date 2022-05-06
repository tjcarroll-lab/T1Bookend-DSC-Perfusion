function [ dellist ] = del_perfResults( bp, flag_del )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    flag_del = 0;
end
%pickid = {'10','11'};
%pickid = eval(['{' sprintf('''%02d'',',1:30) '}']);
pickid = {'04','05','06','12','21','23','26','30'};
%pickid = {'07','13','14','15','17','19','22','24'};
pickdate = {'d02'};
%pickresult = {'Result_MSwcf2'};
%pickresult = {'Result_MSwcf2','DSCanalysis','AIFdata'};
%pickresult = {'AutoSS','Result_MSwcf2'};
pickresult = {'DSCanalysis','Result_MSwcf2'};
%pickresult = {'AIFdata'};
%pickresult = {'Result_MSwcf2','DSCanalysis','AIFdata','T1mapping','AutoSS'};    
%pickperf = {'P001GE_DSC_A.mat','P001GE_M.mat'};%,'P002'};
pickperf = {'P001'};

bp = checkFileSep(bp);
filelist = searchFolders(bp,'.mat',1);
dellist = {};

for ii = 1:length(filelist)
    currpath = checkFileSep(filelist{ii});
    trailingpath = currpath(length(bp)+1:end);
    slashind = strfind(trailingpath,'\');
    if length(slashind) >= 4
        caseid = trailingpath(1:slashind(1)-1);
        casedate = trailingpath(slashind(1)+1:slashind(2)-1);
        caseresult = trailingpath(slashind(2)+1:slashind(3)-1);
        caseperf = trailingpath(slashind(3)+1:slashind(4)-1);
        if ~isempty(regexpi(caseid,'^COLLAT_')) && sum(strcmpi(caseid(8:end),pickid)) &&...
                sum(strcmpi(casedate(1:3),pickdate))
            if sum(strcmpi(caseresult,pickresult)) && ( sum(strcmpi(caseperf,pickperf)) || sum(strcmpi(caseperf(1:4),pickperf)) )
                if flag_del
                    delete(currpath(1:end-1));
                end
                fprintf('%s\n',['Deleted ' currpath(1:end-1)]);
                dellist = [dellist; currpath(1:end-1)];
            end
        end
    end
end

end

