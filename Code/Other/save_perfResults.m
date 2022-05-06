function [ savelist ] = save_perfResults( bp, pickproj, pickid, pickdate, pickresult, pickperf, flag_save )
%load_perfResults load perf results into struct variable
%   v3: more robust, but need to input look-up attributes
%
% Author: Yong Ik Jeong
% Date: 2020-02-27

if ~exist('flag_save','var')
    flag_save = 0;
end
bp = checkFileSep(bp);
obp = bp;
%filelist = searchFolders(bp,'.mat',1);
filelist = dir(fullfile(bp,'**/*.mat'));

filepaths = {filelist.folder};
filenames = {filelist.name};

% Remove basepath bp from filepaths
bp = strrep(bp,'..','');
filepaths = cellfun(@(x) x(strfind(x,bp)+length(bp):end), filepaths, 'UniformOutput', false);

% Pick filepaths only containing look-up attributes
projind = contains(filepaths,pickproj,'IgnoreCase',1);
idind = contains(filepaths,pickid,'IgnoreCase',1);
dateind = contains(filepaths,pickdate,'IgnoreCase',1);
resultind = contains(filepaths,pickresult,'IgnoreCase',1);
perfind = startsWith(filenames,pickperf,'IgnoreCase',1);
allind = projind & idind & dateind & resultind & perfind;
filepaths = filepaths(allind);
filenames = filenames(allind);

if isempty(filepaths)
    fprintf('No files found.\n');
    savelist = [];
    return;
end

% Sort numerically (i.e. 1,2,3,...10,11... instead of 1,10,11,2,3...) 
[~,sortind] = sortFiles(filepaths,'_ \');
filepaths = filepaths(sortind);
filenames = filenames(sortind);

savelist = {};
for ii = 1:length(filepaths)
    
    fileparts = split(filepaths{ii},'\');
    
    % Require 4 levels of folders
    if length(fileparts) ~= 4
        continue;
    else
        currproj = fileparts{1};
        currid = fileparts{2};
        currdate = fileparts{3};
        currresult = fileparts{4};
        currperf = filenames{ii};
        
        criteria = [ startsWith(currproj,pickproj,'IgnoreCase',1) ... %1
                     endsWith(currid,pickid,'IgnoreCase',1) ... %2
                     startsWith(currdate,pickdate,'IgnoreCase',1) ... %3
                     startsWith(currresult,pickresult,'IgnoreCase',1) ... %4
                     startsWith(currperf,pickperf,'IgnoreCase',1) ]; %5
                 
        if sum(criteria) == 5
            
            if flag_save
                savepath = '..\saves\Results mat file';
                if ~exist(savepath,'dir')
                    mkdir(savepath);
                end
                copyfile(fullfile(obp,filepaths{ii},filenames{ii}),fullfile(savepath,[currproj '_' currid '_' currdate '_' currresult '_' currperf]));
            end
            fprintf('%s\n',['Copied ' fullfile(obp,filepaths{ii},filenames{ii})]);
            savelist = [savelist; fullfile(obp,filepaths{ii},filenames{ii})];
        end
    end
end

end



