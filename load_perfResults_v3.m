function [ output ] = load_perfResults_v3( bp, pickproj, pickid, pickdate, pickresult, pickperf )
%load_perfResults load perf results into struct variable
%   v3: more robust, but need to input look-up attributes
%
% Author: Yong Ik Jeong
% Date: 2020-02-27

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
projind = contains(filepaths,pickproj);
idind = contains(filepaths,pickid);
dateind = contains(filepaths,pickdate);
resultind = contains(filepaths,pickresult);
perfind = startsWith(filenames,pickperf);
allind = projind & idind & dateind & resultind & perfind;
filepaths = filepaths(allind);
filenames = filenames(allind);

% Sort numerically (i.e. 1,2,3,...10,11... instead of 1,10,11,2,3...) 
[~,sortind] = sortFiles(filepaths,'_ \');
filepaths = filepaths(sortind);
filenames = filenames(sortind);

count = 0;
output = struct([]);
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

            count = count + 1;
            output(count).id = currid;
            output(count).date = currdate;
            output(count).perf = currperf;
            tmp = load(fullfile(obp,filepaths{ii},filenames{ii}));
            try
                tmp.ROIs = rmfield(tmp.ROIs,'dsc_stack');
                tmp.ROIs = rmfield(tmp.ROIs,'sigmap');
            catch
            end
            output(count).result = tmp;
            output(count).header = dicominfo(fullfile(obp,currproj,currid,currdate,currperf(1:4),'\ep2d_perf\1.dcm'));
        end
    end
end

end

