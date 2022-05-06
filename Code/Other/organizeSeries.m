function organizeSeries(srcbp,dstbp,dsconly,overwrite)
% Organizes cases into PXXX format and LL-EPI and FE-EPI into IR_LL_EPI_PRE
% IR_LL_EPI_POST, and ep2d_perf.
%
% Author: Yong Ik Jeong
% Date: 2016-04-18
%todo add folder structure level count

if ~exist('overwrite','var')
    overwrite = 0;
end
if ~exist('dsconly','var')
    dsconly = 0;
end
srcbp = checkFileSep(srcbp);
finalpaths = searchFolders(srcbp,'.dcm');
prevtopfolder = [];
Pcount = 0;
pre = 0;
post = 0;
perf = 0;
dsc = 0;
copy = 0;
tmppaths = replace(finalpaths,{'AN','Scan','SE','d'},'_');
[~,tmpinds] = sortFiles(tmppaths,'_ \');
finalpaths = finalpaths(tmpinds);

pathlength = length(srcbp);
for ii = 1:length(finalpaths)
    trailingpath = finalpaths{ii}(pathlength:end);
    slashind = strfind(trailingpath,'\');
    currtopfolder = trailingpath(slashind(end-1)+1:slashind(end)-1);
    middlepath = trailingpath(2:slashind(end)-1);
    if ~strcmpi(currtopfolder,prevtopfolder) || isempty(prevtopfolder)
        if Pcount ~= 0 && (perf > 0 && perf < 3)
            if ~dsconly
                error('Could not find all of dsc, t1pre and t1post.');
            end
        end
        Pcount = 1;
        perf = 0;
        pre = 0;
        post = 0;
        dsc = 0;
    elseif dsconly && dsc == 1
        Pcount = Pcount + 1;
        perf = 0;
        pre = 0;
        post = 0;
        dsc = 0;
    elseif perf == 3
        Pcount = Pcount + 1;
        perf = 0;
        pre = 0;
        post = 0;
        dsc = 0;
    end
    prevtopfolder = currtopfolder;
    PN = strcat('P',sprintf('%03d',Pcount));
    currseries = trailingpath(slashind(end)+1:end);    
    if ~isempty(regexpi(currseries,'LONG[_ ]EPI[-_.+]LL[_ ]PRE')) && pre == 0 && ~dsconly
        dstpath = fullfile(dstbp,middlepath,PN,'IR_LL_EPI_PRE');
        pre = 1;
        perf = perf + 1;
        copy = 1;
    elseif ~isempty(regexpi(currseries,'FE[-_.+]EPI'))
        dstpath = fullfile(dstbp,middlepath,PN,'ep2d_perf');
        dsc = 1;
        perf = perf + 1;
        copy = 1;
    elseif ~isempty(regexpi(currseries,'SE[-_.+]EPI'))
        dstpath = fullfile(dstbp,middlepath,PN,'ep2d_se');
        dsc = 1;
        perf = perf + 1;
        copy = 1;
    elseif ~isempty(regexpi(currseries,'LONG[_ ]EPI[-_.+]LL[_ ]POST')) && post == 0 && ~dsconly
        dstpath = fullfile(dstbp,middlepath,PN,'IR_LL_EPI_POST');
        post = 1;
        perf = perf + 1;
        copy = 1;
    end
    if copy
        if ~exist(dstpath,'dir')
            mkdir(dstpath);
            skip = 0;
        else
            tmpdir = dir(dstpath);
            % Check if files already exist
            skip = sum(~[tmpdir.isdir]);
        end
            
        fprintf('%s\n%s\n',fullfile(finalpaths{ii}),fullfile(dstpath));
        
        if skip && ~overwrite
            fprintf('Skipped (files exist).\n\n');
        else
            imagedir = dir(fullfile(finalpaths{ii},'*.dcm'));
            [~,sortind] = sortFiles({imagedir.name},'_.IM');
            imagedir = imagedir(sortind);
            
            filecount = 0;
            for jj = 1:length(imagedir)
                copyfile(fullfile(finalpaths{ii},imagedir(jj).name),fullfile(dstpath,[num2str(jj) '.dcm']));
                filecount = filecount + 1;
                %fprintf('%s\n%s\n\n',fullfile(finalpaths{ii},imagedir(jj).name),fullfile(dstpath,[num2str(jj) '.dcm']));
            end
            fprintf('%d files.\n\n',filecount);
        end
    end
    copy = 0;
end

fprintf('Done.\n');

end
        