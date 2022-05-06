addpath('.\PerfCore');

datasheet = readcell('..\EXCEL FILES\myriad case list v2.xlsx');
cases = datasheet(2:end-2, strcmpi(datasheet(1,:),'caseid'));
pickcases = {'01-001'};
aifcheck = datasheet(2:end-2, strcmpi(datasheet(1,:),'aif'));
vofcheck = datasheet(2:end-2, strcmpi(datasheet(1,:),'vein'));
tmpcheck{1} = 1;
dstbp = 'D:\Users\CarrollLab\Desktop\Yongs Research\MYRIAD TO PROCESS';
currtime = datestr(datetime('now'),'yyyymmdd-HHMMSS');
fh = fopen(fullfile(dstbp,['perflog_' currtime '.txt']),'a');
%fh = fopen(fullfile(dstbp,['test_perflog_' currtime '.txt']),'a');

%%
t0 = tic;
for ii = 1:length(cases)
    if ~contains(cases{ii}, pickcases)
        continue;
    end

    dstbp1 = [dstbp '\' cases{ii}];
    dstdir = dir(dstbp1);

    for jj = 3:length(dstdir)
        dstpath = [dstbp1 '\' dstdir(jj).name];
        if contains(dstdir(jj).name,'')
            if aifcheck{ii} == 1 && vofcheck{ii} == 1 && ...
                    exist(fullfile(dstpath,'AIF_Mask_P001GE_M.mat'),'file') &&...
                    exist(fullfile(dstpath,'Vein_Mask_P001GE_M.mat'),'file')
                disp(dstpath);
                if exist(fullfile(dstpath,'DSCanalysis','P001GE_DSC_A.mat'),'file')
                    fprintf(fh, 'Already done %s\n', [cases{ii} '\' dstdir(jj).name]);
                    continue;
                end
                try
                    Auto_qCBF_Philips_ManualWMmask_MSwcf2_DD_all(dstpath);
                    fprintf(fh, 'Processed %s\n', [cases{ii} '\' dstdir(jj).name]);
                catch err
                    fprintf(fh, '**Failed %s (%s)\n', [cases{ii} '\' dstdir(jj).name], err.getReport('extended','hyperlinks','off'));
                end
            end
        end
    end
end
toc(t0);
fclose(fh);