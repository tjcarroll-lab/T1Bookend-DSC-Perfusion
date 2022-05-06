% Author: Calden Wloka (calden.wloka@alumni.utoronto.ca)
% Function to rename the result files from the Book-End Process with the
% patient MRN.
function renameResults(pNums)

pFiles = cell(length(pNums),1);
for i = 1:length(pNums)
    if pNums(i) < 10
       pFiles{i} = ['P00' num2str(pNums(i))];
    elseif pNums(i)<100
       pFiles{i} = ['P0' num2str(pNums(i))];
    else
       pFiles{i} = ['P' num2str(pNums(i))];
    end
end

resPath = 'Result_MSwcf2\';
for i = 1:length(pNums)
    resName = strcat(pFiles{i}, 'GE_M.mat');
    if(exist(strcat(resPath, resName), 'file') ~= 0)
        dInfo = dicominfo(strcat(pFiles{i}, '\ep2d_perf\1.dcm'), 'dictionary', 'dicom-dict.txt');
        mrn = dInfo.PatientID;
        copyfile(strcat(pFiles{i}, '\ep2d_perf\1.dcm'), strcat('ResultsMRN\', mrn, '.dcm'));
        movefile(strcat(resPath, resName), strcat('ResultsMRN\', mrn, '.mat'));
    else
        warning(strcat('No results for patient ', pFiles{i}, '- patient skipped'))
    end
end