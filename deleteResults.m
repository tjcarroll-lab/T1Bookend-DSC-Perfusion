% Author: Calden Wloka (calden.wloka@alumni.utoronto.ca)
% Function to delete the result files from the Book-End Process for results
% that were no good and have to be redone.
function deleteResults(pNums)

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
        delete(strcat(resPath, resName));
    else
        warning(strcat('No results for patient ', pFiles{i}, '- patient skipped'))
    end
end