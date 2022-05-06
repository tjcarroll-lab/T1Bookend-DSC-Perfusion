% Author: Calden Wloka (calden.wloka@alumni.utoronto.ca)
% Script for Book-End Processing to rename all the sequences without having
% to run in multiple folders and change values within the script.

% run in main directory with FE-EPI and Modified LL folders (all three) in
% a folder named 'R###' where '###' corresponds to the 'P###' folder used
% for the main processing.  The script is capable of operating on multiple
% folders, but Main Directory must not contain any other files or folders
% starting with 'r' or 'R' and consisting of four letters.
function renameAll()

D = dir('R*');

for i = 1:length(D)
    % note that this means no other folder or file in Main Directory can
    % have a name that is length 4 and begins with 'r' or 'R'
    if(length(D(i).name) == 4)
        % rename the folder to the folder name ready for processing:
        % R### -> P###
        subFold = strcat('P',D(i).name(2:4));
        movefile(D(i).name, subFold);
        
        % find the FE-EPI folder
        FEFold = dir(strcat(subFold, '/FE*'));
        FEFold = strcat(subFold, '/', FEFold.name);
        
        % get the file names.  The numbers are always a specific number
        % digits followed by DCM.  Since the number of digits is always 
        % the same, the numeric order of the files is preserved.
        fNames = ls(strcat(FEFold, '/*DCM'));
        
        % convert all the filenames so they run from 1.dcm to 'the maximum
        % number of images'.dcm
        for j = 1:length(fNames)
            movefile(strcat(FEFold, '/', fNames(j,:)), strcat(FEFold, '/', int2str(j), '.dcm'));
        end
        
        % rename the FE-EPI folder appropriately
        movefile(FEFold, strcat(subFold, '/ep2d_perf'));
        
        % find the Modified LL folders and determine which one is PRE and
        % which one POST
        ModLLFold = dir(strcat(subFold, '/Modified*LL*'));
        if(length(ModLLFold) ~= 2)
            error('Wrong number of modified LL sequences');
        end
        ModLLFoldA = strcat(subFold, '/', ModLLFold(1).name);
        ModLLFoldB = strcat(subFold, '/', ModLLFold(2).name);
        
        fLLA = ls(strcat(ModLLFoldA, '/*DCM'));
        fLLB = ls(strcat(ModLLFoldB, '/*DCM'));
        
        % strip off the DCM to give the numerical value for the first file
        % name in each folder
        %adjustment made by grady FOR AVIV
        valA = str2num(fLLA(1,39:end-4));
        valB = str2num(fLLB(1,39:end-4));
        
        if(valA > valB)
            movefile(ModLLFoldA, strcat(subFold, '/IR_LL_EPI_POST'));
            movefile(ModLLFoldB, strcat(subFold, '/IR_LL_EPI_PRE'));
            
            % reassign the ModLL variables to reflect the folder change
            ModLLFoldA = strcat(strcat(subFold, '/IR_LL_EPI_POST'));
            ModLLFoldB = strcat(strcat(subFold, '/IR_LL_EPI_PRE'));
        else
            movefile(ModLLFoldA, strcat(subFold, '/IR_LL_EPI_PRE'));
            movefile(ModLLFoldB, strcat(subFold, '/IR_LL_EPI_POST'));
            
            % reassign the ModLL variables to reflect the folder change
            ModLLFoldA = strcat(strcat(subFold, '/IR_LL_EPI_PRE'));
            ModLLFoldB = strcat(strcat(subFold, '/IR_LL_EPI_POST'));
        end
        
        % now that the Modified LL folders have been properly relabeled,
        % relabel the files within.  The assumption has been made that the
        % two LL sequences have the same number of files in each.  If this
        % is not the case, this for loop will need to be broken into two
        % loops.  A check has been put in just to verify this.
        if(length(fLLA) ~= length(fLLA))
            error('LL sequences of different length');
        end
        for j = 1:length(fLLA)
            movefile(strcat(ModLLFoldA, '/', fLLA(j,:)), strcat(ModLLFoldA, '/', int2str(j), '.dcm'));
            movefile(strcat(ModLLFoldB, '/', fLLB(j,:)), strcat(ModLLFoldB, '/', int2str(j), '.dcm'));
        end        
    end
end