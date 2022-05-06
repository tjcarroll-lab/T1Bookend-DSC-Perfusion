function [curve, matrix] = pull_signal_dsc(locY,locX,slice,nslices,folder)


%load 3d time matrix

dicoms = dir([folder '/*.dcm']);
numfiles = length(dicoms);
timepoints = numfiles/nslices;
for ii = 1:timepoints-1
  %  try
    matrix(:,:,ii) = dicomread([folder '/' num2str((ii-1)*nslices+slice) '.dcm']);
%     catch
%         1;
%     end
    curve(ii) = matrix(locY,locX,ii);
end

plot(curve/curve(1))

end