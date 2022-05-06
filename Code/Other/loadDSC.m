% Load DSC dicoms
function dsc = loadDSC(dscpath)

dscpath = checkFileSep(dscpath);
dcmdir = dir(fullfile(dscpath,'*dcm'));

tmpdsc = [];
tmpdsc_stack = [];
f = waitbar(0,'Loading DSC');
tmphdr = dicominfo(fullfile(dscpath,'1.dcm'));

if isfield(tmphdr,'NumberOfTemporalPositions')
    ntp = tmphdr.NumberOfTemporalPositions;
    N_slices = length(dcmdir)/ntp;
else
    for ii = 1:length(dcmdir)
        tmphdr = dicominfo(fullfile(dscpath, [num2str(ii) '.dcm']));
        slclist(ii) = tmphdr.SliceLocation;
    end
    N_slices = length(unique(slclist));
    ntp = length(dcmdir)/N_slices;
end

%ntp = tmphdr.NumberOfTemporalPositions;
for ii = 1:length(dcmdir)
    tmpdsc = cat(3,tmpdsc,dicomread(fullfile(dscpath,[num2str(ii) '.dcm'])));
    if mod(ii,ntp) == 0
        tmpdsc_stack = cat(4,tmpdsc_stack,tmpdsc);
        tmpdsc = [];
    end
    waitbar(ii/length(dcmdir),f);
end
dsc = tmpdsc_stack;
close(f);