function [mean_qCBF] = make_3_by_3_qCBF (images, image_names, patient_number,path,type)
%-------------------------------------------------------------------%
%                                                                   %
%                                                                   %
%-------------------------------------------------------------------%

raw_qCBF  = images{strmatch('qCBF_nSVD',image_names)};

[n,m]=size(raw_qCBF(:,:,1));
[n,m,nslices]= size(raw_qCBF);

total_pages = ceil(nslices/9);

%-------------------------------------------------------------------%
%  resize to 128 x 128                                              %
%-------------------------------------------------------------------%
qCBF = zeros(128,128,nslices);
if ( n ~= 128)
    for slice = 1:nslices
        rescale = 128/n;
        qCBF(:,:,slice) = imresize(raw_qCBF(:,:,slice),rescale,'bicubic');
    end
else
    qCBF = raw_qCBF;
end

if ( n < 128)
    for slice = 1:nslices
        rescale = 128/n;
        qCBF(:,:,slice) = imresize(raw_qCBF(:,:,slice),rescale,'bicubic');
    end
end

%-----------------------------------------------------------------%
% mean                                                            %
%-----------------------------------------------------------------%
mean_qCBF = 0.0;
n_voxels = 0.0;
for i=1:128;
    for j=1:128
        for k=1:nslices
            if (qCBF(i,j,k) >0.001)
                mean_qCBF = mean_qCBF + qCBF(i,j,k);
                n_voxels = n_voxels + 1;
            end
        end
    end
end
mean_qCBF = mean_qCBF/n_voxels

std_qCBF = 0.0;
n_voxels = 0.0;
for i=1:128;
    for j=1:128
        for k=1:nslices
            if (qCBF(i,j,k) >0.001)
                std_qCBF = std_qCBF + (qCBF(i,j,k)-mean_qCBF)*(qCBF(i,j,k)-mean_qCBF);
                n_voxels = n_voxels + 1;
            end
        end
    end
end
std_qCBF = sqrt(std_qCBF/n_voxels)

image_panel = zeros(384,384);

position.first = [1 129 257 ];
position.last  = [128 256 384];

beenby = 0;
for slice = 1:nslices

    %     page = ceil(slice/9)
    %     row = ceil((slice/9
    %     column =
    %
    %     loc_slice = mod( slice , 10);
    %
    %     if ( slice > 9.5)
    %         loc_slice = loc_slice + 1;
    %     end
    %
    %     column = mod((loc_slice), 3);
    %     if( column == 0)
    %         column = 3;
    %     end
    %
    %     row    = floor( (loc_slice-0.5)/3 )+1;
    %
    %     page = floor(slice/9) + 1;

    page = floor(slice/9)+1;

    row = ceil(round(100*9*(slice/9-floor(slice/9))/100)/3);

    if row ==0;
        row = 3;
    end

    column = round(floor(10*3*(slice/3-floor(slice/3)))/10);

    if column == 0
        column = 3;
    end

    image_panel( position.first(row):position.last(row), position.first(column):position.last(column)) = qCBF(:,:,slice);

    %-----------------------------------------------------------------------
    % print page 1                                                         %
    %-----------------------------------------------------------------------
    if ( page == 2 )
        if (beenby == 0 )
            h = figure(1);
            imshow(image_panel,[0 120]);colorbar, colormap jet; title( [ ' CBF ml/100g/min ', patient_number, ...
                ' Page ', int2str(page-1), ' of ',  int2str(total_pages)]);
            %          filename = [patient_number, ' Page ', int2str(page-1),'.tif'];
            filename = [path '\Maps\' patient_number '\' patient_number type '_qCBF' int2str(page-1) '.tif'];
            text(100, 400, [' Mean ', num2str(mean_qCBF), ' +/- ', num2str(std_qCBF)]);
            saveas(h, filename);
            image_panel = zeros(384,384);
            beenby = 1;
        end
    end

    %-----------------------------------------------------------------------
    % print page 2                                                         %
    %-----------------------------------------------------------------------

    if ( page == 3 )
        if (beenby == 1 )
            h = figure(2);
            imshow(image_panel,[0 120]);colorbar, colormap jet; title( [ ' CBF ml/100g/min ', patient_number,  ...
                ' Page ', int2str(page-1), ' of ',  int2str(total_pages)]);
            % filename = [patient_number, ' Page ', int2str(page),'.tif'];
            filename = [path '\Maps\' patient_number '\' patient_number type '_qCBF' int2str(page-1) '.tif'];
            text(100, 400, [' Mean ', num2str(mean_qCBF), ' +/- ', num2str(std_qCBF)]);
%             filename = [path '\Maps\' patient_number '\' patient_number type '_qCBF' int2str(page) '.tif'];
            saveas(h, filename);
            image_panel = zeros(384,384);
            beenby = 2;
        end
    end

end

%-----------------------------------------------------------------------
% print page 3                                                         %
%-----------------------------------------------------------------------
h = figure(3);
imshow(image_panel,[0 120]);colorbar, colormap jet; title( [ ' CBF ml/100g/min ', patient_number,  ...
    ' Page ', int2str(page), ' of ',  int2str(total_pages)]);
% filename = [patient_number, ' Page ', int2str(page),'.tif'];
filename = [path '\Maps\' patient_number '\' patient_number type '_qCBF' int2str(page) '.tif'];
saveas(h, filename);

return