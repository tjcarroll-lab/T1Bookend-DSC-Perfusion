function [ret] = make_3_by_3_MTT_CVP (images, image_names, patient_number,path,type)
%-------------------------------------------------------------------%
%                                                                   %
%                                                                   %
%-------------------------------------------------------------------%

raw_MTT  = images{strmatch('CMTT_CVP',image_names)};

[n,m]=size(raw_MTT(:,:,1));
[n,m,nslices]= size(raw_MTT);

total_pages = ceil(nslices/9);
    
%-------------------------------------------------------------------%
%  resize to 128 x 128                                              %
%-------------------------------------------------------------------%
if ( n ~= 128)
    qCBV = zeros(128,128,nslices);
    for slice = 1:nslices
        rescale = 128/n;
        MTT(:,:,slice) = imresize(raw_MTT(:,:,slice),rescale,'bicubic');
    end
end
if ( n == 128)
    MTT = raw_MTT;
end


%-----------------------------------------------------------------%
% mean                                                            %
%-----------------------------------------------------------------%
mean_MTT = 0.0;
n_voxels = 0.0;
for i=1:128;
for j=1:128
for k=1:nslices
   if (MTT(i,j,k) > 0.1) && (MTT(i,j,k) < 50)
      mean_MTT = mean_MTT + MTT(i,j,k); 
      n_voxels = n_voxels + 1;
   end
end
end
end
mean_MTT = mean_MTT/n_voxels
            
std_MTT = 0.0;
n_voxels = 0.0;
for i=1:128;
for j=1:128
for k=1:nslices
   if (MTT(i,j,k) >0.1) && (MTT(i,j,k) < 50)
   std_MTT = std_MTT + (MTT(i,j,k)-mean_MTT)*(MTT(i,j,k)-mean_MTT);
   n_voxels = n_voxels + 1;
   end
end
end
end

std_MTT = sqrt(std_MTT/n_voxels)

%-----------------------------------------------------------------%
% calculate the mean value                                        %
%-----------------------------------------------------------------%

image_panel = zeros(384,384);

position.first = [1 129 257 ];
position.last  = [128 256 384];

beenby = 0;
for slice = 1:nslices
    
    page = floor(slice/9)+1;

    row = ceil(round(100*9*(slice/9-floor(slice/9))/100)/3);

    if row ==0;
        row = 3;
    end

    column = round(floor(10*3*(slice/3-floor(slice/3)))/10);

    if column == 0
        column = 3;
    end

   image_panel( position.first(row):position.last(row), position.first(column):position.last(column)) = MTT(:,:,slice);

        %-----------------------------------------------------------------------
    % print page 1                                                         %
    %-----------------------------------------------------------------------
    if ( page == 2 )
        if (beenby == 0 )
            h = figure(1);
            imshow(image_panel,[0 10]);colorbar, colormap jet; title( [ ' MTT sec ', patient_number, ...
                ' Page ', int2str(page-1), ' of ',  int2str(total_pages)]);
            %          filename = [patient_number, ' Page ', int2str(page-1),'.tif'];
            filename = [path '\Maps\' patient_number '\' patient_number type '_MTT' int2str(page-1) '.tif'];
            text(100, 400, [' Mean ', num2str(mean_MTT), ' +/- ', num2str(std_MTT)]);
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
            imshow(image_panel,[0 10]);colorbar, colormap jet; title( [ ' MTT sec ', patient_number,  ...
                ' Page ', int2str(page-1), ' of ',  int2str(total_pages)]);
            % filename = [patient_number, ' Page ', int2str(page),'.tif'];
            filename = [path '\Maps\' patient_number '\' patient_number type '_MTT' int2str(page-1) '.tif'];
            text(100, 400, [' Mean ', num2str(mean_MTT), ' +/- ', num2str(std_MTT)]);
%             filename = [path '\Maps\' patient_number '\' patient_number type '_MTT' int2str(page) '.tif'];
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
imshow(image_panel,[0 10]);colorbar, colormap jet; title( [ ' MTT sec ', patient_number,  ...
    ' Page ', int2str(page), ' of ',  int2str(total_pages)]);
% filename = [patient_number, ' Page ', int2str(page),'.tif'];
filename = [path '\Maps\' patient_number '\' patient_number type '_MTT' int2str(page) '.tif'];
saveas(h, filename);

return