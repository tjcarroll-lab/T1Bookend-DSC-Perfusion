function [ret] = make_3_by_3_qCBV (images, image_names, patient_number,path,type)
%-------------------------------------------------------------------%
%                                                                   %
%                                                                   %
%-------------------------------------------------------------------%

raw_qCBV  = images{strmatch('qCBV_DSC',image_names)};

[n,m]=size(raw_qCBV(:,:,1));
[n,m,nslices]= size(raw_qCBV);

total_pages = ceil(nslices/9);
    
%-------------------------------------------------------------------%
%  resize to 128 x 128                                              %
%-------------------------------------------------------------------%
if ( n > 128)
    qCBV = zeros(128,128,nslices);
    for slice = 1:nslices
        rescale = 128/n;
        qCBV(:,:,slice) = imresize(raw_qCBV(:,:,slice),rescale,'bicubic');
    end
end
if ( n == 128)
    qCBV = raw_qCBV;
end

%-----------------------------------------------------------------%
% mean                                                            %
%-----------------------------------------------------------------%
mean_qCBV = 0.0;
n_voxels = 0.0;
for i=1:128;
for j=1:128
for k=1:nslices
   if (qCBV(i,j,k) >0.001)
      mean_qCBV = mean_qCBV + qCBV(i,j,k);
      n_voxels = n_voxels + 1;
   end
end
end
end

mean_qCBV = mean_qCBV/n_voxels

std_qCBV = 0.0;
n_voxels = 0.0;
for i=1:128;
for j=1:128
for k=1:nslices
   if (qCBV(i,j,k) >0.001)
   std_qCBV = std_qCBV+(qCBV(i,j,k)-mean_qCBV)*(qCBV(i,j,k)-mean_qCBV);
   n_voxels = n_voxels + 1;
   end
end
end
end

std_qCBV = sqrt(std_qCBV/n_voxels)

%-----------------------------------------------------------------%
% calculate the mean value                                        %
%-----------------------------------------------------------------%

image_panel = zeros(384,384);

position.first = [1 129 257 ];
position.last  = [128 256 384];

beenby = 0;
for slice = 1:nslices
    
    loc_slice = mod( slice , 10);
    
    if ( slice > 9.5) 
        loc_slice = loc_slice + 1;
    end
   
    column = mod((loc_slice), 3);
    if( column == 0)
        column = 3;
    end
    row    = floor( (loc_slice-0.5)/3 )+1;

   page = floor(slice/9) + 1;

   image_panel( position.first(row):position.last(row), position.first(column):position.last(column)) = qCBV(:,:,slice);

   %-----------------------------------------------------------------------
   % print page 1                                                         % 
   %-----------------------------------------------------------------------
   if ( page == 2 )
       if (beenby == 0 )
         h = figure(1);
         imshow(image_panel,[0 10]);colorbar, colormap jet; title( [ ' CBV ml/100g ', patient_number, ...
                 ' Page ', int2str(page-1), ' of ',  int2str(total_pages)]); 
         text(100, 400, [' Mean ', num2str(mean_qCBV), ' +/- ', num2str(std_qCBV)]);
%          filename = [patient_number, ' Page ', int2str(page-1),'.tif'];
         filename = [path '\Maps\' patient_number '\' patient_number type '_qCBV' int2str(page-1) '.tif'];
         saveas(h, filename);
         image_panel = zeros(384,384);
         beenby = 1;
       end
   end
end 


%-----------------------------------------------------------------------
% print page 2                                                         % 
%-----------------------------------------------------------------------
h = figure(2);
imshow(image_panel,[0 10]);colorbar, colormap jet; title( [ ' CBV ml/100g ', patient_number, ...
                 ' Page ', int2str(page), ' of ',  int2str(total_pages)]); 
         text(100, 400, [' Mean ', num2str(mean_qCBV), ' +/- ', num2str(std_qCBV)]);
% filename = [patient_number, ' Page ', int2str(page),'.tif'];
filename = [path '\Maps\' patient_number '\' patient_number type '_qCBV' int2str(page) '.tif'];
saveas(h, filename);

return