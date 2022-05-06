function [ret] = make_3_by_3_dBAT (images, image_names, patient_number,path,type)
%-------------------------------------------------------------------%
%                                                                   %
%                                                                   %
%-------------------------------------------------------------------%

raw_dBAT  = images{strmatch('dBATmap',image_names)};

[n,m]=size(raw_dBAT(:,:,1));
[n,m,nslices]= size(raw_dBAT);

total_pages = ceil(nslices/9);

%-------------------------------------------------------------------%
%  resize to 128 x 128                                              %
%-------------------------------------------------------------------%
dBAT = zeros(128,128,nslices);
if ( n ~= 128)  
    for slice = 1:nslices
        rescale = 128/n;
        dBAT(:,:,slice) = imresize(raw_dBAT(:,:,slice),rescale,'bicubic');
    end
else
    dBAT = raw_dBAT;
end

if ( n < 128)  
    for slice = 1:nslices
        rescale = 128/n;
        dBAT(:,:,slice) = imresize(raw_dBAT(:,:,slice),rescale,'bicubic');
    end
end

%-----------------------------------------------------------------%
% mean                                                            %
%-----------------------------------------------------------------%
mean_dBAT = 0.0;
n_voxels = 0.0;
for i=1:128;
for j=1:128
for k=1:nslices
   if (dBAT(i,j,k) >0.001)
      mean_dBAT = mean_dBAT + dBAT(i,j,k);
      n_voxels = n_voxels + 1;
   end
end
end
end
mean_dBAT = mean_dBAT/n_voxels
            
std_dBAT = 0.0;
n_voxels = 0.0;
for i=1:128;
for j=1:128
for k=1:nslices
   if (dBAT(i,j,k) >0.001)
   std_dBAT = std_dBAT + (dBAT(i,j,k)-mean_dBAT)*(dBAT(i,j,k)-mean_dBAT);
   n_voxels = n_voxels + 1;
   end
end
end
end
std_dBAT = sqrt(std_dBAT/n_voxels)

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

   image_panel( position.first(row):position.last(row), position.first(column):position.last(column)) = dBAT(:,:,slice);

   %-----------------------------------------------------------------------
   % print page 1                                                         % 
   %-----------------------------------------------------------------------
   if ( page == 2 )
       if (beenby == 0 )
         h = figure(1);
         imshow(image_panel,[0 4]);colorbar, colormap jet; title( [ ' dBAT sec ', patient_number, ...
                 ' Page ', int2str(page-1), ' of ',  int2str(total_pages)]);  
%          filename = [patient_number, ' Page ', int2str(page-1),'.tif'];
         filename = [path '\Maps\' patient_number '\' patient_number type '_dBAT' int2str(page-1) '.tif'];
         text(100, 400, [' Mean ', num2str(mean_dBAT), ' +/- ', num2str(std_dBAT)]);
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
imshow(image_panel,[0 4]);colorbar, colormap jet; title( [ ' dBAT sec ', patient_number,  ...
                 ' Page ', int2str(page-1), ' of ',  int2str(total_pages)]); 
% filename = [patient_number, ' Page ', int2str(page),'.tif'];
filename = [path '\Maps\' patient_number '\' patient_number type '_dBAT' int2str(page) '.tif'];
saveas(h, filename);

return