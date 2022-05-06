function images_matrix = Cells2Matrix(images_cells);

images_matrix = zeros(size(images_cells{1},1),size(images_cells{1},2),length(images_cells));

for k = 1:length(images_cells)
    images_matrix(:,:,k) = images_cells{k};
end