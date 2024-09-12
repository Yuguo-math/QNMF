function  [Init_Index]= Q_Block_matching( X, par, Neighbor_arr, Num_arr, SelfIndex_arr)
L         =   length(Num_arr);
Init_Index   =  zeros(par.patnum,L);

for  i  =  1 : L
    Patch = X(:, SelfIndex_arr(i),:);
    Neighbors = X(:, Neighbor_arr(1:Num_arr(i),i), :);

    Dist = sum( sum((repmat( Patch, 1, size( Neighbors, 2))- Neighbors).^ 2), 3);
    
    % 只RGB
%     Patch_rgb = Patch;
%     Patch_rgb(:,:,1) = 0;
%     
%     Neighbors_rgb = Neighbors;
%     Neighbors_rgb(:,:,1) = 0;
%     Dist = sum( sum((repmat( Patch_rgb, 1, size( Neighbors_rgb, 2))- Neighbors_rgb).^ 2), 3);

     %只Y
%     Patch_y = Patch;
%     Patch_y(:,:,2:4) = 0.2*Patch_y(:,:,2:4);
%     
%     Neighbors_y = Neighbors;
%     Neighbors_y(:,:,2:4) = 0.2*Neighbors_y(:,:,2:4);
%     Dist = sum( sum((repmat( Patch_y, 1, size( Neighbors_y, 2))- Neighbors_y).^ 2), 3);


    [val, index]= sort( Dist);

    Init_Index(:,i)=Neighbor_arr(index(1: par.patnum), i);
end
