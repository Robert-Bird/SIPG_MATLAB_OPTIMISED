function [no_blocks,nelblo_all]=block_size(no_calcs,bloc_size)
no_blocks = no_calcs/bloc_size;                                              % Block sizes for each block loop.
nelblo_all = zeros(ceil(no_blocks),1);                                       
nelblo_all(ceil(no_blocks)) = int32((no_blocks-floor(no_blocks))*bloc_size); % Last block size
nelblo_all(1:floor(no_blocks)) = bloc_size;                                  % Defines block size up to or including the last loop
no_blocks = ceil(no_blocks);                                                 % Defines number of block loops
end