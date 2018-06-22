function [ idx ] = layerid2idx( l_sizes, layer_num, n_id )
%@(l_sizes, layer_num, n_id) sum(l_sizes(1:layer_num-1)) + n_id;

    % TODO hardcoded, write properly for general case later
    idx = sum(l_sizes(1:layer_num-1)) + n_id;

%     assert(mean(n_id > 0) == 1, 'n_id must be positive int');
%     assert(mean(sum(l_sizes(:)) >= n_id) == 1, 'id %d larger than whole network', n_id);
%     assert(mean(n_id <= l_sizes(layer_num)) == 1, 'id %d larger than the layer %d (%d)', n_id, layer_num, l_sizes(layer_num));
%     
%     sums = [0, cumsum(l_sizes)];
%     idx = sums(layer_num) + n_id;

end

