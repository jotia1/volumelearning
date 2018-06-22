function [ layer_num, n_id ] = idx2layerid( l_sizes, n )
    if n > sum(l_sizes)
        error('idx %d, out of range of layers', n); 
    end
    
    % TODO hardcoded, properly write these function for general case
    layer_num = 1;
    n_id = n;
    if n > l_sizes(1)
        layer_num = 2;
        n_id = n - l_sizes(1);
    end
    
    
%     layer_num = 1;
%     n_id = n;
%     for i = 1 : numel(l_sizes)
%         num_in_layer = l_sizes(i);
%         if n <= num_in_layer
%             break
%         end
%         layer_num = layer_num + 1;
%         n = n - num_in_layer;
%         n_id = n;
%     end
end

