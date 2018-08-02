%% VISUALISEWEIGHTS - Creates visualisation of the current network state
%   Assumes 'usual' network variables are already in the work space and
%   visualises them 


% Visualise input connections
subplot(2, 3, 1);
neuron_to_plot = 1;
[locsx, locsy] = ind2sub([im_size, im_size], pre_hid(neuron_to_plot, :));
drawweightedlocs(locsx, locsy, w_hid(1, :), im_size);
title(sprintf('Dendrite weights of N%d', neuron_to_plot));


% Visualise axon connections
subplot(2, 3, 2);
[locsx, locsy] = ind2sub([im_size, im_size], post_out(neuron_to_plot, :) - N_inp - N_hid);
drawweightedlocs(locsx, locsy, w_out(neuron_to_plot, :), im_size);
title(sprintf('Axon weights of N%d', neuron_to_plot));

% Visualise input variance
subplot(2, 3, 3);


% Visualise output delays

% Visualise 




function [] = drawweightedlocs( locsx, locsy, weightings, im_size )

    pts = linspace(1, im_size + 1, im_size + 1) - 0.5;
    counts = histcounts2(locsx(:), locsy(:), pts, pts );
    weighted = accumarray([locsx(:), locsy(:)], weightings(:));
    imagesc(weighted);
    axis equal;
    set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]));
    colorbar
    
end

function [] = drawlocs( locsx, locsy, im_size )

    pts = linspace(1, im_size + 1, im_size + 1) - 0.5;
    counts = histcounts2(locsx(:), locsy(:), pts, pts );
    imagesc(counts);
    axis equal;
    set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]));
    colorbar

end