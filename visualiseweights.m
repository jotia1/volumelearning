%% VISUALISEWEIGHTS - Creates visualisation of the current network state
%   Assumes 'usual' network variables are already in the work space and
%   visualises them 

neuron_to_plot = 1;

% Visualise network input
subplot(2, 4, 1);
locsx = xs(ts > sec * 1000 & ts <= ((sec + 1) * 1000));
locsy = ys(ts > sec * 1000 & ts <= ((sec + 1) * 1000));
drawweightedlocs(locsx, locsy, ones(size(locsx)), im_size);
title(sprintf('Input at %d seconds', sec));

%% Visualise dendrite weights
subplot(2, 4, 2);
[locsx, locsy] = ind2sub([im_size, im_size], pre_hid(neuron_to_plot, :));
drawweightedlocs(locsx, locsy, w_hid(1, :), im_size);
title(sprintf('Dendrite weights of N%d', neuron_to_plot));

%% Visualise axon weights
subplot(2, 4, 6);
[locsx, locsy] = ind2sub([im_size, im_size], post_out(neuron_to_plot, :) - N_inp - N_hid);
drawweightedlocs(locsx, locsy, w_out(neuron_to_plot, :), im_size);
title(sprintf('Axon weights of N%d', neuron_to_plot));

%% Visualise general connectivity
subplot(4, 4, 9);
[locsx, locsy] = ind2sub([im_size, im_size], pre_hid(neuron_to_plot, :));
drawlocs(locsx, locsy, im_size);
title('Dendrite connections');

subplot(4, 4, 13);
[locsx, locsy] = ind2sub([im_size, im_size], post_out(neuron_to_plot, :) - N_inp - N_hid);
drawlocs(locsx, locsy, im_size);
title('Axon connections');

% Visualise dendrite variance
subplot(2, 4, 8);
[locsx, locsy] = ind2sub([im_size, im_size], post_out(neuron_to_plot, :) - N_inp - N_hid);
drawweightedlocs(locsx, locsy, variance_out(neuron_to_plot, :), im_size);
title(sprintf('Dendrite variance of N%d', neuron_to_plot));

% Visualise axon variance 
subplot(2, 4, 4);
[locsx, locsy] = ind2sub([im_size, im_size], pre_hid(neuron_to_plot, :));
drawweightedlocs(locsx, locsy, variance_hid(neuron_to_plot, :), im_size);
title(sprintf('Dendrite variance of N%d', neuron_to_plot));


%% Visualise dendrite delays
subplot(2, 4, 3);
[locsx, locsy] = ind2sub([im_size, im_size], pre_hid(neuron_to_plot, :));
drawweightedlocs(locsx, locsy, delays_hid(neuron_to_plot, :), im_size);
title(sprintf('Dendrite delays of N%d', neuron_to_plot));

%% Visualise axon delays
subplot(2, 4, 7);
[locsx, locsy] = ind2sub([im_size, im_size], post_out(neuron_to_plot, :) - N_inp - N_hid);
drawweightedlocs(locsx, locsy, delays_out(neuron_to_plot, :), im_size);
title(sprintf('Axon delays of N%d', neuron_to_plot));


%% helper functions
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