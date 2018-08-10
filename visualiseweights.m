%% VISUALISEWEIGHTS - Creates visualisation of the current network state
%   Assumes 'usual' network variables are already in the work space and
%   visualises them 

% network input
num_rows = 3;
num_cols = 4;
subplot(num_rows, 4, 1);
%locsx = xs(ts > sec * 1000 & ts <= ((sec + 1) * 1000));
%locsy = ys(ts > sec * 1000 & ts <= ((sec + 1) * 1000));
[rows, cols] = ind2sub([net.inp_img_size, net.inp_img_size], inp(ts > sec * 1000 & ts <= ((sec + 1) * 1000)));
plot3(cols, rows, ts(ts > sec * 1000 & ts <= ((sec + 1) * 1000)), '.k');
axis([ 0 net.inp_img_size 0 net.inp_img_size -Inf Inf]);
ax = gca;
set(ax, 'Ydir', 'reverse');
xlabel('cols');
ylabel('rows');
zlabel('time (ms)');
%drawweightedlocs(locsx, locsy, ones(size(locsx)), net.inp_img_size);
title(sprintf('Input at %d seconds', sec));

%% dendrite weights
subplot(num_rows, 4, 2);
[locsx, locsy] = ind2sub([net.inp_img_size, net.inp_img_size], pre_dend(neuron_to_plot, :));
drawweightedlocs(locsx, locsy, w_dend(1, :), net.inp_img_size);
title(sprintf('Dendrite weights of N%d', neuron_to_plot));

%% axon weights
subplot(num_rows, 4, 6);
[locsx, locsy] = ind2sub([net.inp_img_size, net.inp_img_size], post_axon(neuron_to_plot, :) - net.N_inp - net.N_hid);
drawweightedlocs(locsx, locsy, w_axon(neuron_to_plot, :), net.inp_img_size);
title(sprintf('Axon weights of N%d', neuron_to_plot));

%%  general connectivity
subplot(num_rows*2, 4, num_cols * 2 + 1);
[locsx, locsy] = ind2sub([net.inp_img_size, net.inp_img_size], pre_dend(neuron_to_plot, :));
drawlocs(locsx, locsy, net.inp_img_size);
title('Dendrite connections');

subplot(num_rows * 2, 4, num_cols * 3 + 1);
[locsx, locsy] = ind2sub([net.inp_img_size, net.inp_img_size], post_axon(neuron_to_plot, :) - net.N_inp - net.N_hid);
drawlocs(locsx, locsy, net.inp_img_size);
title('Axon connections');

% Axon variance
subplot(num_rows, 4, 8);
[locsx, locsy] = ind2sub([net.inp_img_size, net.inp_img_size], post_axon(neuron_to_plot, :) - net.N_inp - net.N_hid);
drawweightedlocs(locsx, locsy, variance_axon(neuron_to_plot, :), net.inp_img_size);
title(sprintf('Axon variance of N%d', neuron_to_plot));

% dendrite variance 
subplot(num_rows, 4, 4);
[locsx, locsy] = ind2sub([net.inp_img_size, net.inp_img_size], pre_dend(neuron_to_plot, :));
drawweightedlocs(locsx, locsy, variance_dend(neuron_to_plot, :), net.inp_img_size);
title(sprintf('Dendrite variance of N%d', neuron_to_plot));


%% dendrite delays
subplot(num_rows, 4, 3);
[locsx, locsy] = ind2sub([net.inp_img_size, net.inp_img_size], pre_dend(neuron_to_plot, :));
drawweightedlocs(locsx, locsy, delays_dend(neuron_to_plot, :), net.inp_img_size);
title(sprintf('Dendrite delays of N%d', neuron_to_plot));

%% axon delays
subplot(num_rows, 4, 7);
[locsx, locsy] = ind2sub([net.inp_img_size, net.inp_img_size], post_axon(neuron_to_plot, :) - net.N_inp - net.N_hid);
drawweightedlocs(locsx, locsy, delays_axon(neuron_to_plot, :), net.inp_img_size);
title(sprintf('Axon delays of N%d', neuron_to_plot));
drawnow;

%% helper functions
function [] = drawweightedlocs( locsx, locsy, weightings, inp_img_size )

    pts = linspace(1, inp_img_size + 1, inp_img_size + 1) - 0.5;
    counts = histcounts2(locsx(:), locsy(:), pts, pts );
    weighted = accumarray([locsx(:), locsy(:)], weightings(:));
    padded = zeros(inp_img_size, inp_img_size);
    padded(find(weighted)) = weighted(find(weighted));
    imagesc(padded);
    axis equal;
    set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]));
    colorbar
    
end

function [] = drawlocs( locsx, locsy, inp_img_size )

    pts = linspace(1, inp_img_size + 1, inp_img_size + 1) - 0.5;
    counts = histcounts2(locsx(:), locsy(:), pts, pts );
    imagesc(counts);
    axis equal;
    set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]));
    colorbar

end
