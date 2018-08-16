%% VISUALIS2DEWEIGHTS - Creates visualisation of the current network state
%   Assumes 'usual' network variables are already in the work space and
%   visualises them 

% network input
num_rows = 4;
num_cols = 5;
subplot(2, num_cols, 1);
[rows, cols] = ind2sub([net.inp_img_size, net.inp_img_size], inp_trimmed);
plot3(cols, rows, ts_trimmed - (sec - 1) * 1000, '.k');
axis([ 0 net.inp_img_size 0 net.inp_img_size 0 ms_per_sec]);
ax = gca;
set(ax, 'Ydir', 'reverse');
xlabel('cols');
ylabel('rows');
zlabel('time (ms)');
%drawweightedlocs(rows, cols, ones(size(rows)), net.inp_img_size);
title(sprintf('Input at %d seconds', sec));

%%  Dendrite connectivity
subplot(num_rows, num_cols, 0 * num_cols + 2);
[rows, cols] = ind2sub([net.inp_img_size, net.inp_img_size], pre_dend(neuron_to_plot, :));
drawlocs(rows, cols, net.inp_img_size);
title('Dendrite connections');

%% Dendrite weights
subplot(num_rows, num_cols, 0 * num_cols + 3);
[rows, cols] = ind2sub([net.inp_img_size, net.inp_img_size], pre_dend(neuron_to_plot, :));
drawweightedlocs(rows, cols, w_dend(neuron_to_plot, :), net.inp_img_size);
title(sprintf('Dendrite weights of N%d', neuron_to_plot));

%% dendrite delays
subplot(num_rows, num_cols, 0 * num_cols + 4);
[rows, cols] = ind2sub([net.inp_img_size, net.inp_img_size], pre_dend(neuron_to_plot, :));
drawweightedlocs(rows, cols, delays_dend(neuron_to_plot, :), net.inp_img_size);
title(sprintf('Dendrite delays of N%d', neuron_to_plot));

% dendrite variance 
subplot(num_rows, num_cols, 0 * num_cols + 5);
[rows, cols] = ind2sub([net.inp_img_size, net.inp_img_size], pre_dend(neuron_to_plot, :));
drawweightedlocs(rows, cols, variance_dend(neuron_to_plot, :), net.inp_img_size);
title(sprintf('Dendrite variance of N%d', neuron_to_plot));

% Axon connections
subplot(num_rows, num_cols, 1 * num_cols + 2);
[rows, cols] = ind2sub([net.inp_img_size, net.inp_img_size], post_axon(neuron_to_plot, :) - net.N_inp - net.N_hid);
drawlocs(rows, cols, net.inp_img_size);
title('Axon connections');

%% Axon weights
subplot(num_rows, num_cols, 1 * num_cols + 3);
[rows, cols] = ind2sub([net.inp_img_size, net.inp_img_size], post_axon(neuron_to_plot, :) - net.N_inp - net.N_hid);
drawweightedlocs(rows, cols, w_axon(neuron_to_plot, :), net.inp_img_size);
title(sprintf('Axon weights of N%d', neuron_to_plot));

%% Axon delays
subplot(num_rows, num_cols, 1 * num_cols + 4);
[rows, cols] = ind2sub([net.inp_img_size, net.inp_img_size], post_axon(neuron_to_plot, :) - net.N_inp - net.N_hid);
drawweightedlocs(rows, cols, delays_axon(neuron_to_plot, :), net.inp_img_size);
title(sprintf('Axon delays of N%d', neuron_to_plot));

% Axon variance
subplot(num_rows, num_cols, 1 * num_cols + 5);
[rows, cols] = ind2sub([net.inp_img_size, net.inp_img_size], post_axon(neuron_to_plot, :) - net.N_inp - net.N_hid);
drawweightedlocs(rows, cols, variance_axon(neuron_to_plot, :), net.inp_img_size);
title(sprintf('Axon variance of N%d', neuron_to_plot));

% Draw spikes
subplot(num_rows, 1, 3);
hold off
offset = (sec - 1) * 1000;
filter = find(spike_times_trace(:,2) <= net.N_inp + net.N_hid & spike_times_trace(:,2) > net.N_inp & spike_times_trace(:,1) > offset);
l2_spike_idxs = spike_times_trace(filter, 2);
l2_spike_times = spike_times_trace(filter, 1);
filter = find(spike_times_trace(:,2) < net.N_inp + net.N_hid & spike_times_trace(:, 1) > offset);
l1_idxs = spike_times_trace(filter, 2);
l1_times = spike_times_trace(filter, 1);
% Note this is the SPIKE TIME (not arrival time)
plot(l1_times - offset, l1_idxs, '.k', 'MarkerSize', 8);
hold on 
plot(l2_spike_times - offset, l2_spike_idxs, '.r', 'MarkerSize', 8)
ax = gca;
i = 0;
%axis([0, 30, -30 N_v + 50]);
while i < numel(l2_spike_times)
    i = i + 1;
    pos = l2_spike_times(i);
    c_idx = mod(l2_spike_idxs(i) - net.N_inp - 1, size(ax.ColorOrder, 1)) + 1;
    colour = ax.ColorOrder(c_idx, :);
    plot( [pos pos] - offset, get( gca, 'Ylim' ), '--', 'Color', colour, 'LineWidth',2);
end
xlabel('Time (ms)');
ylabel('Neuron number');

% Draw voltages
subplot(num_rows, 1, 4);
plot(vt(1:3, :)');
drawnow;


%% helper functions
function [] = drawweightedlocs( rows, cols, weightings, inp_img_size )

    pts = linspace(1, inp_img_size + 1, inp_img_size + 1) - 0.5;
    counts = histcounts2(rows(:), cols(:), pts, pts );
    weighted = accumarray([rows(:), cols(:)], weightings(:), [inp_img_size, inp_img_size]);
    %padded = zeros(inp_img_size, inp_img_size);
    %idxs = find(weighted);
    %padded(idxs(:)) = weighted(idxs(:));
    imagesc(weighted);
    axis equal;
    set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]));
    colorbar
    
end

function [] = drawlocs( rows, cols, inp_img_size )

    pts = linspace(1, inp_img_size + 1, inp_img_size + 1) - 0.5;
    counts = histcounts2(rows(:), cols(:), pts, pts );
    imagesc(counts);
    axis equal;
    set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]));
    colorbar

end
