%% VISUALISE1DWEIGHTS - Visualise network variables in the workspace

num_rows = 4;
num_cols = 5;

subplot(2, num_cols, 1);

idxs = ts > (sec - 1) * 1000 & ts <= (sec * 1000);
plot(inp(idxs), ts(idxs), '.k');

% Dendrite weights
subplot(num_rows, num_cols, 3);
imagesc(w_dend(neuron_to_plot, :)');
title('Weight dendrites');
colorbar;

subplot(num_rows, num_cols, 4);
imagesc(delays_dend(neuron_to_plot, :)');
title('Delays dendrites');
colorbar;

subplot(num_rows, num_cols, 5);
imagesc(variance_dend(neuron_to_plot, :)');
title('Variance dendrite');
colorbar;

subplot(num_rows, num_cols, 8);
imagesc(w_axon(neuron_to_plot, :)');
title('Weight axons');
colorbar;

subplot(num_rows, num_cols, 9);
imagesc(delays_axon(neuron_to_plot, :)');
title('Delays axons');
colorbar;

subplot(num_rows, num_cols, 10);
imagesc(variance_axon(neuron_to_plot, :)');
title('Varaince axons');
colorbar;

% Draw spikes
subplot(num_rows, 1, num_rows - 1);
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
axis([0 1000 1 net.N_v + 1]);
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
subplot(num_rows, 1, num_rows);
plot(vt(1, :)');
legend({'n1'});

drawnow;








