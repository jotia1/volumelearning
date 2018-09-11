%% CLUSTERGLOBALINTEXP

net = getwrightnet(1);

% variable last spike
%last_idxs = find(net.supplied_input == 3);
%net.supplied_ts(last_idxs) = net.supplied_ts(last_idxs) + randi([0, 4], 1, numel(last_idxs)) - 2;

% Supervise 
%net.supervised_seconds = 100;
n2_patts = find(mod(net.supplied_ts, 1000) == 8);
net.supplied_ts(n2_patts) = net.supplied_ts(n2_patts) + 3;
n2_sups = find(net.supplied_input == 4 & mod(net.supplied_ts, 1000) == 13);
net.supplied_input(n2_sups) = 5;


% Enable second neuron
net.w_dend(2, :) = net.w_dend(1, :);

% Get a valid folder to dump results to
count = 0;
exp_code = 'vgi'; % Variable gloabl integral
output_folder = [exp_code, '-comp00'];
while exist(output_folder, 'dir') == 7
    count = count + 1;
    output_folder = sprintf('%s-out%02d', exp_code, count);
end

mkdir(output_folder);
filename_format = 'fgi-%.02f';

for fgi = 4:0.2:19
    net.fgi = fgi;

    out = runsinglelayer(net);
    disp((fgi - 4.0) / (19 - 4));
    removed_dots = strrep(sprintf(['%s/', filename_format], output_folder, fgi), '.', '-');
    save([removed_dots, '.mat'], 'net', 'out');
end

