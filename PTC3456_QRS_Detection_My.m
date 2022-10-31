%%% ECG signal QRS complex detection
close; clear; clc; name = '119m';
fid = fopen([name '.info'], 'rt');
fgetl(fid); fgetl(fid); fgetl(fid);
freqint = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
fs = freqint(1); fgetl(fid);
[row, signal, gain, base, unit] = strread(fgetl(fid),'%d%s%f%f%s','delimiter','\t');
fclose(fid);

%%% Input ECG signal
data = load([name '.mat']);
x = (data.val-base)/gain;
l = length(x);
t = (0:l-1) / fs;
mean_x = mean(x);
std_x = std(x);

%%% MA FIR low-pass filter
MA_fo = 10; % filter order
x_ma = filter(ones(1,MA_fo)/MA_fo, 1, x);

%%% Butterworth high-pass filter
B_fc = 5; B_fo = 4; % cutting frequency and filter order
[H_B_num, H_B_den] = butter(B_fo, 2*B_fc/fs, 'high');
x_ma_b = filter(H_B_num, H_B_den, x_ma);

%%% Hilbert transform
x_ma_b_h = abs(hilbert(x_ma_b));

%%% My detection algorithm
dc = 0; % detection counter
cr = 20; % comparison range + 1
det_v = zeros(1,l); % detected points vector
tshd = 0.35 * max(x_ma_b_h); % threshold
for i = 2:l-1
    if i < cr+1
        c1 = 1;
        c2 = i + cr - 1;
    elseif i < l-cr
        c1 = i - cr + 1;
        c2 = i + cr - 1;
    else
        c1 = i - cr + 1;
        c2 = l;
    end
    svb = max(x_ma_b_h(c1:i-1)); % before
    svn = x_ma_b_h(i); % now
    sva = max(x_ma_b_h(i+1:c2)); % after
    if (svn>svb) && (svn>sva) && (svn>tshd)
        dc = dc + 1;
        det_v(dc) = i;
    end
end
det_v = det_v(1:dc);

% RR interval
IRR = zeros(1,dc-1);
IRR(1) = det_v(1);
for i = 2:dc
    IRR(i) = det_v(i) - det_v(i-1);
end
IRR = IRR / fs;
mean_IRR = mean(IRR);
std_IRR = std(IRR);

%%% Comparison
x = x / norm(x, 'inf'); % x = x/max(x)
mean_xn = mean(x);
std_xn = std(x);
x_ma_b_h = x_ma_b_h / norm(x_ma_b_h, 'inf');

figure(1);
plot(t, x, t, x_ma_b_h);
hold on;
scatter(t(det_v), x_ma_b_h(det_v), '*', 'k');
hold off;
title([ 'File ' name ' \Rightarrow \mu_{input}=' num2str(mean_x)...
    ' | \sigma_{input}=' num2str(std_x) ' \Rightarrow \mu_{norm input}='...
    num2str(mean_xn) ' | \sigma_{norm input}=' num2str(std_xn)]);
xlabel('Time [s]');
ylabel('Infinite normalized amplitude','FontWeight','bold');
legend('Input','Envelope of filtered input','Peak detection',...
    'Location','southeast');
ylim([min(min(x),min(x_ma_b_h)), max(max(x),max(x_ma_b_h))]);
grid minor;

figure(2);
plot(t(det_v), IRR);
title(['File ' name ' \Rightarrow RR interval \Rightarrow \mu_{IRR}='...
    num2str(mean_IRR) ' | \sigma_{IRR}=' num2str(std_IRR)]);
xlabel('Time [s]');
ylabel('IRR [s]');
grid minor;