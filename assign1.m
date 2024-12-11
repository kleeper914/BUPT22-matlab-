T = 1;      % 时间序列1s
% 原始参数
% N = 10000, K = 50
% 若要进行最后一步y(t)的最佳采样，可将N, K更改为
% N = 100, K = 1e5
% 但是这个更改会使其他结果难以观察
% 仅用于验证最后一步，验证前面4个步骤时使用原始即可
N = 10000;   % 采样率
K = 50;     % 符号数
dt = 1/N;   % 采样周期
t = linspace(0, T, N); % 时间序列
%% (1)
% 生成不归零矩阵脉冲序列
rect_pulse = (t >= 0 & t <= T);

% 生成二进制符号序列
% 学号尾号为3
A_k = -3 * ones(1, K);  % 初始全为-3
prob = rand(1, K);      % 生成0-1的随机数
A_k(prob <= 0.2) = 3;   % 概率小于等于0.2设为+3

% 构造信号S(t)
S = zeros(1, length(t)*K);
for k=1:K
    % 每个符号B(k)乘以脉冲波形，并进行时移
    S((k-1)*length(t)+1 : k*length(t)) = A_k(k) * rect_pulse;  % 累加构成S(t)信号
end

% 估计自相关函数
[R, lags] = xcorr(S, 'coeff');

% 功率谱估计
nfft = length(S);
overlap = 2^11;
[PSD, f] = periodogram(S, hann(length(S)), overlap, N);

% 绘制信号 S(t)
figure;
x_1 = linspace(0, K*T, length(S));
plot(x_1, S);
xticks(0:k*T);
ylim([-4, 4])
title('样本波形 S(t)');
xlabel('时间 (秒)');
ylabel('幅度');

% 绘制自相关函数
figure;
plot(lags*dt, R);
title('自相关函数 R_S(t)');
xlabel('时间延迟 (秒)');
ylabel('自相关');

% 绘制功率谱
figure;
plot(f, 10*log10(PSD));
title('功率谱S(f)');
xlabel('频率 (Hz)');
ylabel('功率谱密度 (dB/Hz)');

%% (2)
% 生成正态随机噪声
% 学号倒数第二位为7
N_0 = 0.7;
% 方差为自相关函数在0的取值
n_t = sqrt(N_0/2) * randn(1, K*N);

% 设计带通滤波器
f_high = round(2 * pi);
[b, a] = butter(6, f_high/(N/2), 'low');

% 滤波
n_t_filtered = filter(b, a, n_t);

% 估计自相关函数
[R_1, lags_1] = xcorr(n_t, 'coeff');
[R_2, lags_2] = xcorr(n_t_filtered, 'coeff');

% 估计功率谱
[F_1, f_1] = periodogram(n_t, rectwin(length(n_t)), overlap, N);
[F_2, f_2] = periodogram(n_t_filtered, rectwin(length(n_t_filtered)), overlap, N);


% 绘制噪声信号
figure;
subplot(1, 2, 1);
plot(x_1, n_t);
title('滤波前噪声信号');
ylim([-10, 10]);
xlabel('时间(s)');
ylabel('幅度');
xticks(0:1);
subplot(1, 2, 2);
plot(x_1, n_t_filtered);
ylim([-1, 1]);
title('滤波后噪声信号');
xlabel('时间(s)');
ylabel('幅度');
xticks(0:1);
% 绘制自相关函数
figure;
subplot(1, 2, 1);
plot(lags_1*dt, R_1);
title('滤波前噪声信号自相关函数');
ylim([-0.5, 1]);
subplot(1, 2, 2);
plot(lags_2*dt, R_2);
title('滤波后噪声信号自相关函数');
ylim([-0.5, 1]);
% 绘制功率谱密度
figure;
subplot(1, 2, 1);
plot(f_1, 10*log10(F_1));
title('滤波前噪声信号功率谱密度');
xlabel('频率 (Hz)');
ylabel('功率谱密度dB/Hz');
xlim([0, 1000]);
subplot(1, 2, 2);
plot(f_2, 10*log10(F_2));
title('滤波后噪声信号功率谱密度');
xlabel('频率 (Hz)');
ylabel('功率谱密度dB/Hz');
xlim([0, 50]);

%% (3)
r = S + n_t_filtered;
[R_3, lags_3] = xcorr(r, 'coeff');
[F_3, f_3] = periodogram(r, hann(length(r)), overlap, N);

% 绘制r(t)的样本波形
figure;
plot(x_1, r);
title('r(t)');
ylim([-4, 4]);
xlabel('时间/s');
ylabel('幅度');
% 绘制r(t)的自相关函数
figure;
plot(lags_3/N, R_3);
title('r(t)自相关函数');
% 绘制r(t)的功率谱
figure;
plot(f_3, 10*log10(F_3));
xlabel('频率/Hz');
ylabel('dB/Hz');
title('r(t)的功率谱');

%% (4)
% 匹配滤波器与发传号即+3 * p(t)时匹配
% h1 = 3 * p(Tb - t)
h1 = 3 * rect_pulse;
% y(t) 等于 r(t)与h(t)的卷积
y = zeros(1, length(r));
for k=1:K
    y((k-1)*length(t)+1 : k*length(t)) = conv(r((k-1)*length(t)+1 : k*length(t)), h1, 'same') / N;
end
[R_4, lags_4] = xcorr(y, 'coeff');
[F_4, f_4] = periodogram(y, hann(length(y)), overlap, N);

% 绘制y(t)的样本波形
figure;
plot(x_1, y);
title('y(t)');
xlabel('时间/s');
ylabel('幅度');
%绘制y(t)的自相关函数
figure;
plot(lags_4 / N, R_4);
title('y(t)自相关函数');
xlabel('时间/s');
%绘制y(t)的功率谱密度
figure;
plot(f_4, 10*log10(F_4));
xlabel('频率/Hz');
ylabel('dB/Hz');
title('y(t)的功率谱');

%% (5)
M = 10^5;   % 采样点数
sample_idx = linspace(1, K*T, K); % 生成采样序列，在t=kTs时采样
y_sample = y(sample_idx);

% 画出概率密度函数
figure;
[f, yi] = ksdensity(y_sample);  % 归一化为概率密度
plot(yi, f);
title('y的概率密度函数');
xlabel('幅度');
ylabel('概率密度');
grid on;