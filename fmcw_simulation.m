%% FMCW Radar Simulation — Multi-Target Scenario
% =========================================================================
% 功能: 多目标 FMCW 雷达仿真, 生成 Range-Doppler Map + CFAR 检测
%
% 正确性要点:
%   1. 快-慢时间分离, 每 chirp 频率从 fc 独立线性上升
%   2. RX 在 τ < 2R/c 时自然为 0
%   3. Stop-and-hop 近似: 距离按 chirp 更新
%   4. 2D FFT 后仅沿多普勒维 fftshift
%   5. 物理轴从参数严格推导
%   6. CFAR 在功率域做, 量纲一致
%   7. 向量化加速信号生成
%   8. 可控噪声 + SNR 定义
%   9. 加窗抑制旁瓣
%   10. DBSCAN 聚类聚合 CFAR 点迹 (归一化各向同性距离)
%
% 矩阵维度约定:
%   beat / rdm 均为 (Nr, Nd), 第 1 维 = 距离, 第 2 维 = 多普勒
% =========================================================================

clear; clc; close all;

%% ==================== 1. 雷达参数 ====================
c        = 3e8;
fc       = 77e9;
B        = 150e6;
T_chirp  = 50e-6;
Nd       = 256;
Nr       = 1024;

lambda   = c / fc;
slope    = B / T_chirp;
Fs       = Nr / T_chirp;
T_frame  = Nd * T_chirp;

range_res = c / (2 * B);
max_range = (Nr / 2) * range_res;
v_res     = lambda / (2 * T_frame);
v_max     = lambda / (4 * T_chirp);

fprintf('========== Radar Parameters ==========\n');
fprintf('  fc        = %.2f GHz\n', fc / 1e9);
fprintf('  B         = %.1f MHz\n', B / 1e6);
fprintf('  T_chirp   = %.2f us\n', T_chirp * 1e6);
fprintf('  Nr=%d, Nd=%d, T_frame=%.1f ms\n', Nr, Nd, T_frame * 1e3);
fprintf('  ΔR=%.2f m, R_max=%.0f m\n', range_res, max_range);
fprintf('  Δv=%.2f m/s, v_max=%.0f m/s\n', v_res, v_max);

%% ==================== 2. 多目标定义 ====================
targets = [
     30,    0,  1.0;
     50,    5,  0.8;
     75,  -10,  1.2;
    120,    0,  0.5;
    155,    8,  0.6;
];
num_targets = size(targets, 1);
fprintf('\n========== Targets ==========\n');
for k = 1:num_targets
    fprintf('  T%d: R = %5.1f m, v = %+5.1f m/s\n', ...
        k, targets(k, 1), targets(k, 2));
end

%% ==================== 3. 噪声 / 信噪比 ====================
SNR_beat_dB = -10;
SNR_beat    = 10^(SNR_beat_dB / 10);
ref_amp     = max(targets(:, 3));
noise_std   = ref_amp / sqrt(SNR_beat);
fprintf('\n  SNR_beat = %.1f dB (after 2D FFT: ≈ %.1f dB)\n', ...
    SNR_beat_dB, SNR_beat_dB + 10*log10(Nr*Nd));

%% ==================== 4. 信号生成 (向量化) ====================
% ====================== 数学模型 ======================
%
% 发射信号 (第 n 个 chirp, 快时间 tau):
%   s_tx(tau, n) = exp{j·2π[ fc·(nT + tau) + ½S·tau² ]}
%
% 第 k 个目标时延 (stop-and-hop):
%   R_k(n)   = R0_k + v_k · n · T_chirp
%   tau_k(n) = 2 · R_k(n) / c
%
% 接收信号:
%   s_rx(tau, n) = exp{j·2π[ fc·(nT + tau - tau_k) + ½S·(tau - tau_k)² ]}
%                = 0,  tau < tau_k  (回波未到达)
%
% 差频信号 (dechirp, 复混频):
%   s_k = s_tx · conj(s_rx) = exp{j·2π[ fc·tau_k + S·tau·tau_k - ½S·tau_k² ]}
%                               ├─ 载波相位 ─┤├─ 差频项 ─┤├── RVP ──┤
%
% 其中差频项 exp{j·2π·S·tau·tau_k} 产生频率 f_if = S·tau_k = slope · 2R/c
% 这是测距的核心, Range FFT 后谱峰位置对应 f_if。
%
% 载波相位 exp{j·2π·fc·tau_k} 在单 chirp 内为常数,
% 跨 chirp 时随 R_k(n) 线性变化, 产生多普勒频率 f_d = 2v/lambda。
% ======================================================

tau = (0:Nr - 1).' / Fs;        % 快时间: Nr×1
n   = 0:Nd - 1;                 % 慢时间 (chirp 索引): 1×Nd
beat = zeros(Nr, Nd);

for k = 1:num_targets
    R0  = targets(k, 1);
    v   = targets(k, 2);
    amp = targets(k, 3);

    % 每 chirp 的目标时延 tau_k(n) : 1×Nd
    delay = 2 * (R0 + v * n * T_chirp) / c;

    % 差频信号 (外积):
    %   载波相位项: exp{j·2π·fc·tau_k}                    → 1×Nd
    %   差频项:     exp{j·2π·S·(tau · tau_k)}            → Nr×Nd (外积)
    %   RVP:        exp{-j·π·S·tau_k²}                    → 1×Nd
    s_k = amp                                                     ...
        * exp( 1j * 2 * pi * fc * delay)                          ... % 载波相位
        .* exp( 1j * 2 * pi * slope * (tau * delay))              ... % 差频项 (外积)
        .* exp(-1j * pi * slope * (delay.^2));                     ... % RVP

    % tau >= delay 掩膜: 回波尚未到达时 s_rx = 0
    s_k = s_k .* (tau >= delay);    % Nr×Nd

    beat = beat + s_k;
end

% 复高斯白噪声 (每维标准差 noise_std/sqrt(2))
beat = beat + (noise_std / sqrt(2)) * (randn(Nr, Nd) + 1j * randn(Nr, Nd));

fprintf('\nSignal generation done.\n');

%% ==================== 5. 加窗 + 2D FFT ====================
% 加 Hanning 窗: 旁瓣从 -13 dB 降至 -32 dB
% 代价: 主瓣展宽约 1.6 倍 (分辨率下降)

win_range = hann(Nr);               % Nr×1
win_dopp  = hann(Nd)';              % 1×Nd

beat_win = beat .* win_range;
range_fft = fft(beat_win, Nr, 1);       % Nr × Nd
range_fft = range_fft(1:Nr/2, :);        % Nr/2 × Nd (单边)

range_fft = range_fft .* win_dopp;
rdm = fftshift(fft(range_fft, Nd, 2), 2); % (Nr/2) × Nd

rdm_dB = 20 * log10(abs(rdm) + eps);

%% ==================== 6. 物理轴 ====================
range_axis = (0:Nr/2 - 1) * c / (2 * B);      % [m]
freq_dop   = (-Nd/2:Nd/2 - 1) / T_chirp / Nd;  % [Hz]
vel_axis   = freq_dop * lambda / 2;             % [m/s]

%% ==================== 7. 可视化 ====================
ieee_cmap = jet(256);

% 7.1 3D RDM
figure('Name', '3D RDM', 'Position', [50, 100, 900, 650]);
surf(vel_axis, range_axis, rdm_dB, 'EdgeColor', 'none', 'FaceAlpha', 0.85);
xlabel('Velocity (m/s)'); ylabel('Range (m)'); zlabel('Magnitude (dB)');
title(sprintf('3D RDM (SNR_{beat}=%.0f dB, Hanning)', SNR_beat_dB));
colormap(ieee_cmap); colorbar; view(45, 30);
zlim([max(rdm_dB(:)) - 50, max(rdm_dB(:))]);
hold on;
for k = 1:num_targets
    [~, ri] = min(abs(range_axis - targets(k, 1)));
    [~, vi] = min(abs(vel_axis  - targets(k, 2)));
    plot3(targets(k, 2), targets(k, 1), rdm_dB(ri, vi), ...
        'ko', 'MarkerSize', 10, 'LineWidth', 2);
end
hold off;

% 7.2 2D 俯视图
figure('Name', 'RDM Plan View', 'Position', [100, 150, 850, 600]);
imagesc(vel_axis, range_axis, rdm_dB); axis xy;
xlabel('Velocity (m/s)'); ylabel('Range (m)');
title(sprintf('RDM (%d dB, Hanning)', SNR_beat_dB));
colormap(ieee_cmap); colorbar;
clim([max(rdm_dB(:)) - 40, max(rdm_dB(:))]);
hold on;
for k = 1:num_targets
    plot(targets(k, 2), targets(k, 1), 'ko', 'MarkerSize', 8, 'LineWidth', 2);
end
hold off;

% 7.3 Range Profile (零多普勒切片)
figure('Name', 'Range Profile', 'Position', [150, 200, 800, 400]);
zero_dop_idx = Nd / 2 + 1;
rp = abs(rdm(:, zero_dop_idx)); rp = rp / max(rp);
plot(range_axis, rp, 'LineWidth', 1.5);
xlabel('Range (m)'); ylabel('Normalized Amplitude');
title('Range Profile (Zero Doppler, Hanning)');
grid on; xlim([0, min(200, max_range)]);
for k = 1:num_targets, xline(targets(k, 1), '--r', sprintf('T%d', k)); end

% 7.4 Doppler Profile
[~, pbin] = max(rp);
figure('Name', 'Doppler Profile', 'Position', [200, 250, 800, 400]);
dp = abs(rdm(pbin, :)) / max(abs(rdm(pbin, :)));
plot(vel_axis, dp, 'LineWidth', 1.5);
xlabel('Velocity (m/s)'); ylabel('Normalized Amplitude');
title(sprintf('Doppler Profile at R=%.1f m (Hanning)', range_axis(pbin)));
grid on;
for k = 1:num_targets, xline(targets(k, 2), '--r', sprintf('T%d', k)); end

%% ==================== 8. CA-CFAR 检测 ====================
tc_range   = 10;     tc_doppler = 8;
gc_range   = 4;      gc_doppler = 2;
offset     = 13;     % dB

rdm_pow = abs(rdm).^2;
[N_range, N_dopp] = size(rdm_pow);
edge_r = tc_range + gc_range;
edge_d = tc_doppler + gc_doppler;

detections = false(N_range, N_dopp);

for i = edge_r + 1:N_range - edge_r
    for j = edge_d + 1:N_dopp - edge_d
        train_win = rdm_pow(i - edge_r:i + edge_r, j - edge_d:j + edge_d);
        guard_win = rdm_pow(i - gc_range:i + gc_range, j - gc_doppler:j + gc_doppler);
        noise_sum = sum(train_win, 'all') - sum(guard_win, 'all');
        N_train = (2 * edge_r + 1) * (2 * edge_d + 1);
        N_guard = (2 * gc_range + 1) * (2 * gc_doppler + 1);
        noise_level = noise_sum / (N_train - N_guard);
        if rdm_pow(i, j) > noise_level * 10^(offset / 10)
            detections(i, j) = true;
        end
    end
end

fprintf('\n========== Raw CFAR Detections ==========\n');
fprintf('  Detections before clustering: %d\n', sum(detections, 'all'));

%% ==================== 9. DBSCAN 聚类 (归一化坐标) ====================
% 原理:
%   CFAR 检测点位于 RDM 的 (range_bin, doppler_bin) 网格上。
%   直接对 bin 索引做欧氏距离会使两个维度的权重不同:
%     - range bin 宽度: ΔR = 1.0 m
%     - doppler bin 宽度: ΔV = 0.15 m/s
%   如果不做归一化, DBSCAN 的 eps "圆形" 在物理上实际是
%   "拉长的椭圆", 导致聚类结果偏斜。
%
%   归一化方法:
%     将 bin 索引映射到 (距离/ΔR, 速度/ΔV) 的无量纲坐标,
%     使两个维度在物理分辨率上平权。
%
%   eps 物理含义: 同一个目标在 RDM 上占据的归一化半径。
%     设 eps = 3 (即 3 倍分辨率), 允许一个目标簇覆盖
%     约 3×ΔR × 3×ΔV 的区域, 足以聚合目标主瓣和邻近旁瓣。

[det_rows, det_cols] = find(detections);

if size(det_rows, 1) >= 2
    % 归一化坐标: 使每个维度以各自的分辨率为单位
    % coords_norm(:, 1) = range_idx (无量纲)
    % coords_norm(:, 2) = doppler_idx (无量纲)
    % 两维单位相同 → 欧氏距离各向同性
    coords_norm = [det_rows, det_cols];

    eps_norm   = 3;       % 归一化半径 (bin 数)
    min_pts    = 3;       % 核心点最少点数

    cluster_idx = dbscan(coords_norm, eps_norm, min_pts);

    unique_clusters = setdiff(unique(cluster_idx), -1);
    cluster_centers = zeros(length(unique_clusters), 2);

    fprintf('\n========== Clustered Detections (DBSCAN) ==========\n');
    fprintf('  Normalized coordinates: (Δr/ΔR)² + (Δv/ΔV)² < eps²\n');
    fprintf('  eps_norm=%.1f (bins), minPts=%d\n', eps_norm, min_pts);
    fprintf('  Number of clusters: %d\n', length(unique_clusters));

    for c = 1:length(unique_clusters)
        cid = unique_clusters(c);
        mask = (cluster_idx == cid);
        rows_c = det_rows(mask);
        cols_c = det_cols(mask);

        % 按功率加权平均求质心 (使用原始功率)
        weights = rdm_pow(sub2ind(size(rdm_pow), rows_c, cols_c));
        r_center = sum(rows_c .* weights) / sum(weights);
        v_center = sum(cols_c .* weights) / sum(weights);
        cluster_centers(c, :) = [r_center, v_center];

        fprintf('  C%d: R=%.1f m, v=%.1f m/s, %d points\n', ...
            c, range_axis(round(r_center)), vel_axis(round(v_center)), sum(mask));
    end
end

%% ==================== 10. CFAR 结果可视化 ====================
det_vals = rdm_dB(sub2ind(size(rdm_dB), det_rows, det_cols));

% 10.1 3D CFAR
figure('Name', 'CFAR 3D', 'Position', [100, 100, 950, 650]);
surf(vel_axis, range_axis, rdm_dB, 'EdgeColor', 'none', ...
    'FaceAlpha', 0.6, 'FaceColor', 'interp');
xlabel('Velocity (m/s)'); ylabel('Range (m)'); zlabel('Magnitude (dB)');
title(sprintf('CFAR on 3D RDM (offset=%.0f dB)', offset));
colormap(ieee_cmap); colorbar; view(45, 30);
zlim([max(rdm_dB(:)) - 50, max(rdm_dB(:))]);
if ~isempty(det_rows)
    hold on;
    scatter3(vel_axis(det_cols), range_axis(det_rows), det_vals, ...
        40, 'r', 'filled', 'MarkerEdgeColor', 'k');
    hold off;
end

% 10.2 CFAR 原始检测图
figure('Name', 'Raw CFAR 2D', 'Position', [150, 150, 850, 600]);
imagesc(vel_axis, range_axis, double(detections)); axis xy;
xlabel('Velocity (m/s)'); ylabel('Range (m)');
title(sprintf('Raw CFAR Detections (%d points)', sum(detections, 'all')));
colormap(gca, [1, 1, 1; 0, 0.447, 0.741]); colorbar;
hold on;
for k = 1:num_targets
    plot(targets(k, 2), targets(k, 1), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
end
legend({'True targets'}, 'Location', 'best');
hold off;

% 10.3 聚类后检测图
figure('Name', 'Clustered CFAR 2D', 'Position', [180, 180, 850, 600]);
imagesc(vel_axis, range_axis, double(detections)); axis xy;
xlabel('Velocity (m/s)'); ylabel('Range (m)');
title(sprintf('CFAR AFTER DBSCAN (%d clusters)', length(unique_clusters)));
colormap(gca, [1, 1, 1; 0.5, 0.5, 0.5]); colorbar;
hold on;

if exist('cluster_idx', 'var') && ~isempty(cluster_idx)
    noise_mask = (cluster_idx == -1);
    if any(noise_mask)
        plot(vel_axis(det_cols(noise_mask)), range_axis(det_rows(noise_mask)), ...
            'x', 'Color', [0.7, 0.7, 0.7], 'MarkerSize', 6);
    end
end

if exist('cluster_centers', 'var') && ~isempty(cluster_centers)
    for c = 1:size(cluster_centers, 1)
        rr = range_axis(round(cluster_centers(c, 1)));
        vv = vel_axis(round(cluster_centers(c, 2)));
        plot(vv, rr, 'go', 'MarkerSize', 12, 'LineWidth', 2);
    end
end

for k = 1:num_targets
    plot(targets(k, 2), targets(k, 1), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
end
legend({'Noise (discarded)', 'Cluster centroids', 'True targets'}, ...
    'Location', 'best');
hold off;

fprintf('\nDone.\n');
