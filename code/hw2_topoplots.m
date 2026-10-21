%% hw2_topoplots.m  -  HW2 task 2.2.3: topoplots of channel-summed Fisher scores
% For each subject and session (off, on-1..on-6): D(ch) = sum over the 14 frequencies (4:2:30 Hz) of the Fisher score.
% Input : hw2_channel_sumF.csv (variant, subject, session, channel, chLabel, sumF), ch32Locations.mat
% Output: hw2_figs/topo_*.png, hw2_topo_summary.csv, hw2_topo_console.txt
% Colour scale is shared across sessions within each subject so sessions are comparable.
% Self-contained: uses the local drawTopo function (theta/radius from ch32Locations + griddata 'v4'),
% so no EEGLAB, topoplot.m, readlocs.m, or finputcheck.m is needed.
clear; clc; close all;
if exist('hw2_topo_console.txt','file'), delete('hw2_topo_console.txt'); end
diary('hw2_topo_console.txt'); diary on; cleanupObj = onCleanup(@() diary('off'));
fprintf('hw2_topoplots run at %s\n', datestr(now));

%% ---------- config ----------
EXCLUDE_MASTOIDS = true;   % true: maps drawn without M1/M2 (radius 0.747 shrinks the head); summary table still uses all 32
variants = {'A_win_linear','B_win_log10','C_trial_mean_linear'};
sessLab  = {'off','on-1','on-2','on-3','on-4','on-5','on-6'};
figDir = 'hw2_figs'; if ~exist(figDir,'dir'), mkdir(figDir); end

%% ---------- load ----------
T = readtable('hw2_channel_sumF.csv');
L = load('ch32Locations.mat'); chanlocs = L.ch32Locations; locLabels = {chanlocs.labels}';
assert(numel(chanlocs) == 32, 'Expected 32 channel locations');
fprintf('Loc radius: M1 = %.3f, M2 = %.3f, max over channels = %.3f\n', ...
    chanlocs(strcmp(locLabels,'M1')).radius, chanlocs(strcmp(locLabels,'M2')).radius, max([chanlocs.radius]));

if EXCLUDE_MASTOIDS, plotIdx = find(~ismember(locLabels, {'M1','M2'}))'; else, plotIdx = 1:32; end
thC3 = chanlocs(strcmp(locLabels,'C3')).theta; thC4 = chanlocs(strcmp(locLabels,'C4')).theta;
assert(sin(deg2rad(thC4)) > 0 && sin(deg2rad(thC3)) < 0, ...
    'Left/right orientation is mirrored; flip the sign of x in drawTopo');
fprintf('Orientation check passed: C3 theta = %g, C4 theta = %g\n', thC3, thC4);

% channel groups for the summary table (my grouping, a proposal, not from the HW)
G.motor   = {'C3','Cz','C4','FC1','FC2','CP1','CP2'};
G.latFCCP = {'FC5','FC6','CP5','CP6'};
G.periph  = {'Fp1','Fpz','Fp2','F7','F8','T7','T8','M1','M2'};
G.occ     = {'O1','Oz','O2','POz'};
allGrouped = [G.motor G.latFCCP G.periph G.occ];
G.rest = setdiff(locLabels', allGrouped, 'stable');
assert(numel(allGrouped) + numel(G.rest) == 32, 'Channel groups must cover all 32 channels once');

%% ---------- figures and summary ----------
sumRows = {};
for v = 1:numel(variants)
    V = nan(32, 7, 2);                                   % channel x session(off..on-6) x subject
    for si = 1:2
        for s = 0:6
            rows = strcmp(string(T.variant), variants{v}) & T.subject == si & T.session == s;
            d = T(rows, :); [~, o] = sort(d.channel); d = d(o, :);
            assert(height(d) == 32, 'Missing channel rows for %s S%d s%d', variants{v}, si, s);
            assert(all(strcmpi(string(d.chLabel), string(locLabels))), 'Channel order mismatch with ch32Locations');
            V(:, s+1, si) = d.sumF;
        end
    end

    modes = {'shared'};
    if v == 1, modes = {'shared','normalized'}; end
    for mi = 1:numel(modes)
        W = V;
        if strcmp(modes{mi}, 'normalized')
            for si = 1:2, for s = 1:7, W(:, s, si) = V(:, s, si) / sum(V(:, s, si)); end, end
            cbLabel = 'Fraction of session total';
        else
            cbLabel = 'Summed Fisher score';
        end
        fig = figure('Position', [30 30 1900 600]);
        tiledlayout(2, 7, 'TileSpacing', 'compact', 'Padding', 'compact');
        for si = 1:2
            cmax = max(W(:, :, si), [], 'all');
            fprintf('%s | %s | Subject %d colour limit [0, %.3f]\n', variants{v}, modes{mi}, si, cmax);
            for s = 1:7
                ax = nexttile;
                drawTopo(ax, W(:, s, si), chanlocs, [0 cmax], plotIdx);
                title(ax, sprintf('Subj %d | %s', si, sessLab{s}));
                if s == 7, cb = colorbar(ax); cb.Label.String = cbLabel; end
            end
        end
        sgtitle(sprintf('2.2.3  Channel-summed Fisher scores, 4-30 Hz (%s, %s; one colour scale per subject)', ...
            variants{v}, modes{mi}), 'Interpreter', 'none');
        saveFig(fig, fullfile(figDir, sprintf('topo_%s_%s.png', variants{v}, modes{mi})));
    end

    for si = 1:2
        for s = 0:6
            x = V(:, s+1, si); tot = sum(x);
            [xs, ord] = sort(x, 'descend');
            fr = @(names) sum(x(ismember(locLabels, names))) / tot;
            sumRows(end+1, :) = {variants{v}, si, s, tot, locLabels{ord(1)}, xs(1), locLabels{ord(2)}, xs(2), ...
                locLabels{ord(3)}, xs(3), fr(G.motor), fr(G.latFCCP), fr(G.periph), fr(G.occ), fr(G.rest)}; %#ok<SAGROW>
        end
    end
end

S = cell2table(sumRows, 'VariableNames', {'variant','subject','session','totalF','top1','top1F','top2','top2F', ...
    'top3','top3F','fracMotorCentral','fracLateralFCCP','fracPeripheral','fracOccipital','fracOtherFPz'});
writetable(S, 'hw2_topo_summary.csv');
fprintf('\nVariant A summary:\n'); disp(S(strcmp(S.variant, 'A_win_linear'), :));
fprintf('Groups: motor=%s | lateralFCCP=%s | peripheral=%s | occipital=%s | other=%s\n', ...
    strjoin(G.motor,','), strjoin(G.latFCCP,','), strjoin(G.periph,','), strjoin(G.occ,','), strjoin(G.rest,','));
fprintf('\nSaved figures in %s and hw2_topo_summary.csv\n', figDir);
diary off;

%% ---------- local functions ----------
function drawTopo(ax, vals, chanlocs, lim, useIdx)
    th = deg2rad([chanlocs(useIdx).theta]'); rd = [chanlocs(useIdx).radius]';
    x = rd .* sin(th); y = rd .* cos(th);                      % nose up, positive theta = right
    sq = 0.5 / (max(rd) * 1.02); x = x * sq; y = y * sq;       % head circle has radius 0.5
    v = vals(useIdx); v = v(:);
    [xi, yi] = meshgrid(linspace(-0.5, 0.5, 120));
    zi = griddata(x, y, v, xi, yi, 'v4'); zi(hypot(xi, yi) > 0.5) = NaN;
    hold(ax, 'on');
    surface(ax, xi, yi, zeros(size(zi)), zi, 'EdgeColor', 'none', 'FaceColor', 'interp');
    t = linspace(0, 2*pi, 200); plot(ax, 0.5*cos(t), 0.5*sin(t), 'k', 'LineWidth', 1.8);
    plot(ax, [-0.08 0 0.08], [0.495 0.58 0.495], 'k', 'LineWidth', 1.8);                            % nose
    plot(ax, [-0.5 -0.54 -0.55 -0.54 -0.5], [0.09 0.07 -0.02 -0.10 -0.11], 'k', 'LineWidth', 1.8); % left ear
    plot(ax, [0.5 0.54 0.55 0.54 0.5], [0.09 0.07 -0.02 -0.10 -0.11], 'k', 'LineWidth', 1.8);      % right ear
    scatter(ax, x, y, 10, 'k', 'filled');
    ax.CLim = lim; colormap(ax, jet);
    axis(ax, 'equal'); axis(ax, 'off'); xlim(ax, [-0.65 0.65]); ylim(ax, [-0.55 0.65]);
end

function saveFig(fig, path)
    try, exportgraphics(fig, path, 'Resolution', 200);
    catch, saveas(fig, path); end
end