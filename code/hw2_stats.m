%% hw2_stats.m  -  HW2 Section 2.1 (BCI command delivery): tasks 2.1.1, 2.1.2, 2.1.3
% Input : hw2_run_metrics.csv (36 rows = 2 subjects x 6 online sessions x 3 runs)
% Output: hw2_figs/*.png, hw2_session_means.csv, hw2_stats_trend.csv, hw2_stats_vs_session1.csv,
%         hw2_stats_console.txt
% Units : data points are RUNS (slide 9). accuracy = hits/(hits+misses); pctTimeout = timeouts/trials (0-1).
clear; clc; close all; rng(1);
if exist('hw2_stats_console.txt','file'), delete('hw2_stats_console.txt'); end
diary('hw2_stats_console.txt'); diary on; cleanupObj = onCleanup(@() diary('off'));
fprintf('hw2_stats run at %s\n', datestr(now));

ALPHA = 0.05;
metrics     = {'accuracy','pctTimeout'};
metricLabel = {'Command delivery accuracy','Fraction of timeouts'};
scopeName   = {'Subject 1','Subject 2','Group'};
scopeSubj   = {1, 2, [1 2]};
cols = [0.15 0.35 0.75; 0.85 0.15 0.15];
figDir = 'hw2_figs'; if ~exist(figDir,'dir'), mkdir(figDir); end

%% ---------- load and validate ----------
M = readtable('hw2_run_metrics.csv');
M.subject = double(M.subject); M.session = double(M.session); M.run = double(M.run);
assert(height(M) == 36, 'Expected 36 run rows, found %d', height(M));
assert(isequal(sort(unique(M.session))', 1:6), 'Sessions must be 1..6');
gc = groupcounts(M, {'subject','session','run'}); assert(all(gc.GroupCount == 1), 'Duplicate run rows');
fprintf('Loaded %d run rows | NaN accuracy: %d | NaN timeouts: %d\n', height(M), sum(isnan(M.accuracy)), sum(isnan(M.pctTimeout)));

%% ---------- 2.1.1 run-level values and session means ----------
fprintf('\n===== 2.1.1 Run-level accuracy and timeouts =====\n');
disp(M(:, {'subject','session','run','hits','misses','timeouts','trials','accuracy','pctTimeout'}));
Smean = groupsummary(M, {'subject','session'}, {'mean','std'}, {'accuracy','pctTimeout'});
disp(Smean); writetable(Smean, 'hw2_session_means.csv');

f1 = figure('Position',[100 100 720 480]);
plotSessionLines(M, 'accuracy', 'Command delivery accuracy', cols, true);
title('2.1.1  Accuracy over online sessions (mean of 3 runs, bars = SD, dots = runs)');
saveFig(f1, fullfile(figDir,'fig1_accuracy_by_session.png'));
f2 = figure('Position',[100 100 720 480]);
plotSessionLines(M, 'pctTimeout', 'Fraction of timeouts', cols, false);
title('Timeouts over online sessions (mean of 3 runs, bars = SD, dots = runs)');
saveFig(f2, fullfile(figDir,'fig2_timeouts_by_session.png'));

%% ---------- 2.1.2 trend over sessions ----------
fprintf('\n===== 2.1.2 Trend over sessions (x = session 1..6, runs as data points) =====\n');
trend = struct();
for mi = 1:2
    for sc = 1:3
        d = M(ismember(M.subject, scopeSubj{sc}), :);
        x = d.session; y = d.(metrics{mi}); ok = ~isnan(y); x = x(ok); y = y(ok);
        [r, p]    = corr(x, y);
        [rho, ps] = corr(x, y, 'Type', 'Spearman');
        mdl = fitlm(x, y); ci = coefCI(mdl);
        t = struct('n',numel(x),'r',r,'p',p,'rho',rho,'ps',ps, ...
                   'slope',mdl.Coefficients.Estimate(2),'slopeLo',ci(2,1),'slopeHi',ci(2,2), ...
                   'slopeAdj',NaN,'pAdj',NaN);
        if sc == 3
            tb = table(x, categorical(d.subject(ok)), y, 'VariableNames', {'session','subj','y'});
            m2 = fitlm(tb, 'y ~ session + subj');
            t.slopeAdj = m2.Coefficients{'session','Estimate'};
            t.pAdj     = m2.Coefficients{'session','pValue'};
        end
        trend(mi,sc).t = t;
        fprintf('%-12s | %-9s | n=%2d | Pearson r=%6.3f p=%.4g | Spearman rho=%6.3f p=%.4g | slope=%7.4f [%7.4f, %7.4f]', ...
            metrics{mi}, scopeName{sc}, t.n, t.r, t.p, t.rho, t.ps, t.slope, t.slopeLo, t.slopeHi);
        if sc == 3, fprintf(' | subject-adjusted slope=%7.4f p=%.4g', t.slopeAdj, t.pAdj); end
        fprintf('\n');
    end
end
fprintf('Significant trend (Pearson p < %.2f):\n', ALPHA);
for mi = 1:2
    for sc = 1:3
        fprintf('  %-12s %-9s : %s\n', metrics{mi}, scopeName{sc}, yesno(trend(mi,sc).t.p < ALPHA));
    end
end

%% ---------- 2.1.3 each session vs online session 1 ----------
fprintf('\n===== 2.1.3 Session k (2..6) vs session 1: Wilcoxon rank-sum and Welch t-test, Holm over 5 comparisons =====\n');
vs1 = struct(); rowsVs = {};
for mi = 1:2
    for sc = 1:3
        d = M(ismember(M.subject, scopeSubj{sc}), :);
        base = d.(metrics{mi})(d.session == 1); base = base(~isnan(base));
        pr = nan(1,5); pt = nan(1,5); dm = nan(1,5); mk = nan(1,5); nk = nan(1,5);
        for k = 2:6
            v = d.(metrics{mi})(d.session == k); v = v(~isnan(v));
            pr(k-1) = ranksum(base, v);
            [~, pt(k-1)] = ttest2(base, v, 'Vartype', 'unequal');
            dm(k-1) = mean(v) - mean(base); mk(k-1) = mean(v); nk(k-1) = numel(v);
        end
        prH = holm(pr); ptH = holm(pt);
        vs1(mi,sc).pr = pr; vs1(mi,sc).prH = prH; vs1(mi,sc).pt = pt; vs1(mi,sc).ptH = ptH; vs1(mi,sc).dm = dm;
        minP = 2 / nchoosek(numel(base) + round(mean(nk)), numel(base));
        fprintf('\n-- %s | %s | baseline n=%d mean=%.3f | smallest attainable exact rank-sum p ~ %.3g --\n', ...
            metrics{mi}, scopeName{sc}, numel(base), mean(base), minP);
        Tview = table((2:6)', nk', mk', dm', pr', prH', pt', ptH', 'VariableNames', ...
            {'session','n','mean','diffVsS1','p_ranksum','p_ranksum_holm','p_ttest','p_ttest_holm'});
        disp(Tview);
        for k = 1:5
            rowsVs(end+1,:) = {metrics{mi}, scopeName{sc}, k+1, numel(base), nk(k), mean(base), mk(k), dm(k), ...
                pr(k), prH(k), pt(k), ptH(k)}; %#ok<SAGROW>
        end
    end
end

fprintf('\n===== 2.1.3 First significant session and whether it is sustained =====\n');
crit = {'ranksum raw','prr'; 'ranksum Holm','prH'; 'ttest raw','ptr'; 'ttest Holm','ptH'};
for mi = 1:2
    for sc = 1:3
        pv = {vs1(mi,sc).pr, vs1(mi,sc).prH, vs1(mi,sc).pt, vs1(mi,sc).ptH};
        for c = 1:4
            [first, sus] = firstSustained(pv{c}, ALPHA);
            if isempty(first)
                fprintf('%-12s %-9s %-13s : no session reaches p < %.2f\n', metrics{mi}, scopeName{sc}, crit{c,1}, ALPHA);
            else
                fprintf('%-12s %-9s %-13s : first significant = session %d | sustained through session 6: %s\n', ...
                    metrics{mi}, scopeName{sc}, crit{c,1}, first+1, yesno(sus));
            end
        end
    end
end

%% ---------- figures for 2.1.2 and 2.1.3 ----------
starsG = cell(2,6);
for mi = 1:2
    for k = 2:6, starsG{mi,k} = starString(vs1(mi,3).prH(k-1)); end
end
f3 = figure('Position',[100 100 1250 500]);
for mi = 1:2
    subplot(1,2,mi);
    drawPanel(M, metrics{mi}, metricLabel{mi}, starsG(mi,:), trend(mi,3).t.r, trend(mi,3).t.p);
end
sgtitle('Session means (bars, SD whiskers, run dots), group linear trend, * = group rank-sum vs on-1 (Holm)');
saveFig(f3, fullfile(figDir,'fig3_session_bars_trend_and_tests.png'));

f4 = figure('Position',[100 100 1250 500]);
for mi = 1:2
    subplot(1,2,mi);
    plotTrend(M, metrics{mi}, metricLabel{mi}, cols, trend, mi);
end
sgtitle('2.1.2  Run-level trend over online sessions');
saveFig(f4, fullfile(figDir,'fig4_trend_scatter.png'));

%% ---------- export tables ----------
rowsTr = {};
for mi = 1:2
    for sc = 1:3
        t = trend(mi,sc).t;
        rowsTr(end+1,:) = {metrics{mi}, scopeName{sc}, t.n, t.r, t.p, t.rho, t.ps, t.slope, t.slopeLo, t.slopeHi, t.slopeAdj, t.pAdj}; %#ok<SAGROW>
    end
end
writetable(cell2table(rowsTr, 'VariableNames', {'metric','scope','n','pearson_r','pearson_p','spearman_rho', ...
    'spearman_p','slope','slope_ci_lo','slope_ci_hi','slope_subject_adj','p_subject_adj'}), 'hw2_stats_trend.csv');
writetable(cell2table(rowsVs, 'VariableNames', {'metric','scope','session','n_s1','n_sk','mean_s1','mean_sk', ...
    'diff_vs_s1','p_ranksum','p_ranksum_holm','p_ttest','p_ttest_holm'}), 'hw2_stats_vs_session1.csv');
fprintf('\nSaved figures in %s and CSVs: hw2_session_means, hw2_stats_trend, hw2_stats_vs_session1\n', figDir);
diary off;

%% ---------- local functions ----------
function plotSessionLines(M, metric, ylab, cols, showChance)
    hold on; h = gobjects(1,2);
    for si = 1:2
        d = M(M.subject == si, :); jit = (si - 1.5) * 0.12;
        mu = arrayfun(@(k) mean(d.(metric)(d.session == k), 'omitnan'), 1:6);
        sd = arrayfun(@(k) std(d.(metric)(d.session == k), 'omitnan'), 1:6);
        scatter(d.session + jit, d.(metric), 28, cols(si,:), 'filled', 'MarkerFaceAlpha', 0.35, 'HandleVisibility', 'off');
        h(si) = errorbar((1:6) + jit, mu, sd, '-o', 'Color', cols(si,:), 'LineWidth', 1.8, 'MarkerFaceColor', cols(si,:), 'CapSize', 8);
    end
    if showChance, yline(0.5, ':k', 'two-class chance (0.5)', 'HandleVisibility', 'off'); end
    xlim([0.5 6.5]); ylim([0 1.05]); xticks(1:6); xticklabels(compose('on-%d', 1:6));
    xlabel('Online session'); ylabel(ylab); legend(h, {'Subject 1','Subject 2'}, 'Location', 'best'); box off; grid on;
end

function drawPanel(M, metric, ylab, stars, r, p)
    hold on; off = [-0.2 0.2]; fc = [0.65 0.65 0.65; 0.9 0.45 0.45];
    hb = gobjects(1,2); mu = nan(6,2); sd = nan(6,2);
    for si = 1:2
        for k = 1:6
            v = M.(metric)(M.subject == si & M.session == k);
            mu(k,si) = mean(v, 'omitnan'); sd(k,si) = std(v, 'omitnan');
        end
        hb(si) = bar((1:6) + off(si), mu(:,si), 0.36, 'FaceColor', fc(si,:), 'EdgeColor', 'none');
        errorbar((1:6) + off(si), mu(:,si), sd(:,si), 'k.', 'LineWidth', 1, 'CapSize', 5, 'HandleVisibility', 'off');
        for k = 1:6
            v = M.(metric)(M.subject == si & M.session == k);
            scatter(k + off(si) + 0.06 * (rand(numel(v),1) - 0.5), v, 18, 'k', 'filled', 'MarkerFaceAlpha', 0.6, 'HandleVisibility', 'off');
        end
    end
    ok = ~isnan(M.(metric)); cf = polyfit(M.session(ok), M.(metric)(ok), 1);
    ht = plot(1:6, polyval(cf, 1:6), '-', 'Color', [0.1 0.65 0.2], 'LineWidth', 2.2);
    for k = 2:6
        if ~isempty(stars{k}), text(k, 1.08, stars{k}, 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', 'r', 'FontWeight', 'bold'); end
    end
    text(0.55, 1.17, sprintf('Group trend: R = %.2f, p = %.3g', r, p), 'FontSize', 10);
    ylim([0 1.25]); xlim([0.4 6.6]); xticks(1:6); xticklabels(compose('on-%d', 1:6));
    xlabel('Online session'); ylabel(ylab); box off;
    legend([hb ht], {'Subject 1','Subject 2','Group linear fit'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
end

function plotTrend(M, metric, ylab, cols, trend, mi)
    hold on; h = gobjects(1,3);
    for si = 1:2
        d = M(M.subject == si, :); ok = ~isnan(d.(metric));
        scatter(d.session + 0.08 * (rand(height(d),1) - 0.5), d.(metric), 30, cols(si,:), 'filled', 'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
        cf = polyfit(d.session(ok), d.(metric)(ok), 1);
        h(si) = plot(1:6, polyval(cf, 1:6), '-', 'Color', cols(si,:), 'LineWidth', 2);
    end
    ok = ~isnan(M.(metric)); cf = polyfit(M.session(ok), M.(metric)(ok), 1);
    h(3) = plot(1:6, polyval(cf, 1:6), '--k', 'LineWidth', 2);
    lab = cell(1,3);
    for sc = 1:3
        t = trend(mi,sc).t; nm = {'Subject 1','Subject 2','Group'};
        lab{sc} = sprintf('%s (r = %.2f, p = %.3g, n = %d)', nm{sc}, t.r, t.p, t.n);
    end
    legend(h, lab, 'Location', 'southoutside'); xlim([0.5 6.5]); ylim([0 1.05]);
    xticks(1:6); xticklabels(compose('on-%d', 1:6)); xlabel('Online session'); ylabel(ylab); box off; grid on;
end

function pa = holm(p)
    m = numel(p); [ps, ord] = sort(p); adj = min(cummax((m:-1:1) .* ps), 1);
    pa = nan(size(p)); pa(ord) = adj;
end

function [first, sustained] = firstSustained(p, alpha)
    sig = p < alpha; sig(isnan(p)) = false; first = find(sig, 1, 'first');
    if isempty(first), sustained = false; else, sustained = all(sig(first:end)); end
end

function s = starString(p)
    if isnan(p) || p >= 0.05, s = ''; elseif p < 0.001, s = '***'; elseif p < 0.01, s = '**'; else, s = '*'; end
end

function s = yesno(tf)
    if tf, s = 'YES'; else, s = 'no'; end
end

function saveFig(fig, path)
    try, exportgraphics(fig, path, 'Resolution', 200);
    catch, saveas(fig, path); end
end