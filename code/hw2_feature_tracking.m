%% hw2_feature_tracking.m  -  HW2 tasks 2.2.4 and 2.2.5
% 2.2.4: per subject, pick the single channel/band feature with the highest Fisher score in online session 6 (variant A),
%        track it over online sessions 1-6, test the trend per subject, and compare session 1 vs 6 at group level.
%        Selection on session 6 is circular, so a leave-one-run-out check (CV) is reported next to the nominal result.
% 2.2.5: Subject 2 only: top-10 session-6 features, correlated with mean command-delivery accuracy across online sessions.
% Inputs : hw2_fisher_long.csv, hw2_fisher_perrun.csv, hw2_run_metrics.csv
% Outputs: hw2_feature_tracking_results.csv, hw2_feature_corr_2_2_5.csv, hw2_feature_tracking_series.csv,
%          hw2_figs/feat_*.png, hw2_feature_console.txt
clear; clc; close all; rng(1);
if exist('hw2_feature_console.txt','file'), delete('hw2_feature_console.txt'); end
diary('hw2_feature_console.txt'); diary on; cleanupObj = onCleanup(@() diary('off'));
fprintf('hw2_feature_tracking run at %s\n', datestr(now));

ALPHA = 0.05; nFeat = 448;
variants = {'A_win_linear','B_win_log10','C_trial_mean_linear'};
cols = [0.15 0.35 0.75; 0.85 0.15 0.15];
figDir = 'hw2_figs'; if ~exist(figDir,'dir'), mkdir(figDir); end

%% ---------- load ----------
FL = readtable('hw2_fisher_long.csv'); FR = readtable('hw2_fisher_perrun.csv'); RM = readtable('hw2_run_metrics.csv');
FA = fillSess(FL, variants{1});                                    % subject x session(0..6) x feature
FRun = nan(2,7,3,nFeat);                                           % subject x session(0..6) x run x feature (variant A)
FRun(sub2ind(size(FRun), FR.subject, FR.session+1, FR.run, FR.featureIdx)) = FR.F;
AccRun = nan(2,6,3);                                               % subject x online session x run
AccRun(sub2ind(size(AccRun), RM.subject, RM.session, RM.run)) = RM.accuracy;
assert(~any(isnan(FA(:))) && ~any(isnan(FRun(:))) && ~any(isnan(AccRun(:))), 'Missing values after fill');
m = FL(strcmp(string(FL.variant), variants{1}) & FL.subject == 1 & FL.session == 0, :);
[~, o] = sort(m.featureIdx); featLab = string(m.chLabel(o)); featHz = m.freqHz(o);
assert(numel(featLab) == nFeat, 'Feature metadata incomplete');
fprintf('Loaded Fisher (session-level and per-run) and accuracy for 2 subjects.\n');

resRows = {};      % task, scope, version, n, estimate, p, note

%% ---------- 2.2.4 best session-6 feature ----------
fprintf('\n===== 2.2.4 Best session-6 feature per subject (variant A, window-level linear PSD) =====\n');
R = struct();
for si = 1:2
    [bestF, bf] = max(squeeze(FA(si, 7, :)));
    R(si).bf = bf; R(si).label = featLab(bf); R(si).hz = featHz(bf);
    R(si).sessF = reshape(FA(si, 2:7, bf), 1, 6);                  % session-pooled F, on-1..on-6
    R(si).runF  = reshape(FRun(si, 2:7, :, bf), 6, 3);             % per-run F, sessions x runs
    fprintf('Subject %d: feature #%d = %s at %g Hz, session-6 F = %.3f (%s, %s)\n', si, bf, R(si).label, R(si).hz, bestF, ...
        chGroup(R(si).label), bandName(R(si).hz));
    fprintf('  session-pooled F on-1..on-6: %s\n', mat2str(R(si).sessF, 3));
    fprintf('  per-run F (rows = sessions, cols = runs):\n'); disp(R(si).runF);

    % held-out-run check: choose the feature from the other two session-6 runs, score it on the held-out run
    cv6 = nan(1,3); cvLab = strings(1,3);
    for h = 1:3
        others = setdiff(1:3, h);
        mm = reshape(mean(FRun(si, 7, others, :), 3), [], 1);
        [~, fb] = max(mm); cv6(h) = FRun(si, 7, h, fb); cvLab(h) = sprintf('%s-%gHz', featLab(fb), featHz(fb));
    end
    R(si).cv6 = cv6;
    fprintf('  held-out session-6 values: %s | features chosen per fold: %s\n', mat2str(cv6, 3), strjoin(cvLab, ', '));

    x = repmat((1:6)', 1, 3);
    tN = trendStats(x(:), R(si).runF(:));                          % nominal (includes the selection session)
    yCV = R(si).runF; yCV(6,:) = cv6; tC = trendStats(x(:), yCV(:));
    t5 = trendStats(reshape(x(1:5,:),[],1), reshape(R(si).runF(1:5,:),[],1));   % sessions 1-5 only, independent of selection
    tS = trendStats((1:6)', R(si).sessF');                         % session-pooled values, n = 6 (descriptive)
    R(si).tN = tN; R(si).tC = tC;
    fprintf('  trend nominal (n=%d): r=%.3f p=%.4g | Spearman rho=%.3f p=%.4g | slope=%.4f\n', tN.n, tN.r, tN.p, tN.rho, tN.ps, tN.slope);
    fprintf('  trend held-out s6 (n=%d): r=%.3f p=%.4g | Spearman rho=%.3f p=%.4g\n', tC.n, tC.r, tC.p, tC.rho, tC.ps);
    fprintf('  trend sessions 1-5 only (n=%d): r=%.3f p=%.4g\n', t5.n, t5.r, t5.p);
    fprintf('  trend session-pooled F (n=%d): Spearman rho=%.3f p=%.4g\n', tS.n, tS.rho, tS.ps);
    nm = sprintf('Subject %d (%s %gHz)', si, R(si).label, R(si).hz);
    resRows(end+1,:) = {'2.2.4 trend', nm, 'nominal run-level Pearson', tN.n, tN.r, tN.p, sprintf('slope=%.4f; Spearman p=%.4g', tN.slope, tN.ps)}; %#ok<SAGROW>
    resRows(end+1,:) = {'2.2.4 trend', nm, 'held-out session 6 Pearson', tC.n, tC.r, tC.p, sprintf('Spearman p=%.4g', tC.ps)}; %#ok<SAGROW>
    resRows(end+1,:) = {'2.2.4 trend', nm, 'sessions 1-5 only Pearson', t5.n, t5.r, t5.p, 'unaffected by selection session'}; %#ok<SAGROW>
    resRows(end+1,:) = {'2.2.4 trend', nm, 'session-pooled Spearman', tS.n, tS.rho, tS.ps, 'descriptive, n=6'}; %#ok<SAGROW>
end

% group level: session 1 vs session 6 (runs as data points: 3 per subject per session)
a1  = [R(1).runF(1,:)'; R(2).runF(1,:)'];
a6  = [R(1).runF(6,:)'; R(2).runF(6,:)'];
a6c = [R(1).cv6(:);     R(2).cv6(:)];
gN = grpTests(a1, a6); gC = grpTests(a1, a6c);
fprintf('\nGroup on-1 vs on-6 (6 runs vs 6 runs; each subject''s own feature):\n');
fprintf('  nominal : mean diff = %.3f | rank-sum p = %.4g | Welch p = %.4g | subject-adjusted p = %.4g\n', gN.diff, gN.pRank, gN.pWelch, gN.pAdj);
fprintf('  held-out: mean diff = %.3f | rank-sum p = %.4g | Welch p = %.4g | subject-adjusted p = %.4g\n', gC.diff, gC.pRank, gC.pWelch, gC.pAdj);
for k = 1:2
    g = {gN, gC}; vn = {'nominal', 'held-out s6'};
    resRows(end+1,:) = {'2.2.4 group on-1 vs on-6', 'Group', [vn{k} ' rank-sum'],   12, g{k}.diff, g{k}.pRank,  'n=6 vs 6'}; %#ok<SAGROW>
    resRows(end+1,:) = {'2.2.4 group on-1 vs on-6', 'Group', [vn{k} ' Welch t'],    12, g{k}.diff, g{k}.pWelch, 'n=6 vs 6'}; %#ok<SAGROW>
    resRows(end+1,:) = {'2.2.4 group on-1 vs on-6', 'Group', [vn{k} ' subject-adj'], 12, g{k}.diff, g{k}.pAdj,   'fitlm y ~ post + subject'}; %#ok<SAGROW>
end

% robustness: same selection rule under each Fisher variant (descriptive series, offline included)
serRows = {};
fprintf('\nRobustness: best session-6 feature per Fisher variant\n');
for v = 1:3
    FV = fillSess(FL, variants{v});
    for si = 1:2
        [~, fb] = max(squeeze(FV(si, 7, :)));
        fprintf('  %-22s Subject %d: %s %g Hz | F off..on-6 = %s\n', variants{v}, si, featLab(fb), featHz(fb), mat2str(squeeze(FV(si,:,fb))', 3));
        for s = 0:6
            serRows(end+1,:) = {variants{v}, si, char(featLab(fb)), featHz(fb), s, FV(si, s+1, fb)}; %#ok<SAGROW>
        end
    end
end
writetable(cell2table(serRows, 'VariableNames', {'variant','subject','channel','freqHz','session','F'}), 'hw2_feature_tracking_series.csv');

%% ---------- 2.2.5 Subject 2: top-10 session-6 features vs accuracy ----------
fprintf('\n===== 2.2.5 Subject 2: top-10 session-6 features vs command-delivery accuracy =====\n');
[~, ord] = sort(squeeze(FA(2, 7, :)), 'descend'); top10 = ord(1:10);
accSess = reshape(mean(AccRun(2, :, :), 3), 6, 1);                 % mean of 3 runs per online session
accRun  = reshape(AccRun(2, :, :), 6, 3);
sessIdx = repmat((1:6)', 1, 3);
fprintf('Subject 2 session accuracy (on-1..on-6): %s\n', mat2str(accSess', 3));
cRows = nan(10, 8); cLab = strings(10,1);
for j = 1:10
    f = top10(j); cLab(j) = sprintf('%s %gHz', featLab(f), featHz(f));
    Fs = reshape(FA(2, 2:7, f), 6, 1);                             % session-pooled F
    Fr = reshape(FRun(2, 2:7, :, f), 6, 3);                        % per-run F
    [r6, p6]   = corr(Fs, accSess);
    [rho6, ps6] = corr(Fs, accSess, 'Type', 'Spearman');
    [r18, p18] = corr(Fr(:), accRun(:));
    [rp, pp]   = partialcorr(Fr(:), accRun(:), sessIdx(:));        % controls for the shared session trend
    cRows(j,:) = [r6 p6 rho6 ps6 r18 p18 rp pp];
end
pHolm6  = holm(cRows(:,2)'); pHolm18 = holm(cRows(:,6)'); pHolmP = holm(cRows(:,8)');
C = table((1:10)', top10, featLab(top10), featHz(top10), arrayfun(@(k) chGroup(featLab(k)), top10), arrayfun(@(k) bandName(featHz(k)), top10), ...
    cRows(:,1), cRows(:,2), pHolm6', cRows(:,3), cRows(:,4), cRows(:,5), cRows(:,6), pHolm18', cRows(:,7), cRows(:,8), pHolmP', ...
    'VariableNames', {'rank','featureIdx','channel','freqHz','region','band','r_session','p_session','pHolm_session', ...
    'rho_session','pRho_session','r_run','p_run','pHolm_run','r_partial_run','p_partial','pHolm_partial'});
disp(C); writetable(C, 'hw2_feature_corr_2_2_5.csv');
fprintf('Significant after Holm (alpha = %.2f): session-level n=6: %d | run-level n=18: %d | partial (session removed): %d\n', ...
    ALPHA, sum(pHolm6 < ALPHA), sum(pHolm18 < ALPHA), sum(pHolmP < ALPHA));
fprintf('Uncorrected p < %.2f: session-level: %d | run-level: %d | partial: %d\n', ALPHA, sum(cRows(:,2) < ALPHA), sum(cRows(:,6) < ALPHA), sum(cRows(:,8) < ALPHA));
for j = 1:10
    resRows(end+1,:) = {'2.2.5 correlation', cLab(j), 'run-level Pearson (n=18)', 18, cRows(j,5), cRows(j,6), sprintf('Holm p=%.4g', pHolm18(j))}; %#ok<SAGROW>
    resRows(end+1,:) = {'2.2.5 correlation', cLab(j), 'session-level Pearson (n=6)', 6, cRows(j,1), cRows(j,2), sprintf('Holm p=%.4g', pHolm6(j))}; %#ok<SAGROW>
    resRows(end+1,:) = {'2.2.5 correlation', cLab(j), 'partial, session removed (n=18)', 18, cRows(j,7), cRows(j,8), sprintf('Holm p=%.4g', pHolmP(j))}; %#ok<SAGROW>
end
writetable(cell2table(resRows, 'VariableNames', {'task','scope','test','n','estimate','p','note'}), 'hw2_feature_tracking_results.csv');

%% ---------- figures ----------
f1 = figure('Position', [80 80 1250 480]);
for si = 1:2
    subplot(1, 2, si); hold on;
    for k = 1:3, scatter(1:6, R(si).runF(:,k), 30, cols(si,:), 'filled', 'MarkerFaceAlpha', 0.45, 'HandleVisibility', 'off'); end
    h1 = plot(1:6, mean(R(si).runF, 2), '-o', 'Color', cols(si,:), 'LineWidth', 1.8, 'MarkerFaceColor', cols(si,:));
    h2 = plot(1:6, R(si).sessF, '--s', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.3);
    h3 = scatter(6.15 * ones(1,3), R(si).cv6, 60, 'Marker', 'd', 'MarkerEdgeColor', 'k', 'LineWidth', 1.4);
    xlim([0.5 6.6]); xticks(1:6); xticklabels(compose('on-%d', 1:6)); xlabel('Online session'); ylabel('Fisher score'); box off; grid on;
    title(sprintf('Subject %d: %s %g Hz (r = %.2f, p = %.3g, n = 18)', si, R(si).label, R(si).hz, R(si).tN.r, R(si).tN.p));
    legend([h1 h2 h3], {'Mean of run F', 'Session-pooled F', 'Held-out run F at on-6'}, 'Location', 'northwest');
end
sgtitle('2.2.4  Fisher score of the feature selected at on-6 (selection session is biased upward)');
saveFig(f1, fullfile(figDir, 'feat_fig1_tracking.png'));

f2 = figure('Position', [80 80 1000 460]);
gs = {gN, gC}; b6 = {a6, a6c}; tt = {'Nominal on-6 (selected on same data)', 'Held-out on-6 (leave-one-run-out selection)'};
subj = [1 1 1 2 2 2]';
for k = 1:2
    subplot(1, 2, k); hold on;
    for si = 1:2
        scatter(1 + 0.05*randn(3,1), a1(subj==si), 55, cols(si,:), 'filled');
        scatter(2 + 0.05*randn(3,1), b6{k}(subj==si), 55, cols(si,:), 'filled');
    end
    plot([1 2], [mean(a1) mean(b6{k})], 'k-', 'LineWidth', 2);
    xlim([0.5 2.5]); xticks([1 2]); xticklabels({'on-1','on-6'}); ylabel('Fisher score of selected feature'); box off; grid on;
    title(sprintf('%s\nrank-sum p = %.3g | Welch p = %.3g | subject-adj p = %.3g', tt{k}, gs{k}.pRank, gs{k}.pWelch, gs{k}.pAdj), 'FontSize', 9);
end
saveFig(f2, fullfile(figDir, 'feat_fig2_group_on1_vs_on6.png'));

f3 = figure('Position', [80 80 1100 480]); hold on;
bar(1:10, [C.r_session C.r_run C.r_partial_run], 'grouped');
for j = 1:10
    if pHolm18(j) < ALPHA, text(j, max(C.r_run(j), 0) + 0.05, '*', 'FontSize', 18, 'HorizontalAlignment', 'center'); end
end
xticks(1:10); xticklabels(cLab); xtickangle(35); ylabel('Correlation with accuracy'); ylim([-0.6 1]); box off; grid on;
legend({'Session-level (n=6)', 'Run-level (n=18)', 'Run-level, session removed'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
title('2.2.5  Subject 2: top-10 on-6 features vs accuracy (* = run-level Holm p < 0.05)');
saveFig(f3, fullfile(figDir, 'feat_fig3_corr_2_2_5.png'));

fprintf('\nSaved CSVs and figures in %s\n', figDir);
diary off;

%% ---------- local functions ----------
function F = fillSess(FL, variant)
    r = FL(strcmp(string(FL.variant), variant), :);
    F = nan(2, 7, 448); F(sub2ind(size(F), r.subject, r.session+1, r.featureIdx)) = r.F;
end

function t = trendStats(x, y)
    x = x(:); y = y(:); [r, p] = corr(x, y); [rho, ps] = corr(x, y, 'Type', 'Spearman'); mdl = fitlm(x, y);
    t = struct('n', numel(x), 'r', r, 'p', p, 'rho', rho, 'ps', ps, 'slope', mdl.Coefficients.Estimate(2));
end

function out = grpTests(a, b)       % a, b: 6x1, first 3 = subject 1, last 3 = subject 2
    d = table([a; b], [zeros(6,1); ones(6,1)], categorical([1 1 1 2 2 2 1 1 1 2 2 2]'), 'VariableNames', {'y','post','subj'});
    mdl = fitlm(d, 'y ~ post + subj'); [~, pw] = ttest2(a, b, 'Vartype', 'unequal');
    out = struct('pRank', ranksum(a, b), 'pWelch', pw, 'pAdj', mdl.Coefficients{'post','pValue'}, 'diff', mean(b) - mean(a));
end

function pa = holm(p)
    m = numel(p); [ps, ord] = sort(p); adj = min(cummax((m:-1:1) .* ps), 1); pa = nan(size(p)); pa(ord) = adj;
end

function s = bandName(hz)
    if hz < 8, s = "theta (4-6)"; elseif hz <= 12, s = "mu/alpha (8-12)"; elseif hz <= 20, s = "low beta (14-20)"; else, s = "high beta (22-30)"; end
end

function s = chGroup(lbl)
    lbl = char(lbl);
    if ismember(lbl, {'C3','Cz','C4','FC1','FC2','CP1','CP2'}), s = "central motor";
    elseif ismember(lbl, {'FC5','FC6','CP5','CP6'}), s = "lateral FC/CP";
    elseif ismember(lbl, {'Fp1','Fpz','Fp2','F7','F8','T7','T8','M1','M2'}), s = "frontal/temporal/mastoid";
    elseif ismember(lbl, {'O1','Oz','O2','POz'}), s = "occipital";
    else, s = "other frontal/parietal"; end
end

function saveFig(fig, path)
    axs = findall(fig, 'Type', 'axes');                 % hide axes toolbars so exported PNGs have no toolbar overlay
    for k = 1:numel(axs)
        try, axs(k).Toolbar.Visible = 'off'; catch, end
    end
    drawnow;
    try, exportgraphics(fig, path, 'Resolution', 200);
    catch, saveas(fig, path); end
end