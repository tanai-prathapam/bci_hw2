%% hw2_psd_fisher.m  -  HW2 2.2.1 (PSD features) and 2.2.2 (Fisher scores)
% Features: pwelch, 1 s Hamming-windowed segments shifted 62.5 ms, 0.5 s sub-windows (50% overlap),
% 4:2:30 Hz x 32 channels = 448 features. Feature index = (channel-1)*14 + freqIdx.
% Fisher (HW Eq. 1): F = (mu1-mu2)^2 / (var1+var2), classes LH vs RH.
clear; clc;
if exist('hw2_psd_console.txt','file'), delete('hw2_psd_console.txt'); end
diary('hw2_psd_console.txt'); diary on; cleanupObj = onCleanup(@() diary('off'));

subjFiles = {'subj1.mat','subj2.mat'}; subjVars = {'subj1','subj2'};
P.WIN_S = 1; P.SHIFT_S = 0.0625; P.SUB_S = 0.5; P.SUB_OVL = 0.5; P.FREQS = 4:2:30;
fprintf('hw2_psd_fisher run at %s\n', datestr(now)); disp(P);

T = readtable('hw2_trials.csv');
Lc = load('ch32Locations.mat'); locLabels = {Lc.ch32Locations.labels}';
nCh = 32; freqs = P.FREQS; nF = numel(freqs); nFeat = nCh*nF;
assert(nFeat == 448, 'Expected 448 features');
featCh = repelem((1:nCh)', nF); featF = repmat((1:nF)', nCh, 1); featHz = freqs(featF)';

%% ---------- 2.2.1 PSD features per window ----------
allF = {}; allMeta = {};
for si = 1:2
    A = load(subjFiles{si}); S = A.(subjVars{si}); clear A
    for s = 0:6
        if s == 0, sess = S.offline; else, sess = S.online(s); end
        for r = 1:numel(sess.run)
            R = sess.run(r); X = double(R.eeg); fs = R.header.fs;
            assert(fs == 512 && size(X,2) == 32, 'Unexpected fs or channel count');
            winN = round(P.WIN_S*fs); shiftN = round(P.SHIFT_S*fs);
            subN = round(P.SUB_S*fs);  ovlN = round(P.SUB_OVL*subN);
            for k = find(T.subject == si & T.session == s & T.run == r)'
                p0 = T.taskStartSample(k); p1 = min(T.taskEndSample(k), size(X,1));
                starts = p0:shiftN:(p1 - winN + 1);
                if isempty(starts), continue; end
                F = zeros(numel(starts), nFeat, 'single');
                for w = 1:numel(starts)
                    seg = X(starts(w):starts(w)+winN-1, :);
                    pxx = pwelch(seg, hamming(subN), ovlN, freqs, fs);   % nF x nCh
                    F(w,:) = single(pxx(:))';
                end
                isLH = strcmp(string(T.class(k)), "LH");
                allF{end+1} = F; %#ok<SAGROW>
                allMeta{end+1} = repmat([si s r T.trialIdx(k) isLH], numel(starts), 1); %#ok<SAGROW>
            end
            fprintf('S%d sess %d run %d done\n', si, s, r);
        end
    end
    clear S
end
FEAT = vertcat(allF{:}); META = vertcat(allMeta{:}); clear allF allMeta
N = size(FEAT,1);
fprintf('\nFEAT: %d windows x %d features | any NaN: %d | any Inf: %d | any <=0: %d\n', ...
    N, size(FEAT,2), any(isnan(FEAT(:))), any(isinf(FEAT(:))), any(FEAT(:) <= 0));
fprintf('PSD min/median/max: %.3g / %.3g / %.3g\n', min(FEAT(:)), median(FEAT(:)), max(FEAT(:)));

% trial-mean features (variant C)
[g, ~] = findgroups(META(:,1), META(:,2), META(:,3), META(:,4));
nTr = max(g); Msp = sparse(g, (1:N)', 1, nTr, N); cnt = full(sum(Msp,2));
TM = (Msp * double(FEAT)) ./ cnt;
TMeta = zeros(nTr,5); for j = 1:5, TMeta(:,j) = splitapply(@(x) x(1), META(:,j), g); end

%% ---------- counts ----------
cRows = [];
for si = 1:2, for s = 0:6, for c = [1 0]
    nTrialsAll = sum(T.subject==si & T.session==s & strcmp(string(T.class), ternary(c)));
    m = META(:,1)==si & META(:,2)==s & META(:,5)==c;
    nTrialsWin = sum(TMeta(:,1)==si & TMeta(:,2)==s & TMeta(:,5)==c);
    cRows(end+1,:) = [si s c nTrialsAll nTrialsWin sum(m)]; %#ok<SAGROW>
end, end, end
Ct = array2table(cRows, 'VariableNames', {'subject','session','isLH','trialsTotal','trialsWithWindows','windows'});
writetable(Ct, 'hw2_psd_counts.csv'); disp(Ct);
off = Ct(Ct.session==0,:); fprintf('Offline windows per trial (expect 65): %s\n', mat2str(off.windows ./ off.trialsWithWindows, 4)');

%% ---------- 2.2.2 Fisher scores ----------
vars = {'A_win_linear','B_win_log10','C_trial_mean_linear'};
longT = {}; topT = {}; chT = {};
for v = 1:3
    for si = 1:2
        for s = 0:6
            if v < 3
                rows = META(:,1)==si & META(:,2)==s; X = double(FEAT(rows,:)); lh = META(rows,5)==1;
                if v == 2, X = log10(X); end
            else
                rows = TMeta(:,1)==si & TMeta(:,2)==s; X = TM(rows,:); lh = TMeta(rows,5)==1;
            end
            Fv = fisher(X, lh);
            nBad = sum(~isfinite(Fv)); if nBad > 0, fprintf('WARN %s S%d s%d: %d non-finite F\n', vars{v}, si, s, nBad); end
            longT{end+1} = table(repmat(string(vars{v}),nFeat,1), repmat(si,nFeat,1), repmat(s,nFeat,1), ...
                (1:nFeat)', featCh, string(locLabels(featCh)), featHz, Fv(:), 'VariableNames', ...
                {'variant','subject','session','featureIdx','channel','chLabel','freqHz','F'}); %#ok<SAGROW>
            Fs = Fv; Fs(~isfinite(Fs)) = -Inf; [~, ord] = sort(Fs, 'descend'); ord = ord(1:10);
            topT{end+1} = table(repmat(string(vars{v}),10,1), repmat(si,10,1), repmat(s,10,1), (1:10)', ...
                ord(:), string(locLabels(featCh(ord))), featHz(ord), Fv(ord)', 'VariableNames', ...
                {'variant','subject','session','rank','featureIdx','chLabel','freqHz','F'}); %#ok<SAGROW>
            D = sum(reshape(Fv, nF, nCh), 1)';
            chT{end+1} = table(repmat(string(vars{v}),nCh,1), repmat(si,nCh,1), repmat(s,nCh,1), ...
                (1:nCh)', string(locLabels), D, 'VariableNames', ...
                {'variant','subject','session','channel','chLabel','sumF'}); %#ok<SAGROW>
        end
    end
end
writetable(vertcat(longT{:}), 'hw2_fisher_long.csv');
writetable(vertcat(topT{:}),  'hw2_top10.csv');
writetable(vertcat(chT{:}),   'hw2_channel_sumF.csv');

% per-run Fisher, variant A (replicates for 2.2.4 and 2.2.5)
prT = {};
for si = 1:2, for s = 0:6, for r = 1:3
    rows = META(:,1)==si & META(:,2)==s & META(:,3)==r;
    Fv = fisher(double(FEAT(rows,:)), META(rows,5)==1);
    prT{end+1} = table(repmat(si,nFeat,1), repmat(s,nFeat,1), repmat(r,nFeat,1), (1:nFeat)', Fv(:), ...
        'VariableNames', {'subject','session','run','featureIdx','F'}); %#ok<SAGROW>
end, end, end
writetable(vertcat(prT{:}), 'hw2_fisher_perrun.csv');

fprintf('\nSubject 2 top-5, variant A, per session:\n');
tt = vertcat(topT{:});
disp(tt(tt.variant=="A_win_linear" & tt.subject==2 & tt.rank<=5, :));

save('hw2_psd_features.mat', 'FEAT', 'META', 'TM', 'TMeta', 'P', 'featCh', 'featHz', 'locLabels', '-v7.3');
fprintf('\nSaved CSVs: hw2_psd_counts, hw2_fisher_long, hw2_top10, hw2_channel_sumF, hw2_fisher_perrun\n');
diary off;

%% ---------- local functions ----------
function f = fisher(X, isLH)
    a = X(isLH,:); b = X(~isLH,:);
    f = (mean(a,1) - mean(b,1)).^2 ./ (var(a,0,1) + var(b,0,1));
end
function s = ternary(c)
    if c == 1, s = "LH"; else, s = "RH"; end
end