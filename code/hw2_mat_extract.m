%% hw2_extract.m  -  BCI HW2: extract everything needed from subj1.mat, subj2.mat, ch32Locations.mat
% Outputs: hw2_console.txt (diary), hw2_extracted.mat (local only, large), hw2_trials.csv, hw2_run_metrics.csv
% Place next to subj1.mat, subj2.mat, ch32Locations.mat. Requires MATLAB R2017b+ (mink, local functions in scripts).
clear; clc;
if exist('hw2_console.txt','file'), delete('hw2_console.txt'); end
diary('hw2_console.txt'); diary on;
cleanupObj = onCleanup(@() diary('off'));

%% ---------------- CONFIG ----------------
subjFiles   = {'subj1.mat','subj2.mat'};
locFile     = 'ch32Locations.mat';
TIMEOUT_S   = 7;      % online timeout (HW Table 1)
TIMEOUT_TOL = 0.25;   % s; 7692/7702 with duration >= TIMEOUT_S - TOL counts as timeout. CHECK histogram below.
OFFLINE_S   = 5;      % offline fixed task length
CODE_TRIAL=1000; CODE_FIX=768; CODE_CUE=[769 770]; CODE_START=[7691 7701];
CODE_END_MISS=[7692 7702]; CODE_END_HIT=[7693 7703]; CODE_RUN=32766;
KNOWN = [CODE_TRIAL CODE_FIX CODE_CUE CODE_START CODE_END_MISS CODE_END_HIT CODE_RUN];
fprintf('hw2_extract run at %s\n', datestr(now));
fprintf('Config: TIMEOUT_S=%g TIMEOUT_TOL=%g OFFLINE_S=%g\n', TIMEOUT_S, TIMEOUT_TOL, OFFLINE_S);

%% ---------------- ELECTRODE LOCATIONS ----------------
L = load(locFile);
fn = fieldnames(L);
locVar = '';
for i = 1:numel(fn)
    v = L.(fn{i});
    if isstruct(v) && isfield(v,'labels'), locVar = fn{i}; break; end
end
assert(~isempty(locVar), 'No EEGLAB-style loc struct (with .labels) in %s. Vars: %s', locFile, strjoin(fn',', '));
chanlocs   = L.(locVar);
locLabels  = {chanlocs.labels}';
fprintf('ch32Locations: variable "%s", %d channels\n', locVar, numel(chanlocs));
fprintf('Loc labels: %s\n', strjoin(locLabels', ','));
fprintf('Loc fields: %s\n', strjoin(fieldnames(chanlocs)', ','));
whos('-file', locFile)

%% ---------------- LOAD SUBJECTS ----------------
D = struct();
D.chanlocs = chanlocs; D.locLabels = locLabels;
D.params = struct('TIMEOUT_S',TIMEOUT_S,'TIMEOUT_TOL',TIMEOUT_TOL,'OFFLINE_S',OFFLINE_S);
trialRows = {};   % per-trial rows
qc = {};          % QC messages
fsAll = [];       % every run's fs

for si = 1:numel(subjFiles)
    fprintf('\n===== %s =====\n', subjFiles{si});
    whos('-file', subjFiles{si})
    S = load(subjFiles{si});
    sfn = fieldnames(S);
    subj = [];
    for i = 1:numel(sfn)
        v = S.(sfn{i});
        if isstruct(v) && isfield(v,'offline') && isfield(v,'online')
            subj = v; fprintf('Top-level var: %s\n', sfn{i}); break;
        end
    end
    assert(~isempty(subj), 'No variable with .offline and .online in %s', subjFiles{si});

    nOn = numel(subj.online);
    fprintf('Online sessions found: %d\n', nOn);
    firstRun = subj.offline.run(1);
    fprintf('Offline run(1) fields: %s\n', strjoin(fieldnames(firstRun)', ','));
    if isfield(firstRun,'header') && isstruct(firstRun.header)
        fprintf('Offline run(1).header fields: %s\n', strjoin(fieldnames(firstRun.header)', ','));
    end

    sessList = [0, 1:nOn];                    % 0 = offline
    for sIdx = 1:numel(sessList)
        s = sessList(sIdx);
        if s == 0, sessStruct = subj.offline; else, sessStruct = subj.online(s); end
        nRuns = numel(sessStruct.run);
        for r = 1:nRuns
            R = sessStruct.run(r);
            eeg = double(R.eeg);
            fs      = pick(R, {'fs','header.fs','header.SampleRate'});
            chLab   = pick(R, {'chLabels','chLabel','header.chLabels','header.chLabel'});
            typ = pick(R, {'header.triggers.TYP','header.EVENT.TYP'}); typ = double(typ(:));
            pos = pick(R, {'header.triggers.POS','header.EVENT.POS'}); pos = double(pos(:));
            if size(eeg,2) ~= 32 && size(eeg,1) == 32
                eeg = eeg'; qc{end+1} = sprintf('S%d s%d r%d: EEG transposed', si,s,r); %#ok<SAGROW>
            end
            if iscell(chLab), chLab = chLab(:); else, chLab = cellstr(chLab); end
            fsAll(end+1) = fs; %#ok<SAGROW>

            D.subj(si).sess(sIdx).session = s;
            D.subj(si).sess(sIdx).run(r).eeg = eeg;
            D.subj(si).sess(sIdx).run(r).fs = fs;
            D.subj(si).sess(sIdx).run(r).chLabels = chLab;
            D.subj(si).sess(sIdx).run(r).typ = typ;
            D.subj(si).sess(sIdx).run(r).pos = pos;

            if size(eeg,2) ~= 32, qc{end+1} = sprintf('S%d s%d r%d: %d channels (expected 32)', si,s,r,size(eeg,2)); end %#ok<SAGROW>
            if ~isequal(chLab(:), locLabels(:)), qc{end+1} = sprintf('S%d s%d r%d: chLabels differ from ch32Locations order', si,s,r); end %#ok<SAGROW>
            unk = setdiff(unique(typ), KNOWN);
            if ~isempty(unk), qc{end+1} = sprintf('S%d s%d r%d: unknown codes %s', si,s,r,mat2str(unk')); end %#ok<SAGROW>
            if any(pos<1 | pos>size(eeg,1)), qc{end+1} = sprintf('S%d s%d r%d: POS out of EEG range', si,s,r); end %#ok<SAGROW>
            if numel(typ) ~= numel(pos), qc{end+1} = sprintf('S%d s%d r%d: TYP/POS length mismatch', si,s,r); end %#ok<SAGROW>

            iStart = find(ismember(typ, CODE_START));
            nBad = 0;
            for k = 1:numel(iStart)
                i0 = iStart(k);
                nxt = find(ismember(typ(i0+1:end), [CODE_END_MISS CODE_END_HIT]), 1, 'first');
                if isempty(nxt), nBad = nBad+1; continue; end
                i1 = i0 + nxt;
                if any(ismember(typ(i0+1:i1-1), [CODE_START CODE_TRIAL])), nBad = nBad+1; continue; end
                iCue = find(ismember(typ(1:i0), CODE_CUE), 1, 'last');
                iTr  = find(typ(1:i0) == CODE_TRIAL, 1, 'last');

                cls = 'LH'; if typ(i0) == 7701, cls = 'RH'; end
                endCode = typ(i1);
                p0 = pos(i0); p1 = pos(i1);
                dur = (p1 - p0) / fs;
                if s == 0
                    outcome = 'offline';
                elseif ismember(endCode, CODE_END_HIT)
                    outcome = 'hit';
                elseif dur >= TIMEOUT_S - TIMEOUT_TOL
                    outcome = 'timeout';
                else
                    outcome = 'miss';
                end
                cueP = NaN; if ~isempty(iCue), cueP = pos(iCue); end
                trP  = NaN; if ~isempty(iTr),  trP  = pos(iTr);  end
                trialRows(end+1,:) = {si, s, r, k, cls, double(typ(i0)), double(endCode), outcome, ...
                                      p0, p1, dur, cueP, trP}; %#ok<SAGROW>
            end
            if nBad>0, qc{end+1} = sprintf('S%d s%d r%d: %d incomplete trials skipped', si,s,r,nBad); end %#ok<SAGROW>
            fprintf('S%d sess %d run %d: %d x %d EEG, fs=%g, %d events, %d task starts\n', ...
                    si, s, r, size(eeg,1), size(eeg,2), fs, numel(typ), numel(iStart));
        end
    end
end

%% ---------------- TRIAL TABLE ----------------
T = cell2table(trialRows, 'VariableNames', {'subject','session','run','trialIdx','class', ...
    'startCode','endCode','outcome','taskStartSample','taskEndSample','durationS','cueSample','trialStartSample'});
D.trials = T;
writetable(T, 'hw2_trials.csv');

%% ---------------- QC: TRIAL COUNTS ----------------
fprintf('\nTrials per subject/session/run/class:\n');
disp(groupcounts(T, {'subject','session','run','class'}));
runCnt = groupcounts(T, {'subject','session','run'});
bad = runCnt(runCnt.GroupCount ~= 20,:);
if ~isempty(bad)
    qc{end+1} = sprintf('%d runs do not have exactly 20 trials', height(bad));
    disp(bad);
end
clsCnt = groupcounts(T, {'subject','session','run','class'});
if any(clsCnt.GroupCount ~= 10), qc{end+1} = 'Some run/class groups do not have exactly 10 trials'; end

%% ---------------- ONLINE METRICS (2.1.1) ----------------
On = T(T.session > 0,:);
[g, sj, ss, rr] = findgroups(On.subject, On.session, On.run);   % output sorted by subject, session, run
nH  = splitapply(@(x) sum(strcmp(x,'hit')),     On.outcome, g);
nM  = splitapply(@(x) sum(strcmp(x,'miss')),    On.outcome, g);
nTO = splitapply(@(x) sum(strcmp(x,'timeout')), On.outcome, g);
nTr = splitapply(@numel, On.outcome, g);
acc = nH ./ (nH + nM);                  % NaN if no hits and no misses
pTO = nTO ./ nTr;
M = table(sj, ss, rr, nH, nM, nTO, nTr, acc, pTO, 'VariableNames', ...
          {'subject','session','run','hits','misses','timeouts','trials','accuracy','pctTimeout'});
D.runMetrics = M;
writetable(M, 'hw2_run_metrics.csv');

for si = 1:numel(subjFiles)
    m = M(M.subject==si,:);
    if ~issorted([m.session m.run], 'rows')
        qc{end+1} = sprintf('Subject %d: run metrics not sorted, accMat unreliable', si); %#ok<SAGROW>
    end
    D.accMat{si} = reshape(m.accuracy,   [], max(m.session))';   % rows = sessions, cols = runs
    D.toMat{si}  = reshape(m.pctTimeout, [], max(m.session))';
    fprintf('\nSubject %d accuracy (rows=sessions, cols=runs):\n', si); disp(D.accMat{si});
    fprintf('Subject %d fraction timeout (rows=sessions, cols=runs):\n', si); disp(D.toMat{si});
end

%% ---------------- TASK-PERIOD SEGMENTS (for 2.2 PSD) ----------------
nT = height(T);
seg = repmat(struct('subject',[],'session',[],'run',[],'class',[],'outcome',[],'x',[]), nT, 1);
for k = 1:nT
    si = T.subject(k); s = T.session(k); r = T.run(k);
    sIdx = find([D.subj(si).sess.session] == s);
    X = D.subj(si).sess(sIdx).run(r).eeg;
    i0 = T.taskStartSample(k); i1 = min(T.taskEndSample(k), size(X,1));
    seg(k).subject = si; seg(k).session = s; seg(k).run = r;
    seg(k).class = T.class{k}; seg(k).outcome = T.outcome{k};
    seg(k).x = X(i0:i1, :);
end
D.seg = seg;

%% ---------------- QC: CHANNELS, DURATIONS, FS ----------------
for si = 1:numel(subjFiles)
    allVar = zeros(1,32); n = 0;
    for sIdx = 1:numel(D.subj(si).sess)
        for r = 1:numel(D.subj(si).sess(sIdx).run)
            allVar = allVar + var(D.subj(si).sess(sIdx).run(r).eeg, 0, 1); n = n+1;
        end
    end
    allVar = allVar / n;
    [~, flat] = mink(allVar, 3);
    fprintf('Subject %d lowest-variance channels: %s (mean var %s)\n', si, strjoin(locLabels(flat)',','), mat2str(allVar(flat),3));
    D.chVar{si} = allVar;
end
fprintf('\nUnique fs across all runs: %s\n', mat2str(unique(fsAll)));
if numel(unique(fsAll)) > 1, qc{end+1} = 'More than one sampling rate across runs'; end

fprintf('\nOnline 7692/7702 task durations (s), to choose timeout cutoff:\n');
d = On.durationS(ismember(On.endCode, CODE_END_MISS));
fprintf('  n=%d, min=%.2f, median=%.2f, max=%.2f\n', numel(d), min(d), median(d), max(d));
fprintf('  histcounts over bins [0,1,...,8] s: %s\n', mat2str(histcounts(d, 0:1:8)));
dh = On.durationS(ismember(On.endCode, CODE_END_HIT));
if ~isempty(dh), fprintf('Online hit durations: min=%.2f, median=%.2f, max=%.2f\n', min(dh), median(dh), max(dh)); end
fprintf('Offline durations: median %.2f s (expect %g)\n', median(T.durationS(T.session==0)), OFFLINE_S);

fprintf('\nOutcome totals (online):\n'); disp(groupcounts(On, {'subject','session','outcome'}));

%% ---------------- REPORT & SAVE ----------------
fprintf('\n===== QC FLAGS =====\n');
if isempty(qc), disp('none'); else, cellfun(@(q) fprintf('%s\n',q), qc); end
D.qc = qc;
save('hw2_extracted.mat', 'D', '-v7.3');
fprintf('\nSaved hw2_extracted.mat (local only), hw2_trials.csv, hw2_run_metrics.csv, hw2_console.txt\n');
diary off;

%% ---------------- LOCAL FUNCTION ----------------
function val = pick(s, paths)
    for p = 1:numel(paths)
        parts = strsplit(paths{p}, '.'); v = s; ok = true;
        for j = 1:numel(parts)
            if isstruct(v) && isfield(v, parts{j}), v = v.(parts{j}); else, ok = false; break; end
        end
        if ok, val = v; return; end
    end
    error('None of these fields found: %s', strjoin(paths, ', '));
end