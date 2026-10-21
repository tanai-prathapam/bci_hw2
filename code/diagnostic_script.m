A = load('ch32Locations.mat'); locLab = {A.ch32Locations.labels}';
B1 = load('subj1.mat'); B2 = load('subj2.mat'); subjs = {B1.subj1, B2.subj2};
fid = fopen('hw2_channel_audit.txt','w');
for si = 1:2
    S = subjs{si}; ref = S.offline.run(1).header.chLabels;
    lab = ref; if ~iscell(lab), lab = cellstr(lab); end; lab = strtrim(lab(:));
    fprintf(fid,'S%d chLabels class=%s size=%s\n', si, class(ref), mat2str(size(ref)));
    fprintf(fid,'S%d data labels: %s\n', si, strjoin(lab',','));
    [tf,idx] = ismember(upper(lab), upper(locLab));
    fprintf(fid,'S%d all found in loc: %d | identical order: %d\n', si, all(tf), isequal(idx,(1:32)'));
    fprintf(fid,'S%d data->loc index: %s\n', si, mat2str(idx'));
    same = true;
    for r = 1:3, same = same && isequal(S.offline.run(r).header.chLabels, ref); end
    for s = 1:6, for r = 1:3, same = same && isequal(S.online(s).run(r).header.chLabels, ref); end, end
    fprintf(fid,'S%d labels identical across all 21 runs: %d\n', si, same);
end
R = subjs{2}.online(3).run(2); typ = R.header.triggers.TYP(:); pos = R.header.triggers.POS(:);
b = find(pos < 1 | pos > size(R.eeg,1));
fprintf(fid,'S2 s3 r2 EEG length %d | bad idx %s | TYP %s | POS %s\n', size(R.eeg,1), mat2str(b'), mat2str(typ(b)'), mat2str(pos(b)'));
fclose(fid);