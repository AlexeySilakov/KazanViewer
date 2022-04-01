function [ax,y,dsc] = kv_cstdata(filename)
% load homogenious datasets
ax.x = [];
ax.y = [];
y = [];
dsc = [];

delim1 = '----------------------------------------------------------------------';
sd1 = length(delim1);
delim2 = sprintf('\r\n\r\n');
fid = fopen(filename, 'r')
str = fread(fid, 'char');
fclose(fid);
str = char(str.');


idx1 = strfind(str, delim1)+sd1;
idx2 = strfind(str, delim2);
head = [1, idx2(1:(end-1))];

totnum = 0;
ax.y = [];
for ii = 1:length(idx1)
    ff = sscanf(str(idx1(ii):idx2(ii)), '%f');
    ss = length(ff);
    if ii>1 & totnum~=ss
        error('File contains datasets of variable length, we can not deal with that at this point');
    end
    nums = reshape(ff, 2, ss/2).';
    y(:, ii) = nums(:, 2);
    if ii==1,
        ax.x = nums(:, 1);
    end
    ax.y(end+1) = ii;
    totnum = ss;
    headr = strtrim(str([head(ii):(idx1(ii)-sd1-2)]));
    dsc = setfield(dsc, ['c', num2str(ii)], headr);
end

