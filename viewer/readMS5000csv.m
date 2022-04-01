function varargout = readMS5000csv(filename)

% filename = 'd:\PSU\Data\FeHase\EPR\Patrick\2019_05_10\20190510_181450531_Cb_WT_CO_air_new_40K_10dB_01.csv';

fid = fopen(filename, 'r');
if fid == -1, error('sorry, file not found'); return; end

str = fread(fid, 'char');
fclose(fid);

allascii = char(str)';


fieldnames = {'Recipe', 'Additional', 'Meas'};

linebreak = sprintf('\r\n'); %CR and LF
for ii = 1:numel(fieldnames)
    k(ii) = strfind(allascii, [fieldnames{ii}, linebreak])-1;
end

st = 1;
dsc = struct;

% nlines = strfind(allascii, linebreak)
for ii = 1:(numel(fieldnames))
    lines = allascii(st:k(ii));
    nlines = [0, strfind(lines, linebreak)-1];
    for jj = 2:numel(nlines)
        ss = lines((nlines(jj-1)+2):nlines(jj));
        if isempty(strtrim(ss))
            continue;
        end
        idd = strfind(ss, ';');
        ss(idd) = ' ';
        dsc = setfield(dsc, strtrim(ss(1:idd(1))), ss((idd(1)+1):end));
    end
    st =k(ii)+numel(fieldnames{ii})+2;
end

lines = allascii(st:end);
idx2 = strfind(lines, ';'); lines(idx2) = ' ';
idx = strfind(lines, linebreak); lines(idx) = ' ';lines(idx+1) = ' ';
LL = lines((idx(1)+2):end);
dig = sscanf(LL, '%f');

dig1 = reshape(dig, 2, numel(dig)/2);

y = dig1(2, :).';
ax.x = dig1(1, :).';

ax.type = 'data';
ax.xlabel = 'Magnetic Field, mT';
if isfield(dsc, 'Name')
    ax.title = dsc.Name;
else
    ax.title = '';
end

% assign output depending on number of output arguments
% requested
switch nargout
  case 1
    varargout = {y};
  case 2
    varargout = {ax, y};
  case 3
    varargout = {ax, y, dsc};
end





