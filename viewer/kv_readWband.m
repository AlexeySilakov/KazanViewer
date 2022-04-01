function varargout = kv_readWband(filename)
% Y = kv_readWband(filename);
% [ax, Y] = kv_readWband(filename);
% [ax, Y, dsc] = kv_readWband(filename);
% The script to read W-band files from Berlin
% Designed for Kazan Viewer. 
% Alexey Silakov 12.01.2009
[pathstr_, name_, ext_] = fileparts(filename);

% because the is no structural correlation between different formats, it has to be
% loaded independently


switch ext_
    case '.wtr',
        [ax, Y, dsc] = local_tr(filename);
    case '.w2p',
        [ax, Y, dsc] = local_2p(filename);
    case '.wcw',
        [ax, Y, dsc] = local_cw(filename);
    otherwise
        error('This extention is not supported');
end

%%%% the strings to "description"

fid = fopen(filename, 'r');
if fid == -1, error('sorry, file not found'); return; end
n = 0;
endofcomm =0;
while endofcomm<2
    a = (fgets(fid));
    dsc = setfield(dsc, ['c', num2str(n)], a(1:end-1));
    n = n+1;
    if endofcomm==1; % to get one more string then "% end field"
        endofcomm=2;
    end
    if strcmp(a(1:end-1), '% end field ')
        endofcomm=1;
    end
    if n>300
        disp('N comment strings > 300. No "end field" string within first 300 lines?');
        break;
    end
end
fclose(fid);

% assign output depending on number of output arguments
% requested
switch nargout
  case 1,
    varargout = {Y};
  case 2,
    varargout = {ax, Y};
  case 3,
    varargout = {ax, Y, dsc};
end

function [ax, y, dsc] = local_tr(filename)

% fid = fopen(filename, 'r');
% if fid == -1, error('sorry, file not found'); return; end

% C = textscan(fid, '%f', 'commentStyle', '%');
% Allstuff = C{1};
% fclose(fid);

Allstuff = textread(filename, '%f', 'commentstyle', 'matlab');
%############## time axis ##############
% number of time points 
% ind = find(strcmp(strings, '% number of time points'));
% nTPoints = str2num(strings{ind+1});
nTPoints = Allstuff(1);
% endofdisc = ind-1;

% start time 
% ind = find(strcmp(strings, '% start time'));
% axTS = str2num(strings{ind+1});
axTS = Allstuff(2);
% step time 
% ind = find(strcmp(strings, '% step time'));
% axTStep = str2num(strings{ind+1});
axTStep = Allstuff(3);
% end time 
% ind = find(strcmp(strings, '% end time'));
% axTE = str2num(strings{ind+1});
axTE = Allstuff(4);

axtime = [axTS:axTStep:(axTE-axTStep)].'; % column
%----------------------------------------

%############# field axis ##############
% number of field points 
% ind = find(strcmp(strings, '% number of field points'));
% nFPoints = str2num(strings{ind+1});
nFPoints = Allstuff(5);
% % start field 
% ind = find(strcmp(strings, '% start field'));
% axFS = str2num(strings{ind+1});
axFS = Allstuff(6);
% % step field 
% ind = find(strcmp(strings, '% step field'));
% axFStep = str2num(strings{ind+1});
axFStep = Allstuff(7);
% % end field 
% ind = find(strcmp(strings, '% end field'));
% axFE = str2num(strings{ind+1});
axFE = Allstuff(8);
% dataind = ind+2;

axfield = [axFS:axFStep:axFE].'; % column
%-----------------------------------------

% Data = [];
%############## Data ##############
% for ii = dataind:length(strings)
%     Data(end+1) =str2num(strings{ii});
% end

% y = reshape(Data,nTPoints,nFPoints);

y = reshape(Allstuff(9:end),nTPoints,nFPoints);

ax.type = 'data';
% ax.y = axfield;
ax.y = axfield*108.9+3.4e4;
ax.ylabel = 'Field, G';
ax.x = axtime;
ax.xlabel = 'Time, s';
ax.title = 'Transient';
ax.freq1 = 94.960e+9;
% if fstcol,dsc.FirstColumn='X';else dsc.FirstColumn='Y';end
% dsc.FirstColumn = ['Treated as ', dsc.FirstColumn];
dsc.size_x = num2str(size(y,1));
dsc.size_y = num2str(size(y,2));
%%%% rest of the strings to "description"
% if endofdisc>0
%     for n = 1:(dataind-1)
% %       dsc = setfield(dsc, ['c', num2str(n)], strings{n});
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ax, y, dsc] = local_2p(filename)
% fid = fopen(filename, 'r');
% if fid == -1, error('sorry, file not found'); return; end
% 
% C = textscan(fid, '%f', 'commentStyle', '%');
% Allstuff = C{1};
% 
% fclose(fid);

Allstuff = textread(filename, '%f', 'commentstyle', 'matlab');

%############## time axis ##############

nTPoints = Allstuff(1);

axTS = Allstuff(2);

axTStep = Allstuff(3);

axTE = Allstuff(4);

axtime = [axTS:axTStep:(axTE-axTStep)].'; % column
%----------------------------------------

%############# field axis ##############

nFPoints = Allstuff(5);

axFS = Allstuff(6);

axFStep = Allstuff(7);

axFE = Allstuff(8);

axfield = [axFS:axFStep:axFE].'; % column
%-----------------------------------------

% Data = [];
%############## Data ##############
% for ii = dataind:length(strings)
%     Data(end+1) =str2num(strings{ii});
% end

% y = reshape(Data,nTPoints,nFPoints);
npi = (nTPoints*2+1)*(nFPoints);
ty = reshape(Allstuff(8+[1:npi]),nTPoints*2+1,(nFPoints));
% y = y + i*reshape(Allstuff(8+[1:npi]+npi),nTPoints+1,(nFPoints+1));
y = ty([1:nTPoints]+1, :) + i*ty([1:nTPoints]+nTPoints+1, :);
y = [ty(1, :); y];
ax.type = 'data';
% ax.y = axfield;
ax.y = [axfield*108.9+3.4e4];
ax.ylabel = 'Field, G';
ax.x = [0; axtime];
ax.xlabel = 'Time, s';
ax.title = 'Transient';
ax.freq1 = 94.960e+9;
% if fstcol,dsc.FirstColumn='X';else dsc.FirstColumn='Y';end
% dsc.FirstColumn = ['Treated as ', dsc.FirstColumn];
dsc.size_x = num2str(size(y,1));
dsc.size_y = num2str(size(y,2));
%%%% rest of the strings to "description"
% if endofdisc>0
%     for n = 1:(dataind-1)
% %       dsc = setfield(dsc, ['c', num2str(n)], strings{n});
%     end
% end

function [ax, y, dsc] = local_cw(filename)
% fid = fopen(filename, 'r');
% if fid == -1, error('sorry, file not found'); return; end
% 
% C = textscan(fid, '%f', 'commentStyle', '%');
% Allstuff = C{1};
% fclose(fid);

Allstuff = textread(filename, '%f', 'commentstyle', 'matlab');
%############## field axis ##############

nFPoints = Allstuff(1);

axFS = Allstuff(2);

axFStep = Allstuff(3);

axFE = Allstuff(4);

axfield = [axFS:axFStep:axFE].'; % column
%----------------------------------------

% Data = [];
%############## Data ##############

y = Allstuff(5:end);

ax.type = 'data';
% ax.y = axfield;
ax.x = [axfield*108.9+3.4e4];
ax.xlabel = 'Field, G';
ax.y = [0];
ax.ylabel = 'Time, s';
ax.title = 'Transient';
ax.freq1 = 94.960e+9;
% if fstcol,dsc.FirstColumn='X';else dsc.FirstColumn='Y';end
% dsc.FirstColumn = ['Treated as ', dsc.FirstColumn];
dsc.size_x = num2str(size(y,1));
dsc.size_y = num2str(size(y,2));
%%%% rest of the strings to "description"
% if endofdisc>0
%     for n = 1:(dataind-1)
% %       dsc = setfield(dsc, ['c', num2str(n)], strings{n});
%     end
% end


