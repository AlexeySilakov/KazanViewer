function [ax,y,dsc] = kv_read_RESAktaStart(fname, varargin)


% handles.lscript{1} = ['[ax,y,dsc] = kv_read_RESAktaStart(''', filename, ''',''format'',',str,');'];

fid = fopen(fname, 'r', 'b');
fseek(fid, 0, 'bof');
aa = fread(fid, 'int8');
fseek(fid, 38, 'bof');
Akta = string(fread(fid, 20, 'char').');
RunType = string(fread(fid, 20, 'char').');

fclose(fid);

%num = hex2dec('23521032');%seems like it indicates the beginning/end of a block accept for the end of file
% num = hex2dec('2354'); % big endian "#T"
% num = hex2dec('5423'); % little endian
num1 = hex2dec('23'); % #
num2 = hex2dec('54'); % T
num3 = hex2dec('34'); % 4

cc1 = find(aa==num1);
cc2 = find(aa(cc1+1)==num2);
cc=  cc1(cc2);
cc3 = find(aa(cc+2)==num3);
cc = cc(cc3);
% bl = hex2dec('3030');
bl = hex2dec('30'); %%% somehow the whole deal is upshifted by this number

% figure(53); clf; hold on;

idx = 8;
time = ((aa(cc+idx)-bl)*2^16+ (aa(cc+idx+1)-bl)*2^12 + (aa(cc+idx+2)-bl)*2^6+ aa(cc+idx+3)-bl)/60;%*63.7e-4; % not sure what this is, time?

idx = 14;
ml = ((aa(cc+idx)-bl)*2^16+ (aa(cc+idx+1)-bl)*2^12 + (aa(cc+idx+2)-bl)*2^6+ aa(cc+idx+3)-bl)*1e-4; %ml

idx = 21;
UV = ((aa(cc+idx)-bl)*2^16+ (aa(cc+idx+1)-bl)*2^12 + (aa(cc+idx+2)-bl)*2^6+ aa(cc+idx+3)-bl)*1e-2;

idx = 29;
flow_rate= ((aa(cc+idx)-bl)*2^16+ (aa(cc+idx+1)-bl)*2^12 + (aa(cc+idx+2)-bl)*2^6+ aa(cc+idx+3)-bl)/4e4; % flow rate?

idx = 39;
kPa= ((aa(cc+idx)-bl)*2^16+ (aa(cc+idx+1)-bl)*2^12 + (aa(cc+idx+2)-bl)*2^6+ aa(cc+idx+3)-bl)*1e1; % back pressure lets make it kPa

idx = 45;
cond= ((aa(cc+idx)-bl)*2^16+ (aa(cc+idx+1)-bl)*2^12 + (aa(cc+idx+2)-bl)*2^6+ aa(cc+idx+3)-bl); % condactivity, not sure about units

idx = 51;
nn= ((aa(cc+idx)-bl)*2^16+ (aa(cc+idx+1)-bl)*2^12 + (aa(cc+idx+2)-bl)*2^6+ aa(cc+idx+3)-bl); % condactivity, not sure about units


% idx = 54;
% nn= ((aa(cc+idx)-bl)*2^16+ (aa(cc+idx+1)-bl)*2^12 + (aa(cc+idx+2)-bl)*2^6+ aa(cc+idx+3)-bl); % condactivity, not sure about units
% nn= ((aa(cc+idx)-bl)); % condactivity, not sure about units

% plot(xconv, yconv, 'r')
% plot(ml, UV, 'b')
% plot(ml, cond, 'k')
% plot(ml,flow_rate, 'm');
% plot(ml,kPa, 'c');
% plot(ml,nn, 'm');

if nargin==2
    if ~isempty(strfind(lower(varargin{1}), 'time'))
        ax.x = time;
        ty = ml;
        ax.xlabel = 'Time, min';
        trace5 = 'Volume, mL';
    end
    if ~isempty(strfind(lower(varargin{1}), 'volume'))
        ax.x = ml;
        ty = time;
        ax.xlabel = 'Volume, mL';
        trace5 = 'Time, min';
    end
    if isempty(ty), error('Wrong x-axis parameter, must be either time or volume'); end
else    
    ty = ml;
end
y(:, 1) = UV;
y(:, 2) = kPa*1e-3;
y(:, 3) = flow_rate;
y(:, 4) = cond;
y(:, 5) = ty;    
y(:, 6) = nn;

ax.type = 'data';
dsc.size_x = numel(ax.x);
dsc.size_y = 6;

dsc.RunType = (RunType);
dsc.Akta    = Akta;
dsc.trace1 = 'UV @ 280mnm';
dsc.trace2 = 'Pressure, MPa';
dsc.trace3 = 'Flow Rate, ml/min';
dsc.trace4 = 'Conductance, ?';
dsc.trace5 = trace5;
dsc.trace6 = 'unknown trace';

ax.title = 'Akta Start';
switch nargout
  case 1,
    varargout = {y};
  case 2,
    varargout = {ax, y};
  case 3,
    varargout = {ax, y, dsc};
end
