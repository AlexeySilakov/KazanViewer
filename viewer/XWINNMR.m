% XWINNMR Read data from XWINNMR data files (1D)
%
%   data = XWINNMR(filename)
%   [ax,data] = XWINNMR(filename)
%   [ax,data,dsc] = XWINNMR(filename)
%
%   Reads data and axis information from fid./acqus
%   and 1r./1i./proc spectrometer files. File extension 
%   is insignificant. ax contains fields x,y,z 
%   depends on data dimension and may contain
%   other fields as well
%   dsc is a cell array of description strings

% bme,   4-dec-03, MPI  
% bme,   12-mar-04, MPI  

function varargout = XWINNMR(filename);

[fpath,name,ext] = fileparts(filename);


% Determine the size and type
SizeTD1 = 1;
SizeTD2 = 1;
dsc  = [];
dsc2 = [];

if exist([fpath,'\fid'],'file')
    dsc = XWINNMRpar(fullfile(fpath, 'acqus'));
    SizeTD1 = str2num(dsc.TD);
    DW=1/(2*str2num(dsc.SW_h)); % s
    ax.x = 1E6*DW*[1:SizeTD1/2].';
    ax.xlabel = 'Time, us';
    if str2num(dsc.BYTORDA) > 0.5, byteorder = 'b';
    else, byteorder = 'l';
    end  
    Return = zeros(SizeTD1,SizeTD2);
    sz = SizeTD2*SizeTD1;
    id = fopen(fullfile(fpath, 'fid'), 'r', byteorder);
    [p,count] = fread(id, sz, 'long');
    p = reshape(p,2,SizeTD1/2).';
    data = p(:,1)+j*p(:,2);
    fclose(id);
    ax.reference = str2num(safeget(dsc, 'BF1', '0'))*1E6;
elseif exist([fpath,'\ser'],'file')
    dsc = XWINNMRpar(fullfile(fpath, 'acqus')); 
    SizeTD1 = str2num(dsc.TD)+45*2;
    dsc2 = XWINNMRpar(fullfile(fpath, 'acqu2s')); 
    SizeTD2 = str2num(dsc2.TD);
    DW=1/(2*str2num(dsc.SW_h)); % s
    ax.x = 1E6*DW*[1:SizeTD1/2].';
    ax.xlabel = 'Time, us';
    if str2num(dsc.BYTORDA) > 0.5, byteorder = 'b';
    else, byteorder = 'l';
    end  
    Return = zeros(SizeTD1,SizeTD2);
    sz = SizeTD2*SizeTD1;
    id = fopen(fullfile(fpath, 'ser'), 'r', byteorder);
    if id < 0, error('File not found'); end
    [p,count] = fread(id, sz, 'long');
    fclose(id);
    p = reshape(p,2,sz/2);
    pp = p(1,:)+j*p(2,:);
    data = reshape(pp,SizeTD1/2,SizeTD2);
    if safeget(dsc2,'VDLIST','<>') & exist([fpath,'\vdlist'],'file')
        ax.y = getvdlist(fpath, SizeTD2);
        ax.ylabel = 'Time, s';
    else
        ax.y = [1:SizeTD2]';
        ax.ylabel = 'Points, pts';
    end
    ax.reference = str2num(safeget(dsc, 'BF1', '0'))*1E6;
elseif exist([fpath,'\procs'],'file')
    dsc = XWINNMRpar(fullfile(fpath, 'procs')); 
    SizeTD1 = str2num(dsc.SI);
    
    SW=str2num(dsc.SW_p); % Hz
    SF=str2num(dsc.SF);   % spec freq MHz
    off = str2num(dsc.OFFSET)*SF;
    
    ax.x = off-SW*[0:SizeTD1-1].'/(SizeTD1-1);
    ax.xlabel = 'Frequncy, Hz';
    ax.reference = str2num(safeget(dsc, 'SF', '0'))*1E6;
    
    sz = SizeTD2*SizeTD1;
    if str2num(dsc.BYTORDP) > 0.5, byteorder = 'b';
    else, byteorder = 'l';
    end
    id = fopen(fullfile(fpath, '1r'), 'r', byteorder);
    [Return1,count] = fread(id, sz, 'int32');
    fclose(id);
    
    id = fopen(fullfile(fpath, '1i'), 'r', byteorder);
    [Return2,count] = fread(id, sz, 'int32');
    fclose(id);
    
    % some coefficient ???
    kkk = str2num(safeget(dsc,'NC_proc','0'));
    data = (Return1 + j*Return2) * 2^kkk;
end
ax.type = 'data';
ax.dx = 0;
ax.freq1 = str2num(safeget(dsc, 'BF1', '0'))*1E6;

% try to read the title 
id = fopen(fullfile(fpath, 'title'), 'r'); 
if id >= 0
    ax.title = fgets(id);
    fclose(id);
end


% assign output depending on number of output arguments
% requested
switch nargout
case 1,
    varargout = {data};
case 2,
    varargout = {ax, data};
case 3,
    varargout = {ax, data, dsc};
end

function res = getvdlist(dirpath, nn)
res = [];
id = fopen([dirpath,'\vdlist'], 'r');
if id < 0, error('VDLIST file was not found'); end
while ~feof(id) & length(res) < nn
    ll = fgetl(id);
    res(end+1,1) = kvgetvalue([ll,'s']);
end
fclose(id);