% RMNREAD Read data from rmn files
%
%   data = rmnread(filename, ...)
%   [ax,data] = rmnread(filename, ...);
%   [ax,data,dsc] = rmnread(filename, ...);
%
%   Reads data and axis information from Bruker
%   files DSC/DTA and PAR/SPC. 
%   Aditional arguments are always coming in pairs
%   field, value. Possible fields are: 
%   'type' : '1d' / '2d' / '1dtxt'/ '2dtxt', '1d' by default
%     for '1dtxt' type additional parameters are: 
%     'ncomplex' (default sizeofdata/2) and 
%     'inittime' (0) and 'dwelltime' (20E-6)
%   ax contains fields x,y,z depends on data 
%   dimension and may contain other fields 
%   as well. dsc is a cell array of description 
%   strings. 

% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2003
% Free for non-commercial use. Use the program at your own risk. 
% The authors retains all rights. 
% Contact: epel@mpi-muelheim.mpg.de

% bme,   4-dec-03, MPI  

function varargout = rmnread(filename, varargin)

opt.version = 0;
% create argument structure
if nargin > 1 | mod(nargin,2) == 1
  for ii=2:2:nargin
    opt = setfield(opt, varargin{ii-1}, varargin{ii});
  end
end

if ~isfield(opt, 'type')
  opt.type = '1d';
end

switch(lower(opt.type))
  case '1d'
    % 1D data
    fid = fopen(filename, 'r','ieee-be');
    if fid == 0
      return
    end
    version = fread(fid, 1, 'char');
    opt.ncomplex = num2str(version);
    ax.ncomplex = fread(fid, 1,'int');
    opt.ncomplex = num2str(ax.ncomplex);
    ax.dwelltime = fread(fid, 1,'double');
    opt.dwelltime = num2str(ax.dwelltime);
    ax.inittime = fread(fid, 1,'double');
    opt.inittime = num2str(ax.inittime);
    ax.freq1 = fread(fid, 1,'double');
    opt.freq = num2str(ax.freq1);
    ax.freqoffset = fread(fid, 1,'double');
    opt.freqoffset = num2str(ax.freqoffset);
    opt.comment = fread(fid, 2*512,'uchar');
    
    p = fread(fid, ax.ncomplex*2, 'float');
    k = reshape(p, 2, ax.ncomplex);
    y = (k(1, :)+i*k(2, :)).';
    if nargout > 1 
      ax.x  = ax.inittime + ax.dwelltime*[1:ax.ncomplex].';
    end
    ax.xlabel = '?';
    ax.type = 'data';
    fclose(fid);
    
  case '2d'
    fid = fopen(filename, 'rb','ieee-be');
    if fid == 0, return; end
    opt.version = fread(fid, 1, 'char');
    ax.n2d = fread(fid, 1,'int');
    opt.n2d = num2str(ax.n2d);
    ax.dwell2d = fread(fid, 1,'double');
    opt.dwell2d = num2str(ax.dwell2d);
    ax.init2d = fread(fid, 1,'double');
    ax.freq2d = fread(fid, 1,'double');
    ax.offset2d = fread(fid, 1,'double');
    ax.n1d = fread(fid, 1,'int');
    opt.n1d = num2str(ax.n1d);
    ax.dwell1d = fread(fid, 1,'double');
    opt.dwell1d = num2str(ax.dwell1d);
    ax.init1d = fread(fid, 1,'double');
    ax.freq1d = fread(fid, 1,'double');
    ax.offset1d = fread(fid, 1,'double');
    comment = fread(fid, 2*512,'char');
    pos = min(find(~comment));
    opt.comment = char(comment(1:pos)');
    opt.subtype(end+1:4)='T';
    
    if opt.version ~= 4
        error('Not saupported rmn file version.')
    end
    
    if opt.subtype(3)=='T'
        ax.x = (ax.init1d + ax.dwell1d*[0:ax.n1d].')*1E3;
        ax.xlabel = 'T2, ms';
    else
        ax.x = 1/ax.dwell1d/2/ax.n1d*[-ax.n1d:2:ax.n1d]'*0.001;
        ax.xlabel = 'F2, kHz';
    end
    if opt.subtype(3)=='T'
        ax.y = (ax.init2d + ax.dwell2d*[0:ax.n2d].')*1E3;
        ax.xlabel = 'T1, ms';
    else
        ax.y = 1/ax.dwell2d/2/ax.n2d*[-ax.n2d:2:ax.n2d]'*0.001;
        ax.ylabel = 'F1, kHz';
    end    
    n = (ax.n2d+1)*(ax.n1d+1);
    p = fread(fid, n*2, 'float');
    fclose(fid);
    p = reshape(p, 2, n);
    y = reshape(p(1,:)+j*p(2,:), ax.n1d+1, ax.n2d+1);
    
  case '1dtxt'
    p = load(filename);

    data.inittime = safeget(opt, 'inittime', 0.);
    data.dwelltime = safeget(opt, 'dwelltime',10E-6);
    data.ncomplex = safeget(opt, 'ncomplex', size(p, 1)/2);
    
    k = reshape(p, 2, data.ncomplex);
    y = (k(1, :)+i*k(2, :)).';
    ax.x    = data.inittime + data.dwelltime*[1:data.ncomplex].';
  case '2dtxt'
    data.cfreq = safeget(opt, 'cfreq', 0.);
    data.xinit = safeget(opt, 'xinit', 0.);
    data.xsize = safeget(opt, 'xsize', 1.);
    data.xdwell = safeget(opt, 'xdwell', 10E-6);
    data.yinit = safeget(opt, 'yinit', 0.);
    data.ysize = safeget(opt, 'ysize', 1.);
    data.ydwell = safeget(opt, 'ydwell', 10E-6);
    data.complex = safeget(opt, 'complex', 1) > 0;

    p = load(filename);
    sz = length(p);
    
    if data.complex
        k = reshape(p, 2, sz/2);
        y = (k(1, :)+i*k(2, :)).';
    else
        y = p;
    end
    sz = length(y);
    if data.xsize*data.ysize > sz
        disp(sprintf('Specified sizes %dx%d=%d not equal to array size %d',data.xsize, data.ysize, data.xsize*data.ysize, sz));
        error('!');
        return
    elseif data.xsize*data.ysize < sz
        warning(sprintf('Specified sizes %dx%d=%d not equal to array size %d',data.xsize, data.ysize, data.xsize*data.ysize, sz));
        y = y(1:data.xsize*data.ysize);
    end
    
    y = reshape(y,data.xsize,data.ysize);
    if mod(data.xsize,2) & mod(data.ysize,2)
        ax.x    = 1./data.xdwell/(data.xsize-1)*[-(data.xsize-1)/2:(data.xsize-1)/2].';
        ax.y    = 1./data.ydwell/(data.ysize-1)*[-(data.ysize-1)/2:(data.ysize-1)/2].';
        ax.reference = data.cfreq;
    else
        ax.x    = data.xinit + data.xdwell*[1:data.xsize].';
        ax.y    = data.yinit + data.ydwell*[1:data.ysize].';
    end
    ax.reference = data.cfreq;
end

% assign output depending on number of output arguments
% requested
switch nargout
  case 1,
    varargout = {y};
  case 2,
    varargout = {ax, y};
  case 3,
    varargout = {ax, y, opt};
end

return