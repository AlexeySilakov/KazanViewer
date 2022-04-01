% BRUKERREAD Read data from .dsc/.dta and .par/.spc data files
%
%   data = brukerread(filename)
%   [ax,data] = brukerread(filename);
%   [ax,data,dsc] = brukerread(filename);
%
%   Reads data and axis information from Bruker
%   files DSC/DTA and PAR/SPC. 
%   ax contains fields x,y,z 
%   depends on data dimension and may contain
%   other fields  as well
%   dsc is a cell array of description strings

% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2003
% Free for non-commercial use. Use the program at your own risk. 
% The authors retains all rights. 
% Contact: epel@mpi-muelheim.mpg.de

% bme,   4-dec-03, MPI  

function varargout = brukerread(filename);

[ppath,name,ext] = fileparts(filename);

switch upper(ext)
case {'.DSC', '.DTA'}
    %--------------------------------------------------
    % DTA/DSC file processing (Bruker BES3T format)
    %--------------------------------------------------
    fname = fullfile(ppath,[name,'.dta']);
    dscname = fullfile(ppath,[name,'.dsc']);
case {'.PAR', '.SPC'}
    %--------------------------------------------------
    % SPC/PAR file processing
    %--------------------------------------------------     fname = fullfile(ppath,[name,'.spc']);
    fname = fullfile(ppath,[name,'.spc']);
    dscname = fullfile(ppath,[name,'.par']);
end

dsc = XEPRdsc(dscname);
ax = XEPRpar(dsc);

dims = [size(ax.x, 1), size(ax.y, 1), 1];

Endian = 'ieee-le';

switch upper(ext)
case {'.DSC', '.DTA'}
    if strcmp(safeget(dsc, 'BSEQ', 'BIG'), 'BIG'), Endian = 'ieee-be'; end
    switch safeget(dsc,'IRFMT','D')
    case 'I',Format = 'int32';
    otherwise, Format = 'float64';
    end
    y = getmatrix(fname,dims,1:3,Format,Endian,ax.complex);
case {'.PAR', '.SPC'}
    JSS = str2num(safeget(dsc,'JSS', '2'));
    if JSS >0,
        if str2num(safeget(dsc, 'DOS', '0')), Endian = 'ieee-be'; end
        if ax.complex,
            %       workaround for the problem of complex number for 2D  
            dims(1) = dims(1)*2;
            y = getmatrix(fname,dims,1:3,'int32',Endian,0);     
            sz = size(ax.x, 1);
            y(1:sz, :) = y(1:sz, :) + j * y(sz+1:end, :);
            y = y(1:sz, :);
        else
            % good format indicator is not found yet
            vers = str2num(safeget(dsc, 'VERS', '0'));
            if vers~=769, Format = 'float'; 'int32';
            else, Format = 'float';
            end
            y = getmatrix(fname,dims,1:3,Format,Endian, ax.complex);     
        end;
    end
end

% nonequal spaced axis (Bes3T format)
TYP = {'X'; 'Y'};
for kk = [1,2]
    if isfield(dsc, [TYP{kk},'TYP'])
        axtype = getfield(dsc, [TYP{kk},'TYP']);
        if strcmp(axtype , 'IGD')
            try
                fid = fopen(fullfile(ppath,[name,'.', TYP{kk}, 'GF']), 'r', Endian);
                tmp = fread(fid, length(getfield(ax,lower(TYP{kk}))), 'float64');
                ax = setfield(ax,lower(TYP{kk}), tmp);
                fclose(fid);
            catch, disp('Error during the reading of the axis file');
            end
        end
    end
end

if isempty(ax.x) | size(ax.x, 1)~=size(y, 1)
    ax.x = [1:size(y, 1)].';    
end
if isempty(ax.y) | size(ax.y, 1)~=size(y, 2)
    ax.y = [1:size(y, 2)].';    
end

% assign output depending on number of output arguments
% requested
switch nargout
case 1,
    varargout = {y};
case 2,
    varargout = {ax, y};
case 3,
    varargout = {ax, y, dsc};
end

function out = getmatrix(FileName,Dims,DimOrder,Format,ByteOrder,Complex)

% Format = 'int32';

% open data file
[fid, ErrorMessage] = fopen(FileName,'r',ByteOrder);
error(ErrorMessage);
% calculate expected number of elements and read in
N = ((Complex~=0)+1)*prod(Dims);
[x,effN] = fread(fid,N,Format);
if effN<N
    error('Unable to read all expected data.');
end

% convert to complex
if Complex
    x = x(1:2:end) + i*x(2:2:end);
end

% reshape to matrix and permute dimensions if wanted
out = ipermute(reshape(x(:),Dims(DimOrder)),DimOrder);

% close file
St = fclose(fid);
if St<0, error('Unable to close data file.'); end

return 