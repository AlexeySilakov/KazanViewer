% D01READ Read data from .d01/.exp SpecMan data files
%
%   data = d01read(filename,...)
%   [ax,data] = d01read(filename,...);
%   [ax,data,dsc] = d01read(filename,...);
%
%   Reads data and axis information from Bruker
%   files  DSC/DTA  and  PAR/SPC.   Aditional 
%   arguments  are  always  coming  in pairs 
%   field, value. Fields are:
%     'Dataset', [real dataset, imag dataset] 
%   ax contains fields x,y,z  depends on data 
%   dimension and may contain other fields as 
%   well. dsc is a cell array  of description 
%   strings.

% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2003
% Free for non-commercial use. Use the program at your own risk. 
% The authors retains all rights. 
% Contact: epel@mpi-muelheim.mpg.de

% Igor Gromov, ETH-Hoenggerberg, 17.01.03
% bme,   4-dec-03, MPI  

function varargout=d01read(filename, varargin)

[path,name,ext] = fileparts(filename);
fname = fullfile(path,[name,'.d01']);
dscname = fullfile(path,[name,'.exp']);

fid=fopen(char(fname),'r', 'ieee-le');
if fid<1, error(['File ''',fname,''' can not be open for read.']);end

ndim1=fread(fid, 1,'uint32');       % number of headers, re/im etc.
dformat=fread(fid,1,'uint32');     % format:0-double,1-float
if dformat==1
  sformat='float32';
else
  sformat='double';
end

dstrms = {};
ntotal = 1;
for k=1:ndim1
  ndim2       = fread(fid,1,'int32');
  dstrms{end+1}.dim = fread(fid,4,'int32');
  dstrms{end}.dim(ndim2+1:end) = 1;
  dstrms{end}.first = ntotal;
  dstrms{end}.total = fread(fid,1,'int32');
  ntotal = ntotal + dstrms{end}.total;
end
tmpdat=fread(fid,ntotal,sformat);
fclose(fid);
change_ax_y = 0;
switch(ndim1)
    case 0,
        error('No data present');
    case 2,
        spec=tmpdat(dstrms{1}.first:(dstrms{1}.first+dstrms{1}.total-1))+...
            i*tmpdat((dstrms{2}.first:dstrms{2}.first+dstrms{2}.total-1));
        spec=reshape(spec,dstrms{1}.dim');
    case 1,
        spec=reshape(tmpdat,dstrms{1}.dim');
    otherwise,
        %% find if all data have the same dimensions  
        dim = dstrms{1}.dim;
        isthesame = 1;
        is1D      = (sum(dim~=1) == 1);
        for k=2:ndim1
            if sum(dim == dstrms{k}.dim) ~= 4,
                isthesame = false;
                break;
            elseif is1D & sum(dstrms{k}.dim~=1)~=1 
                is1D = false;    
            end
        end
        
        if isthesame & is1D,
            % read all as columns
            xdim = dim(find(dim~=1));
            spec=reshape(tmpdat, xdim, ndim1);
        elseif isthesame&(sum(dim~=1)==2) % 3D data
            xdim = (dim(find(dim(1))));
            spec=reshape(tmpdat, xdim, ndim1*dim(2));
            change_ax_y = 1;
        else
            spec=tmpdat;
        end
end

dsc = SpecMandsc(dscname);
ax = SpecManpar(dsc);
if change_ax_y
    ax.z = ax.y;
    ax.y = kron(ones(ndim1, 1), ax.z); 
end
if ~isfield(ax, 'x') | size(ax.x, 1)~=size(spec, 1)
  ax.x = [1:size(spec, 1)];
end
ax.type = 'data';
if isfield(dsc,'general_freq1')
    ax.freq1 = kvgetvalue(dsc.general_freq1);
end

% assign output depending on number of output arguments
% requested
switch nargout
 case 1,
   varargout = {spec};
 case 2,
   varargout = {ax, spec};
 case 3,
   varargout = {ax, spec, dsc};
end
