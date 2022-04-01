function varargout = ssrlexafs_ave(filename, varargin)
fid=fopen(filename,'r', 'ieee-le');

try
    while~feof(fid)
        num = fread(fid, inf, 'int32');
%         dat = fread(fid, num/4, 'real*4');
    end
end
fclose(fid);