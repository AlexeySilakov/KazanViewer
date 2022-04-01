function varargout = kvscriptload(fname, varargin)

fid = fopen(fname, 'r');
if fid == -1, error('File not found'); return; end
str = {};
while feof(fid) < 1,
  st = fgetl(fid);
  if ~isempty(st) & st(1) ~= '%',
    eval([st, ';'])
      str{end+1} = st;
  end
end
fclose(fid);

for k=1:min(nargin-1, nargout)
   try, eval(['varargout{k}=',varargin{k},';']); 
    catch disp(lasterr);
    end
end

