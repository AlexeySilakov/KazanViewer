% function scriptload(fname, varargin)
% typical usage: [Exp, Sys, Opt]=scriptload('<fname>','Exp','Sys','Opt') 

function varargout = scriptload(fname, varargin)

fid = fopen(fname, 'r');
if fid == -1, error(sprintf('File ''%s''not found', fname)); return; end
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

