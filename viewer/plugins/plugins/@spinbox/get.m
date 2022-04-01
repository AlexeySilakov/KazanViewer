function out = get(hands, varargin)
% Get function for 'spinbox' object
if ~isa(hands,'spinbox')
   error('Cannot get properties: not a spinbox object!');
end

prop = get(hands.push2, 'userdata');
err = 0;
if nargin == 1
    out = prop;
elseif nargin == 2
    param = varargin{1};
    if ischar(param)
        try 
            out = getfield(prop, param);
        catch
            err = 1;
        end
    elseif iscell(param)
        for ci = 1:length(param)
            try 
                out{ci} = getfield(prop, param{ci});
            catch
                err = 1;
                break;
            end
        end    
    else
        error('get: Invalid property name.');
    end
    if err
        error('get: Invalid property name.');
    end
else
    error('get: Too many input arguments.');
end