function set(hands, varargin)
% SET function for 'spinbox'
if ~isa(hands,'spinbox')
   error('Cannot get properties: not a spinbox object!');
end

spinbox(hands, varargin{1}, varargin{2});