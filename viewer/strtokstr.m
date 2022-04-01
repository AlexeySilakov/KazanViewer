function varargout = strtokstr(varargin)
%  STRTOKSTR Find token in string.
%     STRTOKSTR(S) returns the first token in the string S delimited
%     by "white space".   Any leading white space characters are ignored.
%  
%     STRTOKSTR(S,D) returns the first token delimited by D. 
%     In contrast to "strtok", D can be string with several simbols
%  
%     [T,R] = STRTOKSTR(...) also returns the remainder of the original
%     string, without D. 
%     So, original string can be constracted as [T, D, R].
%
%     If the token is not found in S then R is an empty string and T
%     is same as S. 
%   
%     See also STRTOK, STRFIND

if isempty(varargin),
error('Not enough input arguments.');
return
end

switch nargin
    case 1,
        str = varargin{1};
        delim = ' ';
    case 2
        str = varargin{1};
        delim = varargin{2};
    otherwise
        str = varargin{1};
        delim = ' ';
end

if length(str) < length(delim)
    warning('length(str) < lenght(delim)');
    str1 = str;
    str2 = '';
else
    nnn = strfind(str, delim);
    if ~isempty(nnn), 
        num = nnn(1);
        if ~isempty(num),
            str1 = str(1:num-1);
            len = length(delim);
            str2 = str((num+len):end);
        else
            str1 = str;
            str2 = '';
        end
    else
        str1 = '';
        str2 = str;
    end
end

switch nargout
    case {0, 1}
        varargout{1} = str1;
    case 2
        varargout{1} = str1;
        varargout{2} = str2;
    otherwise
        error('Too many output arguments.');
        return
end
        
    