function s1 = trim(s)
if isempty(s), s1='';
else
  % remove trailing blanks
  [r,c] = find(s ~= ' ' & s ~= 0);
  if isempty(c), s1 = '';
 	else
    s1 = s(:,min(c):max(c));
  end
end