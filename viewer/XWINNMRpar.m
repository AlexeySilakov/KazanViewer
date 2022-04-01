% reading Bruker nmr header of .acq type

function par = XWINNMRpar(filename)
h = fopen(filename);
olda = '';
par = [];
while feof(h)<1
  s = fgetl(h);
  if strfind(s,'##')
    if strfind(s,'$$'), continue; end; 
    s(s=='#') = ' ';
    [fld, s]=strtok(trim(s),'=');
    fld = fld(fld~='$');
    if ~isempty(fld)
      par = setfield(par, fld, trim(s(2:end)));
    end
  else
    if ~isempty(fld)
      par = setfield(par, fld, [getfield(par,fld),'; ',trim(s)]);
    end
  end
end
fclose(h);

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