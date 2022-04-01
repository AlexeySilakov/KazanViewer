function par = SpecMandsc(filename)
h = fopen(filename);
if h<0
  disp('Description file was not found.');
  par = [];
  return;
end
olda = '';
par = [];
forbidden = '~!@#$%^&*()./\:';
section = '';
text = 1; prg = 1; ddsc = 1;
while feof(h)<1
  s = trim(fgetl(h));
  if isempty(s), continue; end;
  sect = find(s=='[' | s==']');
  %   this is a section header   
  if size(sect, 2)==2 & sect(1)==1
    section = s(sect(1)+1:sect(2)-1);
%     par = setfield(par, section, 'section');
  else
    switch section
      case 'text'
        par = setfield(par, ['text', num2str(text)], s);
        text = text + 1;
      case 'program'
        par = setfield(par, ['prg', num2str(prg)], s);
        prg = prg + 1;
      case 'dsc'
          ddsc = ddsc + 1;
          par = setfield(par, ['dsc', num2str(ddsc)], s);
      otherwise
        [a,s]=strtok(s, '=');
        a = trim(a);
        a(a=='/' | a=='\' | a==' ' | a=='.' | a==',')='_';
        par = setfield(par, [section,'_',a], s(2:end));
    end  
  end
end
fclose(h);