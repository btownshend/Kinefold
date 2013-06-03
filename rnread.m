% Read a .rnml or other folding files output from kinefold
function r=rnread(filename)
fd=fopen(filename,'r');
if fd<0
  error('Unable to open file: %s\n', filename);
end
line=fgetl(fd);
r.name=line(4:end);
line=fgetl(fd);
r.seq=line;
while true
  line=fgetl(fd);
  hn=fgetl(fd);
  if isnumeric(line) || isnumeric(hn)
    break;
  end
  delim=find(hn=='H',1);
  
  d1=line(delim:end);
  line=line(1:delim-1);
  d2=hn(delim:end);
  hn=hn(1:delim-1);
  lb=find(line(1:2:end)=='[' | line(1:2:end)=='^');
  rb=find(line(1:2:end)==']' | line(1:2:end)=='^')-1;
  assert(length(lb)==length(rb));
  hvals=zeros(1,floor(length(line)/2));
  for i=1:length(lb)
    id=strrep(strrep(hn(lb(i)*2-1:rb(i)*2),'-',''),' ','');
    hnum=str2num(id);
    if id(end)==''''
      hnum=-hnum;
    end
    %fprintf('%d-%d #%d\n', lb(i), rb(i), hnum);
    hvals(lb(i):rb(i))=hnum;
  end
  f=struct('hvals',hvals);
  f.dot=htodot(hvals);
  C=textscan(d1,' | Structure %d, Free-energy : %f kcal/mol #',1);
  if length(C)~=2 || isempty(C{1}) || isempty(C{2})
    C=textscan(d1,' | %f kcal/mol reached after %f ms, 3''-end : %d over %d bases #');
    if isempty(C{1})||isempty(C{2}) || isempty(C{3}) || isempty(C{4})
      error('Error scanning line: %s\n', d1);
    end
    f.energy=C{1};
    f.time=C{2};
    f.length=C{3};
    f.total=C{4};
  else
    f.structnum=C{1};
    f.energy=C{2};
  end
  f.seq=r.seq(1:length(hvals));
  if ~isfield(r,'fold')
    r.fold=f;
  else
    r.fold=[r.fold,f];
  end
end
