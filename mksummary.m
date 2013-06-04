% Create summary field for r tracking ribozyme and aptamer formation
function r=mksummary(r)
maxtime=min(arrayfun(@(z) z.data.time(end), r.trial));
for fnum=1:2
  s=struct('time',0:maxtime);
  s.fracapt=zeros(size(s.time));
  s.fracribo=zeros(size(s.time));
  for i=1:length(r.trial)
    fold=r.trial(i).kinetic.fold;
    kapt=[];kribo=[];
    for k=1:length(fold)
      kapt(k)=istheo(fold(k).seq,fold(k).dot);
      kribo(k)=isactive(fold(k).seq,fold(k).dot);
    end
    if fnum==1
      s.frac(i)=struct('time',[fold.time],'apt',kapt,'ribo',kribo);
    else
      s.frac(i)=struct('time',r.trial(i).data.time,'apt',prod(r.trial(i).data.frac(:,4:5),2),'ribo',prod(r.trial(i).data.frac(:,1:3),2));
    end
    for t=1:length(s.time)
      sel=find(s.frac(i).time<=s.time(t),1,'last');
      if isempty(sel)
        apt(t)=0;ribo(t)=0;
      else
        apt(t)=s.frac(i).apt(sel);
        ribo(t)=s.frac(i).ribo(sel);
      end
    end
    s.fracapt=s.fracapt+apt;
    s.fracribo=s.fracribo+ribo;
  end
  s.fracribo=s.fracribo/length(r.trial);
  s.fracapt=s.fracapt/length(r.trial);
  if fnum==1
    r.summarycheck=s;
  else
    r.summaryhelix=s;
  end
end

if isfield(r.trial,'frac')
  fprintf('Removing old trial.frac field\n');
  r.trial=rmfield(r.trial,{'frac'});
end
if isfield(r,'summary')
  fprintf('Removing old summary field\n');
  r=rmfield(r,'summary');
end

