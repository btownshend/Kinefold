% Fold a riboswitch, tracing the stem loops and aptamers
function r=ribofold(name, seq,varargin)
defaults=struct('seed',[],'duration',1000,'ntrials',1,'cachedir','~/Dropbox/Synbio/Kinetics/KINEFOLD.CACHE');
args=processargs(defaults,varargin);

seq=upper(seq);
seq(seq=='T')='U';
seq=seq(seq~=' ');

cachefile=sprintf('%s/%s.mat',args.cachedir,seq);
if exist(cachefile,'file')
  r=load(cachefile);
  fprintf('Loaded %d trials from cache file %s\n', length(r.trial), cachefile);
  args.ntrials=max(0,args.ntrials-length(r.trial));
end
if args.ntrials>0
  fprintf('Running kinefold for %d trials\n',args.ntrials);
  % Helix labeling (using L2b12):
  %             1         2         3         4         5         6         7         8         9         0         1         2         3
  %    1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
  %    CUUUUCCGUAUAUCUCGCCAGGCUGUCACCGGAUGUGCUUUCCGGUCUGAUGAGUCCGUUGUCCAUACCAGCAUCGUCUUGAUGCCCUUGGCAGGGACGGGACGGAGGACGAAACAGCGUGGUCCAAGUGAUUCCCAAA';
  % stems                   33333 111111       111111       2222                                                 2222   33333
  %                                                                ((((...((.((((((...))))))....))...))))
  % aptamer                                                        4444bbb   555555   555555      bbb4444
  s3=findhelix(seq,'GCUGU','ACAGC');
  s1=findhelix(seq,'ACCGGA','UCCGGU',s3(1)+s3(3),s3(2)-s3(3));
  s2=findhelix(seq,'GA(GUCC)','(GGAC)GAA',s1(1)+s1(3),s3(2)-s3(3));
  if s2(2)-s2(1) > s1(2)-s1(1)
    aloc=[s2(1)+s2(3),s2(2)-s2(3)];
    fprintf('Aptamer is in stem2');
  else
    aloc=[s1(1)+s1(3),s1(2)-s1(3)];
    fprintf('Aptamer is in stem1');
  end
  fprintf(' with maximum length %d: %s\n',aloc(2)-aloc(1)+1,seq(aloc(1):aloc(2)));
  a1=findhelix(seq,'(....)AUACCAG','UUGG[CA]..(....)',aloc(1)-4,aloc(2)+4)
  if isempty(a1)
    a2=[];
  else
    a2=findhelix(seq,'GCAUC','GAUGC',a1(1)+a1(3),a1(2)-a1(3));
    if isempty(a2)
      a2=findhelix(seq,'AUA.*(...)GUCUU','GUCUU(...).*[CA]AG',a1(1)+a1(3),a1(2)-a1(3));
      if isempty(a2)
        a2=findhelix(seq,'AUA.*(...)AGUC','AGUC(...).*[CA]AG',a1(1)+a1(3),a1(2)-a1(3));
        if isempty(a2)
          a2=findhelix(seq,'AUA.*(...)UCU','UCU(...).*[CA]AG',a1(1)+a1(3),a1(2)-a1(3));
          if isempty(a2)
            error('Unable to locate inner helix of aptamer\n%s[%s]%s[%s]%s',seq(1:a1(1)-1),seq(a1(1):a1(1)+a1(3)-1),seq(a1(1)+a1(3):a1(2)-a1(3)),seq(a1(2)-a1(3)+1:a1(2)),seq(a1(2)+1:end));
          end
        end
      end
    end
  end
  rnew=kinefold(name,seq,'ntrials',1,'trace',[s1;s2;s3;a1;a2],'duration',args.duration,'ntrials',args.ntrials);
  if exist('r','var')
    r.trial=[r.trial,rnew.trial];
  else
    r=rnew;
  end
end

r.name=name;
r.helixlabels={'Stem1','Stem2','Stem3','Theo1','Theo2'};
r=mksummary(r);
save(cachefile,'-struct','r');

function s=findhelix(seq,pt1,pt2,minpos,maxpos)
if nargin<4
  minpos=1;
end
if nargin<5
  maxpos=length(seq);
end
if sum(pt1=='(' | pt1==')')>0
  assert(sum(pt1=='(')==1 && sum(pt1==')')==1);
else
  pt1=['(',pt1,')'];
end
if sum(pt2=='(' | pt2==')')>0
  assert(sum(pt2=='(')==1 && sum(pt2==')')==1);
else
  pt2=['(',pt2,')'];
end
% Find the tokens matching the () part
t1={};t2={};
for i=1:length(seq)-1
  t1p=regexp(seq(i:end),['^',pt1],'tokenExtents');
  if ~isempty(t1p)
    t1{end+1}=t1p{1}+i-1;
  end
  t2p=regexp(seq(i:end),['^',pt2],'tokenExtents');
  if ~isempty(t2p)
    t2{end+1}=t2p{1}+i-1;
  end
end
f1=[];f2=[];len=[];
for i=1:length(t1)
  if (t1{i}(1)<minpos) || (t1{i}(2)>maxpos)
    fprintf('t1{i}=[%d,%d], out of bounds\n', t1{i});
    continue;
  end
  for j=1:length(t2)
    if t2{j}(1)<minpos || t2{j}(2)>maxpos || t2{j}(1)<t1{i}(2);
      fprintf('t2{j}=[%d,%d], out of bounds\n', t2{j});
      continue;
    end
    if diff(t2{j})~=diff(t1{i})
      fprintf('t1{i}=[%d,%d], t2{j}=[%d,%d], different lengths\n', t1{i}, t2{j});
      continue;
    end
    bad=false;
    for k=0:diff(t1{i})
      if ~isbasepair(seq(t1{i}(1)+k),seq(t2{j}(2)-k))
        fprintf('Not a base pairing at position %d of %s/%s\n', k, seq(t1{i}(1):t1{i}(2)),seq(t2{j}(1):t2{j}(2)));
        bad=true;
        break;
      end
    end
    if bad
      continue;
    end
    if any(f1==t1{i}(1) & f2==t2{j}(2) & len==diff(t1{i})+1)
      fprintf('Already have %d,%d\n', t1{i}(1), t2{j}(2));
      continue;
    end
    f1=[f1,t1{i}(1)];
    f2=[f2,t2{j}(2)];
    len=[len,diff(t1{i})+1];
  end
end
if length(f1)<1
  fprintf('Unable to find helix %s  %s in %s\nf1=%s, f2=%s\n',pt1,pt2,seq(minpos:maxpos), sprintf('%d ',f1),sprintf('%d ',f2));
  s=[];
elseif length(f1)>1
  fprintf('Multiple solutions to find helix %s  %s in %s\nf1=%s, f2=%s\n',pt1,pt2,seq(minpos:maxpos), sprintf('%d ',f1),sprintf('%d ',f2));
  s=[f1(1),f2(1),len(1)];
else
  fprintf('Found pt1 %s at %s, pt2 %s at %s:', pt1, sprintf('%d ',f1),pt2, sprintf('%d ',f2));
  fprintf(' helix: %s/%s\n', seq(f1:f1+len-1), seq(f2-len+1:f2));
  s=[f1,f2,len];
end

function y=isbasepair(a,b)
if a=='G'
  y=(b=='C' || b=='U');
elseif a=='C'
  y=b=='G';
elseif a=='A'
  y=b=='U';
else
  y=(b=='G' || b=='A');
end