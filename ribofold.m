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
  if length(r.trial) >= args.ntrials
    return
  end
  args.ntrials=args.ntrials-length(r.trial);
end
fprintf('Running kinefold for %d trials\n',args.ntrials);
% Helix labeling (using L2b12):
%             1         2         3         4         5         6         7         8         9         0         1         2         3
%    1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
%    CUUUUCCGUAUAUCUCGCCAGGCUGUCACCGGAUGUGCUUUCCGGUCUGAUGAGUCCGUUGUCCAUACCAGCAUCGUCUUGAUGCCCUUGGCAGGGACGGGACGGAGGACGAAACAGCGUGGUCCAAGUGAUUCCCAAA';
% stems                   33333 111111       111111       2222                                                 2222   33333
% aptamer                                                      XXXXXXbbbAAAA                AAAAbbbXXXXXX    
s3=findhelix(seq,'GCUGU','ACAGC');
s1=findhelix(seq,'ACCGGA','UCCGGU',s3(1)+s3(3),s3(2)-s3(3));
s2=findhelix(seq,'GA/GUCC/','/GGAC/GAA',s1(1)+s1(3),s3(2)-s3(3));
a1=findhelix(seq,'/..../AUACCAG','UUGG[CA]../..../',s1(1),s2(2));
a2=findhelix(seq,'AUA/CCAG/','/UUGG/[CA]..',a1(1)+a1(3),a1(2)-a1(3));
rnew=kinefold(name,seq,'ntrials',1,'trace',[s1;s2;s3;a1;a2],'duration',args.duration,'ntrials',args.ntrials);
for i=1:length(rnew.trial)
  rnew.trial(i).fracribo=prod(rnew.trial(i).data.frac(:,1:3),2);
  rnew.trial(i).fracapt=prod(rnew.trial(i).data.frac(:,4:5),2);
end
if exist('r','var')
  r.trial=[r.trial,rnew.trial];
else
  r=rnew;
end
maxtime=min(arrayfun(@(z) z.data.time(end), r.trial));
r.summary=struct('time',0:maxtime);
r.summary.fracapt=zeros(size(r.summary.time));
r.summary.fracribo=zeros(size(r.summary.time));
for i=1:length(r.trial)
  apt=interp1(r.trial(i).data.time,r.trial(i).fracapt,r.summary.time);
  r.summary.fracapt=r.summary.fracapt+apt;
  ribo=interp1(r.trial(i).data.time,r.trial(i).fracribo,r.summary.time);
  r.summary.fracribo=r.summary.fracribo+ribo;
end
r.summary.fracribo=r.summary.fracribo/length(r.trial);
r.summary.fracapt=r.summary.fracapt/length(r.trial);

r.helixlabels={'Stem1','Stem2','Stem3','Theo1','Theo2'};
save(cachefile,'-struct','r');

function s=findhelix(seq,pt1,pt2,minpos,maxpos)
if nargin<4
  minpos=1;
end
if nargin<5
  maxpos=length(seq);
end
if sum(pt1=='/')>0
  len=diff(find(pt1=='/'))-1;
  start1=find(pt1=='/',1);
else
  len=length(pt1);
  start1=1;
end
if sum(pt2=='/')>0
  len2=diff(find(pt2=='/'))-1;
  start2=find(pt2=='/',1);
else
  len2=length(pt2);
  start2=1;
end
if len~=len2
  error('Length of helixes different\n');
end
maxpos=maxpos-len+1;
f1=regexp(seq,strrep(pt1,'/',''));
f1=f1(f1>=minpos & f1<=maxpos);
f2=regexp(seq,strrep(pt2,'/',''));
f2=f2(f2>=minpos & f2<=maxpos);
f1=f1(f1+len<max(f2));
f2=f2(f2>min(f1)+len);
if length(f1)~=1 || length(f2)~=1
  error('Error finding helix %s  %s in %s\nf1=%s, f2=%s\n',pt1,pt2,seq(minpos:maxpos), sprintf('%d ',f1),sprintf('%d ',f2));
end
s=[f1+start1-1,f2+start2-1+len-1,len];
