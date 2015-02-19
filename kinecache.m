% Fold a riboswitch, tracing the stem loops and aptamers
function r=kinecache(name, seq,varargin)
defaults=struct('seed',[],'duration',1000,'ntrials',1,'cachedir','~/Dropbox/Synbio/Kinetics/KINEFOLD.CACHE','maxlen',600);
args=processargs(defaults,varargin);

seq=upper(seq);
seq(seq=='T')='U';
seq=seq(seq~=' ');
if length(seq)>args.maxlen
  fprintf('Truncating sequence to first %d/%d nt to allow kinefold to run!\n',args.maxlen,length(seq));
  seq=seq(1:args.maxlen);
end
  
cachefile=sprintf('%s/%s-%d-%d.mat',args.cachedir,name,length(seq),seq2hash(seq));
if exist(cachefile,'file')
  r=load(cachefile);
  fprintf('Loaded %d trials from cache file %s\n', length(r.trial), cachefile);
  args.ntrials=max(0,args.ntrials-length(r.trial));
end
if args.ntrials>0
  fprintf('Running kinefold for %d trials\n',args.ntrials);
  rnew=kinefold(name,seq,'ntrials',1,'duration',args.duration,'ntrials',args.ntrials);
  if exist('r','var')
    r.trial=[r.trial,rnew.trial];
  else
    r=rnew;
  end
end

r.name=name;
save(cachefile,'-struct','r');

