% Kinefold - run kinefold on a sequence using remote server at AWS
% Assumes instance already started with kinefold in AMI (see AWS/startinstance)
function r=kinefold(name,seq,varargin)
defaults=struct('seed',[],'duration',1000,'trace',[],'force',[],'ntrials',1,'cert','~/Dropbox/AWS/keypair20130523.pem');
  %  old cert='~/Documents/Certificates/EC2.pem';
args=processargs(defaults,varargin);
if isempty(args.seed)
  args.seed=randi(10000,1,args.ntrials);
end
if args.ntrials~=length(args.seed)
  error('Specified %d trials, but %d seeds\n', args.ntrials, length(args.seed));
end
[s,host]=system('. ~/Dropbox/AWS/awssetup.sh; ~/Dropbox/AWS/getip');
if s~=0
  error('Failed to get AWS IP address: %s',host)
end
if isempty(host)
  error('AWS instance not running - use AWS/startinstance');
end
if host(end)==10
  host=host(1:end-1);
end
host=['ec2-user@',host];
fprintf('AWS host=<%s>\n',host);
allr=[];
r.seq=seq;
r.name=name;
r.trial=[];
for trial=1:args.ntrials
  tmpdir=sprintf('/tmp/kinefold.%d',args.seed(trial));
  mkdir(tmpdir);
  sshcmd=sprintf('ssh -i %s "%s"',args.cert,host);
  fd=fopen([tmpdir,'/job.req'],'w');
  fprintf(fd,'%d\n',args.seed(trial));
  suffixes={'p','e','rnm','rnms','rnml','rnm2','dat'};
  for i=1:length(suffixes)
    fprintf(fd,'job.%s\n',suffixes{i});
  end
  fprintf(fd,'0		# 0=RNA ; 1=DNA\n');
  fprintf(fd,'6.3460741	# helix minimum free energy in kcal/mol: 6.3460741=10kT\n');
  fprintf(fd,'10000000	# NA\n');
  fprintf(fd,'%d		# folding time requested in msec\n',args.duration);
  fprintf(fd,'1		# pseudoknots   1=yes 0=no\n');
  fprintf(fd,'0		# entanglements	1=yes 0=no\n');
  fprintf(fd,'2 3		# simulation type: 1=renaturation; 2 20=cotrans. @ 20msec/nt\n');
  for i=1:size(args.trace,1)
    fprintf(fd,'T %d %d %d\n',args.trace(i,:));
  end
  for i=1:size(args.force,1)
    fprintf(fd,'F %d %d %d\n',args.trace(i,:));
  end
  fprintf(fd,'job\n');
  fprintf(fd,'job.zip\n');
  fprintf(fd,'<SEQNAME>job_%d\n',args.seed(trial));
  fprintf(fd,'<BASE>job\n');
  fprintf(fd,'<SEQUENCE>%s\n',seq);
  fprintf(fd,'<ZIPFILE>job.zip\n');
  fclose(fd);

  fd=fopen([tmpdir,'/job.dat'],'w');
  fprintf(fd,'< job\n');
  fprintf(fd,'%s\n',seq);
  fclose(fd);

  cmd=sprintf('scp -o StrictHostKeyChecking=no -i %s %s/job.req %s/job.dat %s:',args.cert,tmpdir,tmpdir,host);
  fprintf('Executing %s...',cmd);
  [s,result]=system(cmd);
  if s~=0
    fprintf('Failed %s:\n\t%r\n', cmd, result);
    return;
  end
  fprintf('done\n');

  cmd=[sshcmd,' ./kinefold_long_static job.req'];
  fprintf('Executing %s...',cmd);
  tic
  [s,result]=system(cmd);
  if s~=0
    fprintf('Failed %s:\n\t%r\n', cmd, result);
    return;
  end
  fprintf('done\nElapsed time=%.0f seconds\n',toc);

  cmd=sprintf('scp -i %s "%s:job.*" %s',args.cert,host,tmpdir);
  fprintf('Executing %s...',cmd);
  [s,result]=system(cmd);
  if s~=0
    fprintf('Failed %s:\n\t%r\n', cmd, result);
    return;
  end
  fprintf('done\n');

  dir(tmpdir);
  fd=fopen([tmpdir,'/job.e'],'r');  % Read back trace
  C=textscan(fd,'%s %d %d %d\n',5);
  for i=1:length(C{1})
    t.helixes(i)=struct('TF',C{1}{i},'i',C{2}(i),'j',C{3}(i),'k',C{4}(i));
  end
  C=textscan(fd,'%f\n',1);
  t.a=C{1};
  C=textscan(fd,'%d\n',1);
  t.seed=C{1};
  C=textscan(fd,'%f %f %f %f %f %f %f\n');
  t.data=struct('time',C{1},'energy',C{2},'frac',[C{3:7}]);
  for k=2:size(t.data.frac,2)
    % Kinefold offsets values for each trace to allow offset printing; correct for this
    t.data.frac(:,k)=t.data.frac(:,k)-.02*(k-1);
  end
  fclose(fd);
  t.mfe=rnread([tmpdir,'/job.rnml']);
  t.kinetic=rnread([tmpdir,'/job.rnm']);
  t.subopt=rnread([tmpdir,'/job.rnm2']);
  r.trial=[r.trial,t];
end
