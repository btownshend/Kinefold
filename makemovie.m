% Make a movie of the kinetic folding in r.kinetic and store in r.kinetic.movie
function r=makemovie(r)
k=r.kinetic;
fd=fopen('/tmp/movie.rmm','w');
fprintf(fd,'>%s\n','job');
fprintf(fd,'%s\n',k.seq);
maxsize=[0,0];
for i=1:length(k.fold)
  f=k.fold(i);
  fprintf(fd,'%s\n',f.dot);
end
fclose(fd);
cmd='java -jar ~/Dropbox/SynBio/bin/RNAMovies2.04/RNAMovies2.04.jar -nogui -input /tmp/movie.rmm -output /tmp/movie.gif -gif -size 320';
system(cmd);
return;

% OLD Stuff
maxsize=[0,0];
for i=1:length(k.fold)
  f=k.fold(i);
  im=varna(f.seq,f.dot);
  k.fold(i).image=im;
  mim(1:size(im,1),1:size(im,2),:,i)=im;
  mim(size(im,1):end,size(im,2):end,:,i)=255;
end
k.mov=immovie(mim);
r.kinetic=k;
