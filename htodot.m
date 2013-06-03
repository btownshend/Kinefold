function [dot,n]=htodot(h)
n=h;
uh=unique(abs(h(h~=0)));
while true
  un=unique(abs(n(n~=0)));
  nold=n;
  for i=1:length(un)
    n(nold==un(i))=i;
    n(nold==-un(i))=-i;
  end
  for i=length(uh):-1:2
    %    fprintf('Reducing %d\n', uh(i));
    left=find(h==uh(i),1);
    right=find(h==-uh(i),1,'last');
    curval=n(left);
    if curval>1 && sum(n(left:right)==curval-1) == sum(n(left:right)==-(curval-1))
      n(h==uh(i))=curval-1;
      n(h==-uh(i))=-(curval-1);
      break;
    end
  end
  if all(n==nold)
    break;
  end
end
dot=blanks(length(n));
dot(:)='.';
lsym='([{<';
rsym=')]}>';
over=n>length(lsym)|-n>length(rsym);
if any(over)
  fprintf('Warning: Too many knots (%d)',max(n));
end
dot(n>0 & ~over)=lsym(n(n>0 & ~over));
dot(n<0 & ~over)=rsym(-n(n<0 & ~over));
dot(n>0 & over)='A'+n(n>0 & over)-length(lsym);
dot(n<0 & over)='Z'+n(n<0 & over)+length(rsym);
