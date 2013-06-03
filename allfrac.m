function [a,apt,ribo]= allfrac(r)
for i=1:length(r.trial)
  f=r.trial(i).data.frac;
  for j=1:size(f,2)
    f(:,j)=f(:,j)-min(f(:,j));
  end
  a(i,:)=mean(f);
  apt(i)=mean(f(:,4).*f(:,5));
  ribo(i)=mean(f(:,1).*f(:,2).*f(:,3));
end
