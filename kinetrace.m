function kinetrace(r,trialnum)
clf;
subplot(411);
trial=r.trial(trialnum);
plot(trial.data.time,trial.data.frac);
l={};
for i=1:size(trial.data.frac,2)
  h=trial.helixes(i);
  l{i}=sprintf('%s %c%d,%d,%d',r.helixlabels{i},h.TF,h.i,h.j,h.k);
end
legend(l);
xlabel('Time (msec)');
ylabel('Fraction');

subplot(412);
plot(trial.data.time,trial.data.energy);
xlabel('Time (msec)');
ylabel('Energy (kCal/mole)');

subplot(413);
s=r.summaryhelix;
st=s.frac(trialnum);
plot(st.time,st.ribo,'g');
hold on;
plot(st.time,st.apt,'r');
c=axis;axis([c(1),c(2),-0.1,1.1]);
legend({'Ribozyme Formed','Aptamer Formed'});
xlabel('Time(msec)');
title('Ribozyme/Aptamer Formation');

subplot(414);
plot(s.time,s.fracribo,'g');
hold on;
plot(s.time,s.fracapt,'r');
c=axis;axis([c(1),c(2),-0.1,1.1]);
legend({'Ribozyme Formed','Aptamer Formed'});
xlabel('Time(msec)');
title(sprintf('Ribozyme/Aptamer Formation over %d trials',length(r.trial)));

h=suptitle(sprintf('%s trial %d seed=%d', r.name, trialnum,trial.seed));
set(h,'Interpreter','none');

