% Test with A_L2b12_S
function r=l2b12(ntrials)
seq='CUUUUCCGUAUAUCUCGCCAGGCUGUCACCGGAUGUGCUUUCCGGUCUGAUGAGUCCGUUGUCCAUACCAGCAUCGUCUUGAUGCCCUUGGCAGGGACGGGACGGAGGACGAAACAGCGUGGUCCAAGUGAUUCCCAAA';
% stems                   33333 111111       111111       2222                                                 2222   33333
% aptamer                                                      XXXXXXbbbAAAA                AAAAbbbXXXXXX    
%             1         2         3         4         5         6         7         8         9         0         1         2         3
%    1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
r=ribofold('A_L2b12_S',seq,'ntrials',ntrials);
