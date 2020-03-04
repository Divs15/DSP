%LAB8
%%
% Lab 8.1.1

 
%1.
 hbb = genbankread('hbb_region_chr11(2).gb'); % load hbb_region_chr11.gb.gb file
 % extract CDS field (Coding Sequence from hbb)
 CDS = hbb.CDS;
 % get the range indices of CDS
 CDSrange1 = CDS.indices;
 CDSrange2 = hbb.CDS.indices; % a quicker way to get the range
 % extract the first base from the given sequence
 A=hbb.Sequence(1); 
novds=size(CDS,2)


%Length of CDS is 1422
cds1length=CDS(1).indices(end)-CDS(1).indices(1)+1

%Getting coding and non-coding regions in the sequence
CodingRegion=[]; %coding
for i=1:2:length(CDS(5).indices)
    CodingRegion=[CodingRegion CDS(5).indices(i):CDS(5).indices(i+1)];
end
NonCoding=setdiff(CDS(5).indices(1):CDS(5).indices(end),CodingRegion);%non coding
%2.
A=1.5;
C=0.5;
G=-0.5;
T=-1.5;
coding = hbb.Sequence(CodingRegion);
noncoding=hbb.Sequence(NonCoding);
mapcoding=mapseq(coding);
mapnoncoding=mapseq(noncoding);

%%
% Lab 8.2.1
%1. What is H(z)
% H(z) is the transfer function given by:
% -a(2)*Z^-1-a(3)*Z^-2-a(4)*Z^-3.... a(p+1)*Z^-p
ARcoding=lpc(mapcoding,100);
ARnoncoding=lpc(mapnoncoding,100);

%%
% Lab 8.3.1
est_coding = filter([0 -ARcoding(2:end)],1,mapcoding);
est_noncoding=filter([0 -ARnoncoding(2:end)],1,mapnoncoding);
figure
plot(mapcoding(200:300));
hold on
plot(est_coding(200:300));
title('coding')
legend('original','predicted')
figure
plot(mapnoncoding(200:300));
hold on
plot(est_noncoding(200:300));
title('noncoding')
legend('original','predicted')

format
MSEcoding=(sum((abs(mapcoding-est_coding)).^2))/length(mapcoding);
MSEnoncoding=(sum((abs(mapnoncoding-est_noncoding)).^2))/length(mapnoncoding);
