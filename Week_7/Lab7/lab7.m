%Lab7
%%

x=getgenbank('NC_001416','SequenceOnly',true)
figure(1)
ntdensity(x,'window',2000)

figure(2)
ntdensity(x,'window',3000)
figure(2)
ntdensity(x,'window',4000)
hold on
%2
Transition=rand(2,2)
Emission=rand(2,4)
%normalize them
Transition=Transition./sum(Transition,2)
Emission=Emission./sum(Emission,2)


seq=nt2int(x)
[esT,estE]=hmmtrain(seq,Transition,Emission,'Maxiterations',5000)
Estimated_states=hmmviterbi(seq,esT,estE)
plot(Estimated_states)

