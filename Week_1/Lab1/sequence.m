clc;
close all;
clear;
%%
% Retrieve sequence information from GenBank database
getgenbank('nm_000520','toFile', 'nm_000520.txt');
% load in Matlab
s = genbankread('nm_000520.txt');
% Load Sequence only
seq = getgenbank('nm_000520','SequenceOnly',true)
%%
% Format long sequence output for easy viewing
seqdisp(s.Sequence)
% Count nucleotides in sequence
[seq_counts] = basecount(s.Sequence)% Plot density of nucleotides along sequence
figure(1)
seq_density = ntdensity(s.Sequence)
% Count dimers in nucleotide sequence
figure(2)
[Dimers, Percent] = dimercount(s.Sequence, 'chart', 'pie')
% Count 3-mer in nucleotide sequence
trimer = nmercount(s.Sequence,3)
%%
%lab1.3.1
 % returns codon counts for the first reading frame 
 % and plot the results in a heat map
 %First 3 frames
 figure(3)
 r1codons = codoncount(s.Sequence,'frame',1,'figure',true);
 figure(4)
 r2codons = codoncount(s.Sequence,'frame',2,'figure',true);
 figure(5)
 r3codons = codoncount(s.Sequence,'frame',3,'figure',true);
 figure(6)
 %Reversing the sequence to get last 3 frames
 r4codons = codoncount(s.Sequence,'Reverse','true','frame',1,'figure',true);
 figure(7)
 r5codons = codoncount(s.Sequence,'Reverse','true','frame',2,'figure',true);
 figure(8)
 r6codons = codoncount(s.Sequence,'Reverse','true','frame',3,'figure',true);
 
 %%
 %Lab1.3.2
  %ORFs of minimum length 50 in Frame 1 
  orfminlength_50=seqshoworfs(s.Sequence,'Frames',1,'MinimumLength',50);
  %ORFs of minimum length 500 in Frame 1 
  orfminlength_500=seqshoworfs(s.Sequence,'Frames',1,'MinimumLength',500);
  %%
 %Lab 1.3.3
 % Estimate P(stop) from the sequence nm_000520 and determine the threshold given ? = 0.05
 % Find the total of stop codons in each frame 
frame1_Pstop = r1codons.TAA + r1codons.TAG + r1codons.TGA
frame2_Pstop = r2codons.TAA + r2codons.TAG + r2codons.TGA
frame3_Pstop = r3codons.TAA + r3codons.TAG + r3codons.TGA
frame4_Pstop = r4codons.TAA + r4codons.TAG + r4codons.TGA
frame5_Pstop = r5codons.TAA + r5codons.TAG + r5codons.TGA
frame6_Pstop = r6codons.TAA + r6codons.TAG + r6codons.TGA

%sum of all frames above

Frametotal = frame1_Pstop +  frame2_Pstop + frame3_Pstop + frame4_Pstop + frame5_Pstop + frame6_Pstop

%Total number of codons in sequence
% calculate the value in each codon
c1=struct2cell(r1codons);
c2=struct2cell(r2codons);
c3=struct2cell(r3codons);
c4=struct2cell(r4codons);
c5=struct2cell(r5codons);
c6=struct2cell(r6codons);
%total codons in sequence
Total = sum([c1{:}])+ sum([c2{:}]) + sum([c3{:}])+ sum([c4{:}])+ sum([c5{:}]) + sum([c6{:}]) 
%probability of stop codons
Pstop = (Frametotal./Total)

%%
%probability of k non stop 
%Given confident level alpha
alpha= 0.05 ;
%%
% k is the length of the sequence
k = log(0.05) / log(1-(Pstop));
%optimal length
knew = k + 1 + 1
