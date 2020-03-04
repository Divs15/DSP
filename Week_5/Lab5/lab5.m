%%
% LAB5
%%
% 5.1.1
%%

%{
Background
1.SARS (Severe Acute Respiratory Syndrome) was an illness that first appeared
in the second half of 2002 in Guangdong Province (China). The disease is caused by the SARS coronavirus.
2. Maximum likelihood phylogenetic tree
The SARS virus was recognized early on as a coronavirus,having the same
genes in the same order as other known coronaviruses.
Phylogenetic tree was generated using the CIPRES Gateway.
SARS_Data_Part_1.fasta file was used as input and further converted to
MAFFT (alignment) then formatting was done using Readseq tool further tree building(RAxML) was done. 
From there we take the data to iTOL for visualization
Output from iTOL:
%}
%
figure(1)
imshow('RAxMLpart1.jpg');
title('ML tree for SARS Data part 1');
%coronavirus Palm Civet is most closely related to human SARS.

%%
% 5.2.1
%%
% Created a maximum likelihood tree using the pipeline described in Part I
% and input data as SARS Data part 2.fasta file
% Output from iTOL:

figure(2)
imshow('RAxMLpart2.jpg');
title('ML tree for SARS Data part 2');


% 2. Epidemic likely began in Zhongshan Guangdong on 26th December 2002.
% 3. Importing sequence
[headers , sequences]= fastaread('SARS_Data_Part_2.fasta')
% 4. Creating a distance matrix for the sequence data using the jukes-cantor method
dist = seqpdist(sequences,'Alphabet','NT','indels','pair');
% 5. Using the UPGMA algorithm to create a tree. 
Tree= seqlinkage(dist,'average',headers)
view(Tree)

%%

%{
6. On analysing the tree, we can see that the Zhongshan Guangdong epidemic
is branched along with Guangzhou Guangdong implying that they shared a
common ancestor.
The difference between the RAxML tree and the UPGMA model tree is that  in
UPGMA the rate of mutations is constant over time for all entries in the
tree.Whereas, the RAxML shows sequences changing, mutating and varying length of times not utilizing the basic assumption of molecular hypothesis.
 %}
