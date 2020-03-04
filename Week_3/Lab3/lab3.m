%% 
%Lab3
%%
%Lab3.1.1
 clc;close all;clear; 
 % read in txt file 
 Mys_Seq = textread('mystery_sequence.txt', '%q'); 
 % extract sequence information 
Seq = strcat(Mys_Seq{4:11}); 
% write the sequence into a fasta file 
fastawrite('MysSeq.fasta', 'Sequence', Seq); 
% read in the fasta file 
Seq = fastaread('MysSeq.fasta');
 % Perform a blastn and retrieve the report from NCBI 
 %[Data_blastn, RTOE_n] = blastncbi('MysSeq.fasta', 'blastn', 'database', 'nr') 
 %blast_n = getblast(Data_blastn,'ToFile','blast_n.rpt')


%%
%Lab3.2.1
%BLAST mystery_sequence.txt online
%%
%Lab 3.2.2
%had to Run BLAST online...was taking a lot of time to run BLAST in matlab
% Answers provided in comments
 %Perform a blastn and retrieve the report from NCBI
 % read in txt file 
 Mys_Seq = textread('dinoDNA.txt', '%q'); 
 % extract sequence information 
Seq = strcat(Mys_Seq{4:11}); 
% write the sequence into a fasta file 
fastawrite('dinoDNA.fasta', 'Sequence', Seq); 
% read in the fasta file 
Seq = fastaread('dinoDNA.fasta');
 %[Data_blastn, RTOE_n] = blastncbi('dinoDNA.fasta', 'blastn', 'database', 'nr') 
  %blast_n = getblast(Data_blastn,'ToFile','blast_n.rpt')


  %Lab3.2.3
   % Perform a tblastx for the first 40 base pairs 
   %Seq_40 = Seq(1:40); 
   %fastawrite('Seq_40.fasta', 'Sequence', Seq_40); 
   %Seq_40 = fastaread('Seq_40.fasta'); 
   %[Data_tblastx_40, RTOE_40] = blastncbi('Seq_40.fasta', 'tblastx', 'database', 'nr') 
   % Perform a tblastx for the first 300 base pairs 
   %Seq_300 = Seq(1:300); 
   
% 1. What is the difference between blastn and tblastx
%    blastn is a nucleotide to nucleotide function that finds regions of similarity between
%    nucleotide sequences whereas, the purpose of tblastx is to find very distant relationships between nucleotide sequences It'll use one 


% 2. Why were shorter lengths chosen?
%    BLAST is optimal when ran with sequences broken down into smaller % segments

% 3. For the full length BLAST, what organisms were closely related? What taxa do these belong to?
%    Niveispirillum cyanobacteriorum strain TH16 chromosome eg_1, complete sequence was the entry that first came up. The genus of this organism is Niveispirillum and family would be Rhodospirillaceae
% 4. Did this sequence share homology (similarity) with any known 
%    and functionally annotated genes?
%    This sequence shared similarity with Niveispirillum cyanobacteriorum
%    within two  regions: 3142583 to 3143083 and 3142824 to 3143039

%5 What was the di?erence between the tblastx searches performed by the full-length sequence vs. 300 bp sequence vs. 40 bp sequence?
%  Both results returned with different matches. 
%  For 300 bp, Indioceanicola profundi strain SCSIO 08040 chromosome, complete genome



   %fastawrite('Seq_300.fasta', 'Sequence', Seq_300); 
   %Seq_300 = fastaread('Seq_300.fasta'); 
   %[Data_tblastx_300, RTOE_300] = blastncbi('Seq_300.fasta', 'tblastx', 'database', 'nr')

  
  %Results obtained online 
  %1  What are the top hits of this BLAST search?
  %a  Gallus gallus GATA binding protein 1 (globin transcription factor 1) (GATA1), mRNA
  %b  X.laevis GATA-binding protein (XGATA-2) gene, complete cds
  %c  PREDICTED: Xenopus laevis GATA binding protein 2 L homeolog (gata2.L), transcript variant X3, mRNA 
  %d  Xenopus laevis GATA binding protein 2 L homeolog (gata2.L), mRNA 
  
  %2. What organism do you think Mark used for his dinoDNA sequence. 
  %   So result a, 'Gallus Gallus' is some breed of chicken. It is known that birds are recent relatives to the dinosuars. The second result is a breed
  %   of frogs, which was noted was mixed in the extracted dinosaur DNA.
  %3. Are the top hits for the blastx search the same as the ones you saw for blastn? Why might this be the case?
  %   The top hit was gallus gallus again, and the second hit was a zebra fish.
  %   The third hit was the breed of frog.
  %   The dinoDNA must be containing similar amino acid sequences from these three animals.
  %4. What is the hidden message that Mark put in the sequence?
  %   Mark
  

  %%
  %Lab3.3.1
   %for NT
   
   gene1 = getgenpept('AAD01939');
   gene2 = getgenpept('AAQ67266'); 
   % Sequence Dot Plot 5
   seqdotplot(gene1,gene2) 
   axis square 
    title('NT') 
    % nw alignment with the scoring sapce and wining path 
     
    % for AA
     
    % Sequence Dot Plot 
    seqdotplot(gene1.Sequence,gene2.Sequence) 
    axis square 
    title('AA') 
    % nw alignment with the scoring space and wining path 
    [score, alignment] = nwalign(gene1, gene2, 'Alphabet', 'AA', 'showscore', true)
