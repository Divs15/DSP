%%
%Homework 1
%%
%Question 1
%1A
seq = getgenpept('O29042','PartialSeq', [141 147], 'toFile', 'O29042.txt');
% sequence of amino acids (in single letter representation) from positions 141-147
%accessing the sequence
seq.Sequence
%1B
%Function created 
%{
    function possible_sequence = nucleotide_seq(desired_sequence)
   
    map = revgeneticcode(1, 'Alphabet', 'DNA');

A = desired_sequence

m = length(A);

result=1;

for n = 1:1:m

   i = A(n);

   result=length(map.(i))*result;

end

result
possible_sequence = result;        
    

end
%}


%1C
possible_seq = nucleotide_seq('HARVARD')


%1d
possible_cat = nt2aa_cgbiased('HARVARD')
%{
function possible_cat = nt2aa_cgbiased(desired_sequence)
    

map = revgeneticcode(1, 'Alphabet', 'DNA');
A= desired_sequence


m = length(A);

z=[];

possible_cat=[];

%g='';

result=1;

for n = 1:1:m

   i = A(n);

   b=(map.(i));

   j=length(b);

   for k=1:1:j

       x=length(strfind((char(b(k))),'C'));

       y=length(strfind((char(b(k))),'G'));

       z=[[z] (x+y)];

   end

   [e f]=max(z);

   z=[];

   possible_cat=[[possible_cat] cell2mat(b(f))];
end
len=length(possible_cat)

end
%}
%%
%Question 4
%Part A
%Blasted it online and got these results
% E. Coli bacterial genomes NC 002655 and NC 002695
% Max Score 9.921e+05 	
% Total Score 1.420e+07 	
% Query value 99% 	
% Ident 99% 	
% Accession Number NC_002695.2
% Its not possible to take global alignment of the sequence , because the
% sequences are too long 
% But if segmentation of sequence is done then we can find the global
% alignment

%Part B
g1= genbankread('gene1.gb'); % load *.gb file as a struct 
x1=g1.Sequence
g2 = genbankread('gene2.gb'); % load *.gb file as a struct 7 CDS = hbb.CDS
x2=g2.Sequence



[Score, Alignment] = nwalign(x1,x2,'Alphabet','NT','Showscore',true)
%region of dense dissimilarity
% According to the heat map, the blue region shows dissimilarity between
% the two sequences. 


%Part C
[Score, Alignment] = nwalign(x1,x2,'Alphabet','AA','Showscore',true)
% nucleotide regions of dissimilarity do not yield a high di?erence in the amino
% acid sequence but they yield a slight difference in the region close to
% the line of optimal alignment (yellow region)


%seqdotplot(x1,x2)
