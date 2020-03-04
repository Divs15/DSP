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