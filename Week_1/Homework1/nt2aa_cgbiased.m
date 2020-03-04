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


