function highest_content = cg_value(aa_cell)
 cg=0;
 position=0;
 
 for n=1:length(aa_cell)
     c=count(aa_cell(n),{'G','C'})
     
     if(c>=cg)
         cg=c
         position=n
     end
 end
 highest_content = aa_cell(position);
end

         