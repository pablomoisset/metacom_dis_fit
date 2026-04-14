function fullM = fill_extinct(M,e)

fullM = zeros(numel(e), size(M,2)) ;
current_row = 1 ;
for i = 1:numel(e)
     if e(i)==0
          fullM(i,:) = M(current_row,:) ;
          current_row = current_row+1 ;
      end
end

