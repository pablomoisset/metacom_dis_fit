bmp = bmp_orig ;

for i=1:6
    test_bmp = bmp{i};
    test_bmp(16,7)=(test_bmp(16,6)+test_bmp(16,8))/2 ;
    test_bmp(25,2)=(test_bmp(25,1)+test_bmp(25,3))/2 ;
    bmp{i} = test_bmp ;
end
