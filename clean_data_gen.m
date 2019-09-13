function clean_data = clean_data_gen(img_real_line,rand_num, num_of_delta_alpha,k_max)

clean_data = [];

for i = 1:1:2*k_max+1      
    index = (i-num_of_delta_alpha:1:i+num_of_delta_alpha)-1;
    tmp = (img_real_line(:,mod(index,2*k_max+1)+1));
    clean_data = [clean_data,repmat(tmp,1,rand_num(i))];    
end


