function [estimated_coeff,estimated_pmf,tmp] = coeff_pmf_recon(x_init,p_init,index_coeff_full,k_max,coeff,zeropad,pmf_theta)

%% esitimated coeff (with global rotation)
x_cell = mat2cell(x_init,index_coeff_full,1);
x_cell = cellfun(@(x,y) (x+conj(y))./2,x_cell,flipud(x_cell),'UniformOutput',false);
hat_a = x_cell(k_max+1:2*k_max+1,1);
tmp = hat_a;

if real(hat_a{1}(1))<0
    hat_a = cellfun(@(x) -1.*x,hat_a,'UniformOutput',false);
end

%% estimated FFT pmf (with global rotation)
hat_FFT_pmf = [1;(p_init+flipud(conj(p_init)))./2];


%% find the global rotation
 for i = -k_max:1:k_max    
     if i<0     
         temp(i+k_max+1,1) = conj(hat_a{abs(i)+1}).'*coeff{abs(i)+1};    
     else        
         temp(i+k_max+1,1) = hat_a{i+1}.'*conj(coeff{i+1});    
     end    
 end
temp = [temp;zeros(zeropad-size(temp,1),1)]; %pad 0's to condcut FFT for analysis
[row,~] = find(abs(fft(temp)) == max(abs(fft(temp))));
globalrotation = (row(1,1)-1)./zeropad*2*pi;

%% rotate the coeff and pmf
for i =1:1:k_max+1      
    estimated_coeff{i,1}=hat_a{i}.*exp(-sqrt(-1)*(i-1)*globalrotation);       
end

estimated_pmf = ifft(hat_FFT_pmf);
estimated_pmf = circshift(estimated_pmf,round((row-1)./zeropad*size(pmf_theta,1)));


end