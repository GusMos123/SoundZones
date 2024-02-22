function [output] = multiply_w_bins_safe(x,h,q,nfft)

%output är 4096 lång. 

numbins=121;
bin_length=ceil(nfft/numbins); 
output=zeros(1,nfft);
x=x';


for bin=1:numbins-1
    indexes=(bin-1)*bin_length+1:bin*bin_length;
    for speaker = 1:16
        output(indexes)= output(indexes)+ x(indexes).*h(speaker,indexes).*q(bin);
    end  
end

end