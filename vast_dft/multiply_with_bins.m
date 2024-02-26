function [output] = multiply_with_bins(x,h,q,nfft)

%output är 4096/2 lång. 
n=nfft/2;

numbins=121;
bin_length=floor(n/numbins); 
output=zeros(1,n);
x=x';
x=x(1:n);
h=h(:,1:n);


for bin=1:numbins-1
    indexes=(bin-1)*bin_length+1:bin*bin_length;
  
    for speaker = 1:16
        output(indexes)= output(indexes)+ x(indexes).*h(speaker,indexes).*q(bin);
    end  
end

end