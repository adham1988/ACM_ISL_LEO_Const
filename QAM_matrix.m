clear all
close all
clc

M = 2^4

p = '000001';

MQAM_encode(p,M);


function x = MQAM_encode(p,M)
L = sqrt(M);
k = log2(M);
q = quantizer([k 0],'Mode','ufixed');

%Make QAM matrix
for j = 1:L
    for i = 1:L
        a(j,i) = -L + 1 + 2*(i-1);
        b(i,j) = -L + 1 + 2*(i-1);
    end
end

figure
hold on
box on

%Make symbol vectors
for i = 1:M
    j = bin2gray(i-1,'qam',M)+1;
    c(:,j) = [a(i) b(i)];
    plot(c(1,j),c(2,j),'b.')
    text(c(1,j)+0.1,c(2,j),sprintf("%i",j))
end

%Outputs the symbol vector 'x' for a given 'p'
p = bin2num(q,p);
for i = 1:M
    if p == i
        x = c(:,i)
        plot(c(1,i),c(2,i),'r*')
    end
end


axis equal
axis([-L L+1 -L L])
xline(0)
yline(0)
end

function c_hat = MQAM_decode(y,M)








end




















