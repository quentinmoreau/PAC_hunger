function [M_RN, l2_l1, M_H, M_S]=tfr_measures(faslt)

M_RN=sum(sum(((abs(faslt)).^4)))/(sum(sum(abs(faslt).^2)))^2;

l2_l1= sqrt(sum(sum(abs(faslt).^2)))/sum(sum(abs(faslt)));

H=size(faslt, 1)*size(faslt, 2);

M_H=(sqrt(H)-(sum(sum(abs(faslt))))/(sqrt(sum(sum(abs(faslt).^2)))))*(sqrt(H)-1)^-1;

%1\H

M_S=(1/H)*(sum(sum( sqrt(abs((faslt)/(sum(sum((faslt)))) ) ) )))^2;

%{
figure;
plot(1, M_RN, 'o');
hold on
plot(2, l2_l1, 'o');
hold on
plot(3, M_H, 'o');
ylim([-0.1, 2*max([M_RN l2_l1 M_H])])
xlim([0 4])
%hold on
%plot(4, M_S, 'o');
%}

%R=renyi(abs(faslt))

end