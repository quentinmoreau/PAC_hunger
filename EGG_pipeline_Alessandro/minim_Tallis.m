function [en,a]=minim_Tallis(sig, n, plots)

en=zeros(n, 1);
for i=1:n
    en(i)=timeE_Tsallis_entropy(sig, i/n);
end

[a, b]=min(en);

if strcmp(plots, 'on')

figure;
plot((1:n)/n, en, 'o')
%plot(en)
title(['The minimus Tsallis entropy is: ', num2str(a), ' occuring at q= ', num2str(b/n)])
xlabel('q')
ylabel('Normalized Tsallis entropy')
xline(b/n, '--r', 'LineWidth', 2)
end

end