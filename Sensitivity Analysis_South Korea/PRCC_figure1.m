load PRCC_HepaB_total.mat

%figure('Unit','centimeters','Position',[10,10,15,10])
%tiledlayout(1,1,'padding','compact','TileSpacing','Compact')
nexttile(2)
len_s = length(s);

for r=1:len_s
    c1=['PRCCs at time = ' num2str(s(r))];
    a=find(uncorrected_sign(r,:)<alpha);
    ['Significant PRCCs'];
    a;
    PRCC_var(a);
    prcc(r,a);
    b=num2str(prcc(r,a));
    sign_label_struct.index{r}=a;
    sign_label_struct.label{r}=PRCC_var(a);
    sign_label_struct.value{r}=b;
    %% Plots of PRCCs and their p-values
    hold on
    %Title(c1),set(gca,'XTickLabel',PRCC_var,'XTick',[1:k]),Title('PRCCs');
end

b=bar(PRCCs(:,:)');
for i=1:len_s
    b(i).FaceColor = (len_s-i)/len_s*[.3 .5 .9]+i/len_s*[.9 .5 .3];
end
title('PRCC','FontSize',12);
set(gca,'FontName','Times','XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14],...
    'XTickLabel',{'b','\mu','d','\beta','\epsilon','k','r_1','r_2','\sigma','w','p','q','v','dummy'} ...
    ,'FontSize',12,'YGrid','on','YMinorGrid','on')
