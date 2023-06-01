
% Antal bolag med fullständig data för FF
for i=1:height(mcap)

    mask_all(i,:)=sum((isnan(mcap(i,:)) | isnan(returns(i,:)) | isnan(bm(ceil(i/12),:))...
        | isnan(op(ceil(i/12),:)) |isnan(inve(ceil(i/12),:)))==0);        
end

% Plot observations used
figure;
plot(month, mask_all, LineWidth=2, Color="b")
hold on
%plot(month, mask_ret, 'r-', LineWidth=2)
ylabel('Number of stock observations')
ylim([1000 2000])
xlabel('Year')
legend('Post-screening for constructing the FF-factors.', 'Location', 'northwest')
set(gca, 'FontSize', 16);

%Includes following varibles Inv, OP, BM, MC, Returns
% max(mask_all) ger högsta värdet på observationer
%% Antal bolag med fullständig data för returns 
for i=1:height(returns)
    mask_ret(i,:)= sum((isnan(returns(i,:)))==0);        
end
max(mask_all)
min(mask_all)
%%
%% Antal bolag med esg data för  
for i=1:height(esg)
    mask_esg(i,:)= sum((isnan(esg(i,:)))==0);        
end
max(mask_esg)
min(mask_esg)

% Plot ESG rating
figure;  hold all;
plot(month, mask_esg, 'color', [31/255 183/255 56/255],'LineWidth',2)
legend('Evolution of ESG rated firms', 'Location','northwest')
xlabel('Year')
ylabel('ESG-firms listed on the NYSE and AMEX')
set(gca, 'FontSize', 16);

%%
% end
for i=1:height(mcap)

    mask_bm(i,:)=sum((isnan(mcap(i,:)) | isnan(returns(i,:)) | isnan(bm(i,:)))==0);        
end

for i=1:height(mcap)

    mask_op(i,:)=sum((isnan(mcap(i,:)) | isnan(returns(i,:)) | isnan(op(i,:)))==0);        
end

for i=1:height(mcap)

    mask_inv(i,:)=sum((isnan(mcap(i,:)) | isnan(returns(i,:)) | isnan(inve(i,:)))==0);        
end

figure; hold all;
plot(mask_all);
plot(mask_bm);
plot(mask_op);
plot(mask_inv);
legend("all", "bm", "op", "inv");