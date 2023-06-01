% PLOT Carhart Low/high, SP500, FF LOW/HIGH and Market


%% Plot port/market/sp500 with adjusted Y axis to % percentage-points
% SP500 returns
sp500data=readmatrix("sp500returns.xlsx", "Sheet", "Blad1");

cum_port=cumsum((portfolio_ex)); %portfolio_ex = high portfolio
cum_mark=cumsum((mkt_e));
cum_low_portfolio=cumsum(low_portfolio_ex);

sp500_excess= sp500data-rf
sp500_ex=log(sp500_excess+1)
cum_sp500=cumsum((sp500_ex));

avg_ret_sp500=mean(sp500_ex);
variance_sp500=var(sp500_ex);
Volatility_sp500=std(sp500_ex);
sharpe_sp500=avg_ret_low/Volatility_low;


figure;
hold all;
plot(month, 100*(cum_port, linewidth2));
plot(month, 100*(cum_low_portfolio));
plot(month, 100*(cum_mark));
plot(month, 100*(cum_sp500), 'color', 'black');
legend('High portfolio, 10% best ESG score','Low portfolio, 10% lowest ESG score', 'Market', 'S&P 500, Risk adjusted');
xlabel('Year');
ylabel('Sum of logarithmic returns');
ytickformat('%g%%'); % Set y-axis tick labels to percentage points
