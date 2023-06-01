Carhart=readmatrix("FF2.xlsx", "sheet", "Carhart");
Carhart=Carhart/100;
rf=Carhart(:,5);
Carhart=log(Carhart+1);
%% SP500 returns
sp500returns=readmatrix("sp500returns.xlsx", "Sheet", "Blad1");

%%
month_matrix = reshape(month, 12, 19);
years=year(month_matrix(1,:))';

ESG=readmatrix("ESG.xlsx", "Sheet", "ESG");
ESG=ESG(2:20,2:end);




%% Sorting 

for i=1:width(ESG)
    if all(ESG(1:18,i) == ESG(19,i)) 
        ESG(1:18,i) = NaN;
    end

end

for i=1:height(ESG)
    
    start_month=(i-1)*12+1;
    end_month=start_month+11;

    missing = any(isnan(returns(start_month:end_month,:)))... 
        | any(isnan(mcap(start_month:end_month,:)));

    ESG(i, missing)= NaN;
end

for i=1:height(ESG)

    start_month=(i-1)*12+1;
    end_month=start_month+11;

    esg(start_month:end_month,:)=repmat(ESG(i,:),end_month-start_month+1,1);
    market_cap(start_month:end_month,:)=repmat(cap(i,:),end_month-start_month+1,1);
end


for i=1:height(ESG)
   
    start_month=(i-1)*12+1;
    end_month=start_month+11;

    nan_mask=any(~isnan(esg(start_month:end_month,:)),1)...
        | any(~isnan(returns(start_month:end_month,:)),1)...
        | any(~isnan(mcap(start_month:end_month,:)),1)...
        | any(~isnan(market_cap(start_month:end_month,:)),1);

    mask(start_month:end_month,:)=repmat(nan_mask,end_month-start_month+1,1);
    mask_control(start_month:end_month,:)=sum(mask(start_month,:));
end

esg(mask==0) = NaN;
returns(mask == 0) = NaN;
mcap(mask == 0) = NaN;
market_cap(mask == 0) = NaN;


%% Portfolio Construction High ESG

[T,N]=size(returns);
portfolio_returns=zeros(T,N);
portfolio=zeros(T,1);
w=zeros(size(returns));
portfolio_control=zeros(T,1);
w_control=zeros(T,1);

for i=1:height(returns)
    
    [ESG_sorted, index] = sort(esg(i,:),'descend', MissingPlacement='last');
    
    ESG_quantile = prctile(ESG_sorted, 90 , 2);

    mask = esg(i,:) >= ESG_quantile & ~isnan(returns(i,:)) & ~isnan(mcap(i,:));
    
    portfolio_control(i,:)=sum(mask);


    portfolio_returns(i,:)=returns(i,:).*mask;
   
    portfolio_returns(i,mask == 0) = NaN;

    
 % w(i,mask)=market_cap(i,mask)/sum(market_cap(i,mask)); 
  w(i,mask)=mcap(i,mask)/sum(mcap(i,mask)); 
   % w(i,mask)=1/sum(mask);
    portfolio_returns(isnan(portfolio_returns)) = 0;
    w(isnan(w)) = 0;
    w_control(i,:)=sum(w(i,:));

    portfolio_control(i,1)=sum(mask);
    portfolio_control(i,2)=sum(portfolio_returns(i,:) ~=0);
    portfolio_control(i,3)=sum(w(i,:) ~= 0);
    portfolio_control(i,4)=sum(portfolio_returns(i,mask) == 0);


    portfolio(i,1)=sum(w(i,:).*portfolio_returns(i,:));
end


portfolio_exess = portfolio - rf;
portfolio_ex = log(portfolio_exess +1);
%% Sharpe Ratio

avg_ret=mean(portfolio_ex);
variance=var(portfolio_ex);
Volatility=std(portfolio_ex);
sharpe=avg_ret/Volatility;
                                    % "Generally, the higher the Sharpe ratio,
                                    % the more attractive the risk-adjusted
                                    % return."
                                    % Från investopedia.
                                    % dvs våra portföljer har sämre
                                    % riskjusterad avkastning

avg_ret_m=mean(mkt_e);
variance_m=var(mkt_e);
Volatility_m=std(mkt_e);
sharpe_m=avg_ret_m/Volatility_m;

%% Full Time-horizon Regression ( Nummer 1) * High

cum_port=cumsum((portfolio_ex));
cum_mark=cumsum((mkt_e));
sp500_log=log(sp500returns+1);
cum_sp500=cumsum((sp500_log));

% Plot port/market/sp500
%figure;  hold all;
%plot(month, cum_port);
%plot(month, cum_mark);
%plot(month, cum_sp500, 'color', 'black');
%legend('Portfolio', 'Market', 'S&P 500');
%xlabel('Year');
%ylabel('Return')
%ytickformat('%g%%'); % Set y-axis tick labels to percentage points

% Plot port/market/sp500 with adjusted Y axis to % percentage-points
%figure;
%hold all;
%plot(month, 100*(cum_port));
%plot(month, 100*(cum_mark));
%plot(month, 100*(cum_sp500), 'color', 'black');
%legend('Portfolio', 'Market', 'S&P 500, Inflation Adjusted');
%xlabel('Year');
%ylabel('Cum. Logarithimic returns, Carhart Four Factor High ESG Portfolio');
%ytickformat('%g%%'); % Set y-axis tick labels to percentage points



% Regression using OLS
Reg_table=table(mkt_e, smb, hml, wml, portfolio_ex, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'WML', 'Portfolio'} );
Regression_Monthly_C = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + WML'); %,...

disp(Regression_Monthly_C)

% Test for Homoscedasticity - Engle's Arch Test  
[h, Pvalue, stat, cValue] = archtest(Regression_Monthly_C.Residuals.Raw, Lags=[1 3]) 
% h=0 -> H0 not rejected = homo, men p= 0.8862. Hac for safety
% med 3 laggar, h=1 -> H0 rejected = hetero, p= 0.0257. Hac

% Test for Autocorrelation
[h, Pvalue, stat, cValue] = lbqtest(Regression_Monthly_C.Residuals.Raw)
% h=0 -> H0 not rejected = ingen autokorrelation, men p=0.575 Hac for safety

% Regression using Hac
[EstCoeffCov, se, coeff] = hac(Regression_Monthly_C)
t_hac = coeff./se
p_hac = erfc(0.7071*abs(t_hac))

HAC_table_C = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'WML'});

disp(HAC_table_C)

% Regression using OLS Talwar weighting
Regression_Monthly_Talwar_C = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + WML', ...
    'RobustOpts', 'talwar'); 
disp(Regression_Monthly_Talwar_C)

% Annualizing
Annualized_alpha_ols_hac = Regression_Monthly_C.Coefficients.Estimate(1) *12
Annualized_alpha_talwar = Regression_Monthly_Talwar_C.Coefficients.Estimate(1) *12
%Alpha_M = alpha * 10000          %Pricing error in Basis Points

% Transforming into basis points
Alpha_M_ols_hac_C = Annualized_alpha_ols_hac * 10000          %Pricing error in Basis Points
Alpha_M_Talwar_C = Annualized_alpha_talwar * 10000          %Pricing error in Basis Points

outliers=find(isoutlier(Regression_Monthly_C.Residuals.Raw))
outliers=find(isoutlier(Regression_Monthly_C_Talwar.Residuals.Raw))

anova(Regression_FF_Monthly)
anova(Regression_FF_Monthly_Talwar)

%% Monthly Regression 2003 - 2011 

pt_cum_port=cumprod((portfolio_ex(1:108,:))+1);
pt_cum_mark=cumprod(mkt_e(1:108,:)+1);

% Graph on performance 2010 - 2021
figure;  hold all;
plot(pt_cum_port);
plot(pt_cum_mark);

% Part time Sharpe

avg_ret_pt=mean(portfolio_ex(1:108,:));
variance_pt=var(portfolio_ex(1:108,:));
Volatility_pt=std(portfolio_ex(1:108,:));
sharpe_pt=avg_ret_pt/Volatility_pt;


avg_ret_pt_mkt=mean(mkt_e(1:108,:));
variance_pt_mkt=var(mkt_e(1:108,:));
Volatility_pt_mkt=std(mkt_e(1:108,:));
sharpe_pt_mkt=avg_ret_pt_mkt/Volatility_pt_mkt;

% Regression
pt_mkt_e = mkt_e(1:108,:);
pt_smb = smb(1:108,:);
pt_hml = hml(1:108,:);
pt_wml = wml(1:108,:);
pt_portfolio_ex = portfolio_ex(1:108,:);


% Regression using OLS
Reg_table=table(pt_mkt_e, pt_smb, pt_hml, pt_wml, pt_portfolio_ex, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'WML', 'Portfolio'} );
Regression_Part_Time_Monthly_2003_C = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + WML');

disp(Regression_Part_Time_Monthly_2003_C)

% Regression using OLS w HAC standard errors
[EstCoeffCov, se, coeff] = hac(Regression_Part_Time_Monthly_2003_C);
t_hac = coeff./se;
p_hac = erfc(0.7071*abs(t_hac));

HAC_table_Part_Time2003_C = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'WML'});
disp(HAC_table_Part_Time2003_C)

% Regression using OLS talwAR
Regression_Part_Time_Monthly_Talwar_2003_C = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + WML', ...
    'RobustOpts', 'talwar'); 
disp(Regression_Part_Time_Monthly_Talwar_2003_C)

% Annualizing
Annualized_alpha_ols_hac_pt = Regression_Part_Time_Monthly_2003_C.Coefficients.Estimate(1) *12
Annualized_alpha_talwar_pt = Regression_Part_Time_Monthly_Talwar_2003_C.Coefficients.Estimate(1) *12

% Transforming into annual basis points       
Alpha_M_ols_hac_pt2003_C = Annualized_alpha_ols_hac_pt * 10000          %Pricing error in Basis Points
Alpha_M_Talwar_pt2003_C = Annualized_alpha_talwar_pt * 10000 



%% Monthly Regression 2012 - 2021

pt_cum_port=cumprod((portfolio_ex(109:end,:))+1);
pt_cum_mark=cumprod(mkt_e(109:end,:)+1);

% Graph on performance 2012 - 2021
figure;  hold all;
plot(pt_cum_port);
plot(pt_cum_mark);

% Part time Sharpe

avg_ret_pt=mean(portfolio_ex(109:end,:));
variance_pt=var(portfolio_ex(109:end,:));
Volatility_pt=std(portfolio_ex(109:end,:));
sharpe_pt=avg_ret_pt/Volatility_pt;


avg_ret_pt_mkt=mean(mkt_e(109:end,:));
variance_pt_mkt=var(mkt_e(109:end,:));
Volatility_pt_mkt=std(mkt_e(109:end,:));
sharpe_pt_mkt=avg_ret_pt_mkt/Volatility_pt_mkt;

% Regression
pt_mkt_e = mkt_e(109:end,:);
pt_smb = smb(109:end,:);
pt_hml = hml(109:end,:);
pt_wml = wml(109:end,:);
pt_portfolio_ex = portfolio_ex(109:end,:);

% Regression using OLS
Reg_table=table(pt_mkt_e, pt_smb, pt_hml, pt_wml, pt_portfolio_ex, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'WML', 'Portfolio'} );
Regression_Part_Time_Monthly_2021_C = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + WML');

disp(Regression_Part_Time_Monthly_2021_C)

% Regression using OLS w HAC standard errors
[EstCoeffCov, se, coeff] = hac(Regression_Part_Time_Monthly_2021_C);
t_hac = coeff./se;
p_hac = erfc(0.7071*abs(t_hac));

HAC_table_Part_Time_2021_C = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'WML'});
disp(HAC_table_Part_Time_2021_C)

% Regression using OLS talwAR
Regression_Part_Time_Monthly_Talwar_2021_C = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + WML', ...
    'RobustOpts', 'talwar'); 
disp(Regression_Part_Time_Monthly_Talwar_2021_C)

% Annualizing
Annualized_alpha_ols_hac_pt = Regression_Part_Time_Monthly_2021_C.Coefficients.Estimate(1) *12;
Annualized_alpha_talwar_pt = Regression_Part_Time_Monthly_Talwar_2021_C.Coefficients.Estimate(1) *12;

% Transforming into annual basis points       
Alpha_M_ols_hac_pt2021_C = Annualized_alpha_ols_hac_pt * 10000          %Pricing error in Basis Points
Alpha_M_Talwar_pt2021_C = Annualized_alpha_talwar_pt * 10000 


%% Regression Using Data from Fama-French-Kenneth

cum_port=cumsum((portfolio_ex));
cum_C=cumsum((Carhart(:,1)));

figure;  hold all;
plot(cum_port);
plot(cum_C);

C_mkt = Carhart(:,1);
C_smb = Carhart(:,2);
C_hml = Carhart(:,3);
C_wml = Carhart(:,4);


Reg_table=table(C_mkt, C_smb, C_hml, C_wml, portfolio_ex, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'WML', 'Portfolio'} );
Regression_FF_Monthly_C = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + WML');

disp(Regression_FF_Monthly_C)

% Test for Homoscedasticity - Engle's Arch Test  
[h, Pvalue, stat, cValue] = archtest(Regression_FF_Monthly_C.Residuals.Raw, Lags=[1 3]) 
% h=0 -> H0 not rejected = homo, men p= 0.616. Hac for safety
% med 3 laggar, h=0 -> H0 not rejected = homo, p= 0.3601. Hac

% Test for Autocorrelation
[h, Pvalue, stat, cValue] = lbqtest(Regression_FF_Monthly_C.Residuals.Raw)
% h=0 -> H0 not rejected = ingen autokorrelation, men p=0.3745 Hac for safety

[EstCoeffCov, se, coeff] = hac(Regression_FF_Monthly_C)
t_hac = coeff./se
p_hac = erfc(0.7071*abs(t_hac))

HAC_table_FF_C = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'WML'});

disp(HAC_table_FF_C)

Regression_FF_Monthly_Talwar_C = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + WML', ...
    'RobustOpts', 'talwar'); 
disp(Regression_FF_Monthly_Talwar_C)


Annualized_alpha_ols_hac = Regression_FF_Monthly_C.Coefficients.Estimate(1) *12
Annualized_alpha_talwar = Regression_FF_Monthly_Talwar_C.Coefficients.Estimate(1) *12
%Alpha_M = alpha * 10000          %Pricing error in Basis Points

% Transforming into annual basis points       
Alpha_M_ols_hac_FF_C = Annualized_alpha_ols_hac_pt * 10000      %Pricing error in Basis Points
Alpha_M_Talwar_FF_C = Annualized_alpha_talwar_pt * 10000 
       


%% Low portfolio 10% (Nummer 2*)

[T,N]=size(returns);
low_portfolio_returns=zeros(T,N);
low_portfolio=zeros(T,1);
w=zeros(size(returns));
low_portfolio_control=zeros(T,1);
w_control=zeros(T,1);
for i=1:height(returns)
    
    [ESG_sorted, index] = sort(esg(i,:),'ascend', MissingPlacement='last');
    
    
    ESG_quantile = prctile(ESG_sorted, 10, 2); %10 for low-portfolio

    mask = (esg(i,:) <= ESG_quantile) & ~isnan(returns(i,:)) & ~isnan(mcap(i,:));
    
    

    low_portfolio_control(i,:)=sum(mask);

    low_portfolio_returns(i,:)=returns(i,:).*mask;
   
    low_portfolio_returns(i,mask == 0) = NaN;

    
  %w(i,mask)=market_cap(i,mask)/sum(market_cap(i,mask)); 
  w(i,mask)=mcap(i,mask)/sum(mcap(i,mask)); 
   % w(i,mask)=1/sum(mask);
    low_portfolio_returns(isnan(low_portfolio_returns)) = 0;
    w(isnan(w)) = 0;
    w_control(i,:)=sum(w(i,:));

    portfolio_control(i,1)=sum(mask);
    portfolio_control(i,2)=sum(low_portfolio_returns(i,:) ~=0);
    portfolio_control(i,3)=sum(w(i,:) ~= 0);
    portfolio_control(i,4)=sum(low_portfolio_returns(i,:)==0 & mask == 1);

    
    low_portfolio(i,1)=sum(w(i,:).*low_portfolio_returns(i,:));
end

low_portfolio_excess = low_portfolio-rf;
low_portfolio_ex = log(low_portfolio_excess+1);

%% Low portfolio Graphs


low_cum_port=cumsum((low_portfolio_ex));
cum_mark=cumsum((mkt_e));

figure;  hold all;
plot(low_cum_port);
plot(cum_mark);

% Plot port/market/sp500
figure;  hold all;
plot(month, cum_port);
plot(month, cum_mark);
plot(month, cum_sp500, 'color', 'black');
legend('Portfolio', 'Market', 'S&P 500');
xlabel('Year');
ylabel('Return')

% Plot port/market/sp500 with adjusted Y axis to % percentage-points
figure;
hold all;
plot(month, 100*(low_cum_port));
plot(month, 100*(cum_mark));
plot(month, 100*(cum_sp500), 'color', 'black');
legend('Portfolio', 'Market', 'S&P 500, Inflation Adjusted');
xlabel('Year');
ylabel('Returns, Carhart Four Factor Low ESG Portfolio');
ytickformat('%g%%'); % Set y-axis tick labels to percentage points



% Low portfolio Sharpe
avg_ret_low=mean(low_portfolio_ex);
variance_low=var(low_portfolio_ex);
Volatility_low=std(low_portfolio_ex);
sharpe_low=avg_ret_low/Volatility_low;

%% Low portfolio Full Time-horizon Regression

Reg_table=table(mkt_e, smb, hml, wml, low_portfolio_ex, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'WML', 'Portfolio'} );
Regression_Monthly_Low_C = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + WML'); %,...

disp(Regression_Monthly_Low_C)

% Test for Homoscedasticity - Engle's Arch Test  
[h, Pvalue, stat, cValue] = archtest(Regression_Monthly_Low_C.Residuals.Raw, Lags=[1 3]) 
% h=0 -> H0 not rejected = homo, men p= 0.5441. Hac for safety
% med 3 laggar, h=0 -> H0 not rejected = homo, p= 0.813. Hac

% Test for Autocorrelation
[h, Pvalue, stat, cValue] = lbqtest(Regression_Monthly_Low_C.Residuals.Raw)
% h=0 -> H0 not rejected = ingen autokorrelation, men p=0.5023 Hac for safety

% Regression using Hac
[EstCoeffCov, se, coeff] = hac(Regression_Monthly_Low_C);
t_hac = coeff./se;
p_hac = erfc(0.7071*abs(t_hac));

HAC_table_Low_C = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'WML'});

disp(HAC_table_Low_C)

% Regression using OLS Talwar weighting
Regression_Monthly_Talwar_Low_C = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + WML', ...
    'RobustOpts', 'talwar'); 
disp(Regression_Monthly_Talwar_Low_C)

% Annualizing
Annualized_alpha_ols_hac = Regression_Monthly_Low_C.Coefficients.Estimate(1) *12
Annualized_alpha_talwar = Regression_Monthly_Talwar_Low_C.Coefficients.Estimate(1) *12
%Alpha_M = alpha * 10000          %Pricing error in Basis Points

% Transforming into annual basis points       
Alpha_M_ols_hac_Low_C = Annualized_alpha_ols_hac * 10000      %Pricing error in Basis Points
Alpha_M_Talwar_Low_C = Annualized_alpha_talwar * 10000 

%% Monthly Part Time regression Low Portfolio 2003 - 2011
pt_cum_port=cumprod((low_portfolio_ex(1:108,:))+1);
pt_cum_mark=cumprod(mkt_e(1:108,:)+1);

% Graph on performance 2003 - 2011
figure;  hold all;
plot(pt_cum_port);
plot(pt_cum_mark);

% Part time Sharpe

avg_ret_pt_low=mean(low_portfolio_ex(1:108,:));
variance_pt_low=var(low_portfolio_ex(1:108,:));
Volatility_pt_low=std(low_portfolio_ex(1:108,:));
sharpe_pt_low=avg_ret_pt_low/Volatility_pt_low;

% Regression
pt_mkt_e = mkt_e(1:108,:);
pt_smb = smb(1:108,:);
pt_hml = hml(1:108,:);
pt_wml = wml(1:108,:);
pt_portfolio_ex = low_portfolio_ex(1:108,:);



% Regression using OLS
Reg_table=table(pt_mkt_e, pt_smb, pt_hml, pt_wml, pt_portfolio_ex, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'WML', 'Portfolio'} );
Regression_Part_Time_Monthly_Low_2003_C = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + WML');

disp(Regression_Part_Time_Monthly_Low_2003_C)

% Regression using OLS w HAC standard errors
[EstCoeffCov, se, coeff] = hac(Regression_Part_Time_Monthly_Low_2003_C);
t_hac = coeff./se;
p_hac = erfc(0.7071*abs(t_hac));

HAC_table_Part_Time_Low2003_C = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'WML'});
disp(HAC_table_Part_Time_Low2003_C)

% Regression using OLS talwAR
Regression_Part_Time_Monthly_Low_Talwar_2003_C = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + WML', ...
    'RobustOpts', 'talwar'); 
disp(Regression_Part_Time_Monthly_Low_Talwar_2003_C)

% Annualizing
Annualized_alpha_ols_hac_pt = Regression_Part_Time_Monthly_Low_2003_C.Coefficients.Estimate(1) *12
Annualized_alpha_talwar_pt = Regression_Part_Time_Monthly_Low_Talwar_2003_C.Coefficients.Estimate(1) *12

% Transforming into annual basis points       
Alpha_M_ols_hac_pt_Low2003_C = Annualized_alpha_ols_hac_pt * 10000          %Pricing error in Basis Points
Alpha_M_Talwar_pt_Low2003_C = Annualized_alpha_talwar_pt * 10000 




%% Monthly Part Time regression Low Portfolio 2012 - 2021
pt_cum_port=cumsum((low_portfolio_ex(109:end,:)));
pt_cum_mark=cumsum(mkt_e(109:end,:));

% Graph on performance 2012 - 2021
figure;  hold all;
plot(pt_cum_port);
plot(pt_cum_mark);

% Part time Sharpe

avg_ret_pt_low=mean(low_portfolio_ex(109:end,:));
variance_pt_low=var(low_portfolio_ex(109:end,:));
Volatility_pt_low=std(low_portfolio_ex(109:end,:));
sharpe_pt_low=avg_ret_pt_low/Volatility_pt_low;

% Regression
pt_mkt_e = mkt_e(109:end,:);
pt_smb = smb(109:end,:);
pt_hml = hml(109:end,:);
pt_wml = wml(109:end,:);
pt_portfolio_ex = low_portfolio_ex(109:end,:);

% Regression using OLS
Reg_table=table(pt_mkt_e, pt_smb, pt_hml, pt_wml, pt_portfolio_ex, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'WML', 'Portfolio'} );
Regression_Part_Time_Monthly_Low_2021_C = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + WML');

disp(Regression_Part_Time_Monthly_Low_2021_C)

% Regression using OLS w HAC standard errors
[EstCoeffCov, se, coeff] = hac(Regression_Part_Time_Monthly_Low_2021_C);
t_hac = coeff./se;
p_hac = erfc(0.7071*abs(t_hac));

HAC_table_Part_Time_Low2021_C = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'WML'});
disp(HAC_table_Part_Time_Low2021_C)

% Regression using OLS talwAR
Regression_Part_Time_Monthly_Low_Talwar_2021_C = fitlm(Reg_table,...
    'Portfolio ~ Market + SMB + HML + WML', 'RobustOpts', 'talwar'); 
disp(Regression_Part_Time_Monthly_Low_Talwar_2021_C)

% Annualizing
Annualized_alpha_ols_hac_pt = Regression_Part_Time_Monthly_Low_2021_C.Coefficients.Estimate(1) *12
Annualized_alpha_talwar_pt = Regression_Part_Time_Monthly_Low_Talwar_2021_C.Coefficients.Estimate(1) *12

% Transforming into annual basis points       
Alpha_M_ols_hac_pt_Low2021_C = Annualized_alpha_ols_hac_pt * 10000          %Pricing error in Basis Points
Alpha_M_Talwar_pt_Low2021_C = Annualized_alpha_talwar_pt * 10000 




%% Regression Using Data from Fama-French-kenneth Low Portfolio

cum_port=cumsum((low_portfolio_ex));
cum_C=cumsum((Carhart(:,1)));

figure;  hold all;
plot(cum_port);
plot(cum_C);


Reg_table=table(C_mkt, C_smb, C_hml, C_wml, low_portfolio_ex, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'WML', 'Portfolio'} );
Regression_Monthly_Low_FF_C = fitlm(Reg_table,...
    'Portfolio ~ Market + SMB + HML + WML'); %,...

disp(Regression_Monthly_Low_FF_C)

% Test for Homoscedasticity - Engle's Arch Test  
[h, Pvalue, stat, cValue] = archtest(Regression_Monthly_Low_FF_C.Residuals.Raw, Lags=[1 3]) 
% h=0 -> H0 not rejected = homo, men p= 0.93. Hac for safety
% med 3 laggar, h=0 -> H0 not rejected = homo, p= 0.7333. Hac

% Test for Autocorrelation
[h, Pvalue, stat, cValue] = lbqtest(Regression_Monthly_Low_FF_C.Residuals.Raw)
% h=0 -> H0 not rejected = ingen autokorrelation, men p=0.5152 Hac for safety

% Regression using Hac
[EstCoeffCov, se, coeff] = hac(Regression_Monthly_Low_FF_C)
t_hac = coeff./se
p_hac = erfc(0.7071*abs(t_hac))

HAC_table_Low_FF_C = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'WML'});

disp(HAC_table_Low_FF_C)

% Regression using OLS Talwar weighting
Regression_Monthly_Talwar_Low_FF_C = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + WML', ...
    'RobustOpts', 'talwar'); 
disp(Regression_Monthly_Talwar_Low_FF_C)

% Annualizing
Annualized_alpha_ols_hac = Regression_Monthly_Low_FF_C.Coefficients.Estimate(1) *12
Annualized_alpha_talwar = Regression_Monthly_Talwar_Low_FF_C.Coefficients.Estimate(1) *12
%Alpha_M = alpha * 10000          %Pricing error in Basis Points

% Transforming into annual basis points       
Alpha_M_ols_hac_Low_FF_C = Annualized_alpha_ols_hac * 10000      %Pricing error in Basis Points
Alpha_M_Talwar_Low_FF_C = Annualized_alpha_talwar * 10000 



%% Long - Short Portfolio

ls_portfolio = portfolio - low_portfolio; 

cum_ls_port=cumsum((ls_portfolio));
cum_mark=cumsum((mkt_e));

figure;  hold all;
plot(month, cum_ls_port, 'Linewidth', 2);
plot(month, cum_mark, 'Linewidth', 2);
%% Long - Short Portfolio

sl_portfolio = low_portfolio - portfolio; 

cum_sl_port=cumsum((sl_portfolio));
cum_mark=cumsum((mkt_e));

figure;  hold all;
plot(month, cum_ls_port, 'Linewidth', 2);
plot(month, cum_mark, 'Linewidth', 2);
%%
% Long - Short Sharpe

avg_ret_ls=mean(ls_portfolio);
variance_ls=var(ls_portfolio);
Volatility_ls=std(ls_portfolio);
sharpe_ls=avg_ret_ls/Volatility_ls;

%% Long Short Monthly Regression (hög-låg portfolio) " Nummer 3"

Reg_table=table(mkt_e, smb, hml, wml, ls_portfolio, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'WML', 'Portfolio'} );
Regression_Monthly_LS_C = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + WML'); %,...

disp(Regression_Monthly_LS_C)

% Test for Homoscedasticity - Engle's Arch Test  
[h, Pvalue, stat, cValue] = archtest(Regression_Monthly_LS_C.Residuals.Raw, Lags=[1 3]) 
% h=0 -> H0 not rejected = homo, men p= 0.5653. Hac for safety
% med 3 laggar, h=0 -> H0 not rejected = homo, p= 0.7017. Hac

% Test for Autocorrelation
[h, Pvalue, stat, cValue] = lbqtest(Regression_Monthly_LS_C.Residuals.Raw)
% h=0 -> H0 not rejected = ingen autokorrelation, men p=0.2518 Hac for safety

% Regression using Hac
[EstCoeffCov, se, coeff] = hac(Regression_Monthly_LS_C);
t_hac = coeff./se;
p_hac = erfc(0.7071*abs(t_hac));

HAC_table_LS_C = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'WML'});

disp(HAC_table_LS_C)

% Regression using OLS Talwar weighting
Regression_Monthly_Talwar_LS_C = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + WML', ...
    'RobustOpts', 'talwar'); 
disp(Regression_Monthly_Talwar_LS_C)

% Annualizing
Annualized_alpha_ols_hac = Regression_Monthly_LS_C.Coefficients.Estimate(1) *12
Annualized_alpha_talwar = Regression_Monthly_Talwar_LS_C.Coefficients.Estimate(1) *12
%Alpha_M = alpha * 10000          %Pricing error in Basis Points

% Transforming into annual basis points       
Alpha_M_ols_hac_LS_C = Annualized_alpha_ols_hac * 10000      %Pricing error in Basis Points
Alpha_M_Talwar_LS_C = Annualized_alpha_talwar * 10000 


%% Part Time Regression Long Short Portfolio Monthly 2003 - 2011
pt_cum_port=cumprod((ls_portfolio(1:108,:))+1);

pt_month = month(1:108,:);
% Graph on performance 2003 - 2011
figure;  hold all;
plot(pt_month, pt_cum_port);
plot(pt_month, pt_cum_mark);

% Part time Sharpe

avg_ret_pt_ls=mean(ls_portfolio(1:108,:));
variance_pt_ls=var(ls_portfolio(1:108,:));
Volatility_pt_ls=std(ls_portfolio(1:108,:));
sharpe_pt_ls=avg_ret_pt_ls/Volatility_pt_ls;

% Regression
pt_mkt_e = mkt_e(1:108,:);
pt_smb = smb(1:108,:);
pt_hml = hml(1:108,:);
pt_wml = rmw(1:108,:);
pt_portfolio_ls = ls_portfolio(1:108,:);

Reg_table=table(pt_mkt_e, pt_smb, pt_hml, pt_wml, pt_portfolio_ls, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'WML', 'Portfolio'} );
Regression_Monthly_PT_LS2003_C = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + WML'); %,...

disp(Regression_Monthly_PT_LS2003_C)


% Regression using Hac
[EstCoeffCov, se, coeff] = hac(Regression_Monthly_PT_LS2003_C);
t_hac = coeff./se;
p_hac = erfc(0.7071*abs(t_hac));

HAC_table_PT_LS2003_C = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'WML'});

disp(HAC_table_PT_LS2003_C)

% Regression using OLS Talwar weighting
Regression_Monthly_Talwar_PT_LS2003_C = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + WML', ...
    'RobustOpts', 'talwar'); 
disp(Regression_Monthly_Talwar_PT_LS2003_C)

% Annualizing
Annualized_alpha_ols_hac = Regression_Monthly_PT_LS2003_C.Coefficients.Estimate(1) *12
Annualized_alpha_talwar = Regression_Monthly_Talwar_PT_LS2003_C.Coefficients.Estimate(1) *12
%Alpha_M = alpha * 10000          %Pricing error in Basis Points

% Transforming into annual basis points       
Alpha_M_ols_hac_PT_LS2003_C = Annualized_alpha_ols_hac * 10000      %Pricing error in Basis Points
Alpha_M_Talwar_PT_LS2003_C = Annualized_alpha_talwar * 10000 



%% Part Time Regression Long Short Portfolio Monthly 2012 - 2021
pt_cum_port=cumprod((ls_portfolio(109:end,:))+1);

pt_month = month(109:end,:);
% Graph on performance 2012 - 2021
figure;  hold all;
plot(pt_month, pt_cum_port);
plot(pt_month, pt_cum_mark);

% Part time Sharpe

avg_ret_pt_ls=mean(ls_portfolio(109:end,:));
variance_pt_ls=var(ls_portfolio(109:end,:));
Volatility_pt_ls=std(ls_portfolio(109:end,:));
sharpe_pt_ls=avg_ret_pt_ls/Volatility_pt_ls;

% Regression
pt_mkt_e = mkt_e(109:end,:);
pt_smb = smb(109:end,:);
pt_hml = hml(109:end,:);
pt_wml = rmw(109:end,:);
pt_portfolio_ls = ls_portfolio(109:end,:);

Reg_table=table(pt_mkt_e, pt_smb, pt_hml, pt_wml, pt_portfolio_ls, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'WML', 'Portfolio'} );
Regression_Monthly_PT_LS2021_C = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + WML'); %,...

disp(Regression_Monthly_PT_LS2021_C)


% Regression using Hac
[EstCoeffCov, se, coeff] = hac(Regression_Monthly_PT_LS2021_C);
t_hac = coeff./se;
p_hac = erfc(0.7071*abs(t_hac));

HAC_table_PT_LS2021_C = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'WML'});

disp(HAC_table_PT_LS2021_C)

% Regression using OLS Talwar weighting
Regression_Monthly_Talwar_PT_LS2021_C = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + WML', ...
    'RobustOpts', 'talwar'); 
disp(Regression_Monthly_Talwar_PT_LS2021_C)

% Annualizing
Annualized_alpha_ols_hac = Regression_Monthly_PT_LS2021_C.Coefficients.Estimate(1) *12
Annualized_alpha_talwar = Regression_Monthly_Talwar_PT_LS2021_C.Coefficients.Estimate(1) *12
%Alpha_M = alpha * 10000          %Pricing error in Basis Points

% Transforming into annual basis points       
Alpha_M_ols_hac_PT_LS2021_C = Annualized_alpha_ols_hac * 10000      %Pricing error in Basis Points
Alpha_M_Talwar_PT_LS2021_C = Annualized_alpha_talwar * 10000 



%% Regression Using Data from Fama-French Long Short Portfolio

cum_port=cumsum((ls_portfolio));
cum_C=cumsum((FF(:,1)));

figure;  hold all;
plot(month, cum_port, 'Linewidth', 2);
plot(month, cum_C, 'Linewidth', 2);


Reg_table=table(C_mkt, C_smb, C_hml, C_wml, ls_portfolio, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'WML', 'Portfolio'} );
Regression_FF_ls_Monthly_C = fitlm(Reg_table,...
    'Portfolio ~ Market + SMB + HML + WML'); %,...

disp(Regression_FF_ls_Monthly_C)



% Regression using Hac
[EstCoeffCov, se, coeff] = hac(Regression_FF_ls_Monthly_C)
t_hac = coeff./se
p_hac = erfc(0.7071*abs(t_hac))

HAC_table_LS_FF_C = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'WML'});

disp(HAC_table_LS_FF_C)

% Regression using OLS Talwar weighting
Regression_Monthly_Talwar_LS_FF_C = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + WML', ...
    'RobustOpts', 'talwar'); 
disp(Regression_Monthly_Talwar_LS_FF_C)

% Annualizing
Annualized_alpha_ols_hac = Regression_FF_ls_Monthly_C.Coefficients.Estimate(1) *12
Annualized_alpha_talwar = Regression_Monthly_Talwar_LS_FF_C.Coefficients.Estimate(1) *12
%Alpha_M = alpha * 10000          %Pricing error in Basis Points

% Transforming into annual basis points       
Alpha_M_ols_hac_LS_FF_C = Annualized_alpha_ols_hac * 10000      %Pricing error in Basis Points
Alpha_M_Talwar_LS_FF_C = Annualized_alpha_talwar * 10000 




%%
% clearvars -except mkt_e portfolio_ex portfolio sharpe sharpe_m Regression_Monthly Alpha_M...
%     Regression_Yearly Alpha_Y sharpe_pt sharpe_pt_mkt Regression_Part_Time_Monthly ...
%     Alpha_PT_M Regression_Part_Time_Yearly Alpha_PT_Y Regression_FF_Monthly...
%     Alpha_FF_M Regression_FF_Yearly Alpha_FF_Y low_portfolio_control low_portfolio ...
%     low_portfolio_ex sharpe_low avg_ret_low avg_ret_m avg_ret Regression_LOW_Monthly...
%     Alpha_LOW_M yearly_low_ex_returns Regression_LOW_Yearly Alpha_LOW_Y sharpe_pt_low ...
%     Regression_LOW_Part_Time_Monthly Alpha_LOW_PT_M Regression_LOW_Part_Time_Yearly ...
%     Alpha_PT_Y Regression_FF_LOW_Monthly Alpha_FF_LOW_M Regression_FF_LOW_Yearly ...
%     Alpha_FF_LOW_Y ls_portfolio avg_ret_ls sharpe_ls Regression_ls_Monthly Alpha_ls_M ...
%     yearly_ls_returns Regression_ls_Yearly Alpha_ls_Y sharpe_pt_ls ... 
%     Regression_ls_Part_Time_Monthly Alpha_ls_PT_M Regression_ls_Part_Time_Yearly...
%     Alpha_ls_PT_Y Regression_FF_ls_Monthly Alpha_FF_ls_M Regression_FF_ls_Yearly ...
%     Alpha_FF_ls_Y
