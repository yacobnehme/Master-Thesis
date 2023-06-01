clear;
clc;
clearvars -except returns cap mcap smb hml rmw cma rf ESG FF NYSE month market_ret;

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


%% Portfolio Construction High

[T,N]=size(returns);
portfolio_returns=zeros(T,N);
portfolio=zeros(T,1);
w=zeros(size(returns));
portfolio_control=zeros(T,1);
w_control=zeros(T,1);

for i=1:height(returns)
    
    [ESG_sorted, index] = sort(esg(i,:),'descend', MissingPlacement='last');
    
    ESG_quantile = prctile(ESG_sorted, 90 , 2);     %%% KOLLA HÄR %%%
    

    mask = esg(i,:) >= ESG_quantile & ~isnan(returns(i,:)) & ~isnan(mcap(i,:)); 


    portfolio_returns(i,:)=returns(i,:).*mask;
   

    
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
    
    portfolio(i,1)=sum(w(i,mask).*portfolio_returns(i,mask));
end


portfolio_excess = portfolio - rf;
portfolio_ex=log(portfolio_excess +1);


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


%% Full Time-horizon Regression

cum_port=cumsum((portfolio_ex));
cum_mark=cumsum((mkt_e));

figure;  hold all;
plot(month,cum_port);
plot(month,cum_mark);

T=length(month);
X=[mkt_e, smb , hml, rmw, cma];
Y=portfolio_ex;

% Regression using OLS
Reg_table=table(mkt_e, smb, hml, rmw, cma, portfolio_ex, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'RMW', 'CMA', 'Portfolio'} );
Regression_Monthly = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + RMW + CMA'); %,...

disp(Regression_Monthly)

% Test for Homoscedasticity - Engle's Arch Test  
[h, Pvalue, stat, cValue] = archtest(Regression_Monthly.Residuals.Raw, Lags=[1 3]) 
% h=0 -> H0 not rejected = homo, men p= 0.46. Hac for safety
% med 3 laggar, h=1 -> H0 rejected = hetero, p= 0.0169. Hac

% Test for Autocorrelation
[h, Pvalue, stat, cValue] = lbqtest(Regression_Monthly.Residuals.Raw)
% h=0 -> H0 not rejected = ingen autokorrelation, men p=0.59 Hac for safety


% Regression using Hac
[EstCoeffCov, se, coeff] = hac(Regression_Monthly)
t_hac = coeff./se
p_hac = erfc(0.7071*abs(t_hac))




HAC_table = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'RMW', 'CMA'});

disp(HAC_table)


% Regression using OLS Talwar weighting
Regression_Monthly_Talwar = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + RMW + CMA', ...
    'RobustOpts', 'talwar'); 
disp(Regression_Monthly_Talwar)




% Annualizing
Annualized_alpha_ols_hac = Regression_Monthly.Coefficients.Estimate(1) *12
Annualized_alpha_talwar = Regression_Monthly_Talwar.Coefficients.Estimate(1) *12
%Alpha_M = alpha * 10000          %Pricing error in Basis Points

% Transforming into basis points
Alpha_M_ols_hac = Annualized_alpha_ols_hac * 10000          %Pricing error in Basis Points
Alpha_M_Talwar = Annualized_alpha_talwar * 10000          %Pricing error in Basis Points


outlier = find(isoutlier(Regression_Monthly.Residuals.Raw))
Regression_Monthly_Talwar.Residuals.Raw(73)
Regression_Monthly_Talwar.Residuals.Raw(10)

anova(Regression_Monthly)
anova(Regression_Monthly_Talwar)


%% Monthly Regression 2003 - 2011

pt_cum_port=cumsum((portfolio_ex(1:108,:)));
pt_cum_mark=cumsum(mkt_e(1:108,:));

% Graph on performance 2012 - 2021
figure;  hold all;
plot(month(1:108),pt_cum_port);
plot(month(1:108),pt_cum_mark);

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
pt_rmw = rmw(1:108,:);
pt_cma = cma(1:108,:);
pt_portfolio_ex = portfolio_ex(1:108,:);

% Regression using OLS
Reg_table=table(pt_mkt_e, pt_smb, pt_hml, pt_rmw, pt_cma, pt_portfolio_ex, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'RMW', 'CMA', 'Portfolio'} );
Regression_Part_Time_Monthly_2003 = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + RMW + CMA');

disp(Regression_Part_Time_Monthly_2003)

% Regression using OLS w HAC standard errors
[EstCoeffCov, se, coeff] = hac(Regression_Part_Time_Monthly_2003);
t_hac = coeff./se;
p_hac = erfc(0.7071*abs(t_hac));

HAC_table_Part_Time2003 = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'RMW', 'CMA'});
disp(HAC_table_Part_Time2003)

% Regression using OLS talwAR
Regression_Part_Time_Monthly_Talwar_2003 = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + RMW + CMA', ...
    'RobustOpts', 'talwar'); 
disp(Regression_Part_Time_Monthly_Talwar_2003)

outlier = find(isoutlier(Regression_Part_Time_Monthly_Talwar_2003.Residuals.Raw))

% Annualizing
Annualized_alpha_ols_hac_pt = Regression_Part_Time_Monthly_2003.Coefficients.Estimate(1) *12
Annualized_alpha_talwar_pt = Regression_Part_Time_Monthly_Talwar_2003.Coefficients.Estimate(1) *12

% Transforming into annual basis points       
Alpha_M_ols_hac_pt2003 = Annualized_alpha_ols_hac_pt * 10000          %Pricing error in Basis Points
Alpha_M_Talwar_pt2003 = Annualized_alpha_talwar_pt * 10000 


%% Monthly Regression 2012 - 2021

pt_cum_port=cumsum((portfolio_ex(109:end,:)));
pt_cum_mark=cumsum(mkt_e(109:end,:));

% Graph on performance 2012 - 2021
figure;  hold all;
plot(month(109:end,:), pt_cum_port);
plot(month(109:end,:), pt_cum_mark);

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
pt_rmw = rmw(109:end,:);
pt_cma = cma(109:end,:);
pt_portfolio_ex = portfolio_ex(109:end,:);


% Regression using OLS
Reg_table=table(pt_mkt_e, pt_smb, pt_hml, pt_rmw, pt_cma, pt_portfolio_ex, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'RMW', 'CMA', 'Portfolio'} );
Regression_Part_Time_Monthly_2021 = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + RMW + CMA');

disp(Regression_Part_Time_Monthly_2021)

% Regression using OLS w HAC standard errors
[EstCoeffCov, se, coeff] = hac(Regression_Part_Time_Monthly_2021);
t_hac = coeff./se;
p_hac = erfc(0.7071*abs(t_hac));

HAC_table_Part_Time_2021 = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'RMW', 'CMA'});
disp(HAC_table_Part_Time_2021)

% Regression using OLS talwAR
Regression_Part_Time_Monthly_Talwar_2021 = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + RMW + CMA', ...
    'RobustOpts', 'talwar'); 
disp(Regression_Part_Time_Monthly_Talwar_2021)

% Annualizing
Annualized_alpha_ols_hac_pt = Regression_Part_Time_Monthly_2021.Coefficients.Estimate(1) *12;
Annualized_alpha_talwar_pt = Regression_Part_Time_Monthly_Talwar_2021.Coefficients.Estimate(1) *12;

% Transforming into annual basis points       
Alpha_M_ols_hac_pt2021 = Annualized_alpha_ols_hac_pt * 10000          %Pricing error in Basis Points
Alpha_M_Talwar_pt2021 = Annualized_alpha_talwar_pt * 10000 


%% Regression Using Data from Fama-French

cum_port=cumsum((portfolio_ex));
cum_FF=cumsum((FF(:,1)));

figure;  hold all;
plot(cum_port);
plot(cum_FF);

FF_mkt = FF(:,1);
FF_smb = FF(:,2);
FF_hml = FF(:,3);
FF_rmw = FF(:,4);
FF_cma = FF(:,5);


Reg_table=table(FF_mkt, FF_smb, FF_hml, FF_rmw, FF_cma, portfolio_ex, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'RMW', 'CMA', 'Portfolio'} );
Regression_FF_Monthly = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + RMW + CMA');

disp(Regression_FF_Monthly)

% Test for Homoscedasticity - Engle's Arch Test  
[h, Pvalue, stat, cValue] = archtest(Regression_FF_Monthly.Residuals.Raw, Lags=[1 4]) 
% h=0 -> H0 not rejected = homo, men p= 0.73. Hac for safety
% med 3 laggar, h=0 -> H0 not rejected = homo, p= 0.4565. Hac

% Test for Autocorrelation
[h, Pvalue, stat, cValue] = lbqtest(Regression_FF_Monthly.Residuals.Raw)
% h=0 -> H0 not rejected = ingen autokorrelation, men p=0.2658 Hac for safety

[EstCoeffCov, se, coeff] = hac(Regression_FF_Monthly)
t_hac = coeff./se
p_hac = erfc(0.7071*abs(t_hac))

HAC_table_FF = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'RMW', 'CMA'});

disp(HAC_table_FF)

Regression_FF_Monthly_Talwar = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + RMW + CMA', ...
    'RobustOpts', 'talwar'); 
disp(Regression_FF_Monthly_Talwar)


Annualized_alpha_ols_hac = Regression_FF_Monthly.Coefficients.Estimate(1) *12
Annualized_alpha_talwar = Regression_FF_Monthly_Talwar.Coefficients.Estimate(1) *12
%Alpha_M = alpha * 10000          %Pricing error in Basis Points

% Transforming into annual basis points       
Alpha_M_ols_hac_FF = Annualized_alpha_ols_hac_pt * 10000      %Pricing error in Basis Points
Alpha_M_Talwar_FF = Annualized_alpha_talwar_pt * 10000 

outliers=find(isoutlier(Regression_FF_Monthly.Residuals.Raw))
outliers=find(isoutlier(Regression_FF_Monthly_Talwar.Residuals.Raw))

anova(Regression_FF_Monthly)
anova(Regression_FF_Monthly_Talwar)


%% Low portfolio

[T,N]=size(returns);
low_portfolio_returns=zeros(T,N);
low_portfolio=zeros(T,1);
w=zeros(size(returns));
low_portfolio_control=zeros(T,1);
w_control=zeros(T,1);
for i=1:height(returns)
    
    [ESG_sorted, index] = sort(esg(i,:),'descend', MissingPlacement='last');
    
    ESG_quantile = prctile(ESG_sorted, 10 , 2);         %%% KOLLA HÄR %%%

    mask = (esg(i,:) <= ESG_quantile) & ~isnan(returns(i,:)) & ~isnan(mcap(i,:));
    
    low_portfolio_control(i,:)=sum(mask);


    low_portfolio_returns(i,:)=returns(i,:).*mask;
   
    

    
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
    
    low_portfolio(i,1)=sum(w(i,mask).*low_portfolio_returns(i,mask));
end

%% Low portfolio Graphs
low_portfolio_excess = low_portfolio-rf;
low_portfolio_ex = log(low_portfolio_excess+1);

%low_cum_port=cumprod((low_portfolio_ex)+1);
%cum_mark=cumprod((mkt_e)+1);


figure;  hold all;
plot(low_cum_port);
plot(cum_mark);

% Low portfolio Sharpe
avg_ret_low=mean(low_portfolio_ex);
variance_low=var(low_portfolio_ex);
Volatility_low=std(low_portfolio_ex);
sharpe_low=avg_ret_low/Volatility_low;

%% Low portfolio Full Time-horizon Regression

% Regression using OLS
Reg_table=table(mkt_e, smb, hml, rmw, cma, low_portfolio_ex, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'RMW', 'CMA', 'Portfolio'} );
Regression_Monthly_Low = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + RMW + CMA'); %,...

disp(Regression_Monthly_Low)

% Test for Homoscedasticity - Engle's Arch Test  
[h, Pvalue, stat, cValue] = archtest(Regression_Monthly_Low.Residuals.Raw, Lags=[1 3]) 
% h=0 -> H0 not rejected = homo, men p= 0.32. Hac for safety
% med 3 laggar, h=0 -> H0 not rejected = homo, p= 0.52. Hac

% Test for Autocorrelation
[h, Pvalue, stat, cValue] = lbqtest(Regression_Monthly_Low.Residuals.Raw)
% h=0 -> H0 not rejected = ingen autokorrelation, men p=0.54 Hac for safety

% Regression using Hac
[EstCoeffCov, se, coeff] = hac(Regression_Monthly_Low)
t_hac = coeff./se
p_hac = erfc(0.7071*abs(t_hac))

HAC_table_Low = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'RMW', 'CMA'});

disp(HAC_table_Low)

% Regression using OLS Talwar weighting
Regression_Monthly_Talwar_Low = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + RMW + CMA', ...
    'RobustOpts', 'talwar'); 
disp(Regression_Monthly_Talwar_Low)

% Annualizing
Annualized_alpha_ols_hac = Regression_Monthly_Low.Coefficients.Estimate(1) *12
Annualized_alpha_talwar = Regression_Monthly_Talwar_Low.Coefficients.Estimate(1) *12
%Alpha_M = alpha * 10000          %Pricing error in Basis Points

% Transforming into annual basis points       
Alpha_M_ols_hac_Low = Annualized_alpha_ols_hac * 10000      %Pricing error in Basis Points
Alpha_M_Talwar_Low = Annualized_alpha_talwar * 10000 




%% Monthly Part Time regression Low Portfolio 2003 - 2021
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
pt_rmw = rmw(1:108,:);
pt_cma = cma(1:108,:);
pt_portfolio_ex = low_portfolio_ex(1:108,:);

% Regression using OLS
Reg_table=table(pt_mkt_e, pt_smb, pt_hml, pt_rmw, pt_cma, pt_portfolio_ex, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'RMW', 'CMA', 'Portfolio'} );
Regression_Part_Time_Monthly_Low_2003 = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + RMW + CMA');

disp(Regression_Part_Time_Monthly_Low_2003)

% Regression using OLS w HAC standard errors
[EstCoeffCov, se, coeff] = hac(Regression_Part_Time_Monthly_Low_2003);
t_hac = coeff./se;
p_hac = erfc(0.7071*abs(t_hac));

HAC_table_Part_Time_Low2003 = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'RMW', 'CMA'});
disp(HAC_table_Part_Time_Low2003)

% Regression using OLS talwAR
Regression_Part_Time_Monthly_Low_Talwar_2003 = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + RMW + CMA', ...
    'RobustOpts', 'talwar'); 
disp(Regression_Part_Time_Monthly_Low_Talwar_2003)

% Annualizing
Annualized_alpha_ols_hac_pt = Regression_Part_Time_Monthly_Low_2003.Coefficients.Estimate(1) *12
Annualized_alpha_talwar_pt = Regression_Part_Time_Monthly_Low_Talwar_2003.Coefficients.Estimate(1) *12

% Transforming into annual basis points       
Alpha_M_ols_hac_pt_Low2003 = Annualized_alpha_ols_hac_pt * 10000          %Pricing error in Basis Points
Alpha_M_Talwar_pt_Low2003 = Annualized_alpha_talwar_pt * 10000 


%% Monthly Part Time regression Low Portfolio 2012-2021
pt_cum_port=cumprod((low_portfolio_ex(109:end,:))+1);
pt_cum_mark=cumprod(mkt_e(109:end,:)+1);

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
pt_rmw = rmw(109:end,:);
pt_cma = cma(109:end,:);
pt_portfolio_ex = low_portfolio_ex(109:end,:);

% Regression using OLS
Reg_table=table(pt_mkt_e, pt_smb, pt_hml, pt_rmw, pt_cma, pt_portfolio_ex, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'RMW', 'CMA', 'Portfolio'} );
Regression_Part_Time_Monthly_Low_2021 = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + RMW + CMA');

disp(Regression_Part_Time_Monthly_Low_2021)

% Regression using OLS w HAC standard errors
[EstCoeffCov, se, coeff] = hac(Regression_Part_Time_Monthly_Low_2021);
t_hac = coeff./se;
p_hac = erfc(0.7071*abs(t_hac));

HAC_table_Part_Time_Low2021 = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'RMW', 'CMA'});
disp(HAC_table_Part_Time_Low2021)

% Regression using OLS talwAR
Regression_Part_Time_Monthly_Low_Talwar_2021 = fitlm(Reg_table,...
    'Portfolio ~ Market + SMB + HML + RMW + CMA', 'RobustOpts', 'talwar'); 
disp(Regression_Part_Time_Monthly_Low_Talwar_2021)

% Annualizing
Annualized_alpha_ols_hac_pt = Regression_Part_Time_Monthly_Low_2021.Coefficients.Estimate(1) *12
Annualized_alpha_talwar_pt = Regression_Part_Time_Monthly_Low_Talwar_2021.Coefficients.Estimate(1) *12

% Transforming into annual basis points       
Alpha_M_ols_hac_pt_Low2021 = Annualized_alpha_ols_hac_pt * 10000          %Pricing error in Basis Points
Alpha_M_Talwar_pt_Low2021 = Annualized_alpha_talwar_pt * 10000 





%% Regression Using Data from Fama-French Low Portfolio

cum_port=cumprod((low_portfolio_ex)+1);
cum_FF=cumprod((FF(:,1))+1);

figure;  hold all;
plot(cum_port);
plot(cum_FF);


% Regression using OLS
Reg_table=table(FF_mkt, FF_smb, FF_hml, FF_rmw, FF_cma, low_portfolio_ex, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'RMW', 'CMA', 'Portfolio'} );
Regression_Monthly_Low_FF = fitlm(Reg_table,...
    'Portfolio ~ Market + SMB + HML + RMW + CMA'); %,...

disp(Regression_Monthly_Low_FF)

% Test for Homoscedasticity - Engle's Arch Test  
[h, Pvalue, stat, cValue] = archtest(Regression_Monthly_Low_FF.Residuals.Raw, Lags=[1 3]) 
% h=0 -> H0 not rejected = homo, men p= 0.91. Hac for safety
% med 3 laggar, h=0 -> H0 not rejected = homo, p= 0.60. Hac

% Test for Autocorrelation
[h, Pvalue, stat, cValue] = lbqtest(Regression_Monthly_Low_FF.Residuals.Raw)
% h=0 -> H0 not rejected = ingen autokorrelation, men p=0.41 Hac for safety

% Regression using Hac
[EstCoeffCov, se, coeff] = hac(Regression_Monthly_Low_FF)
t_hac = coeff./se
p_hac = erfc(0.7071*abs(t_hac))

HAC_table_Low_FF = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'RMW', 'CMA'});

disp(HAC_table_Low_FF)

% Regression using OLS Talwar weighting
Regression_Monthly_Talwar_Low_FF = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + RMW + CMA', ...
    'RobustOpts', 'talwar'); 
disp(Regression_Monthly_Talwar_Low_FF)

% Annualizing
Annualized_alpha_ols_hac = Regression_Monthly_Low_FF.Coefficients.Estimate(1) *12
Annualized_alpha_talwar = Regression_Monthly_Talwar_Low_FF.Coefficients.Estimate(1) *12
%Alpha_M = alpha * 10000          %Pricing error in Basis Points

% Transforming into annual basis points       
Alpha_M_ols_hac_Low_FF = Annualized_alpha_ols_hac * 10000      %Pricing error in Basis Points
Alpha_M_Talwar_Low_FF = Annualized_alpha_talwar * 10000 


%% Long - Short Portfolio

ls_portfolio = portfolio - low_portfolio; 

cum_ls_port=cumprod((ls_portfolio)+1);
cum_mark=cumprod((mkt_e)+1);

figure;  hold all;
plot(month, cum_ls_port, 'Linewidth', 2);
plot(month, cum_mark, 'Linewidth', 2);

% Long - Short Sharpe

avg_ret_ls=mean(ls_portfolio);
variance_ls=var(ls_portfolio);
Volatility_ls=std(ls_portfolio);
sharpe_ls=avg_ret_ls/Volatility_ls;

%% Long Short Monthly Regression

% Regression using OLS
Reg_table=table(mkt_e, smb, hml, rmw, cma, ls_portfolio, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'RMW', 'CMA', 'Portfolio'} );
Regression_Monthly_LS = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + RMW + CMA'); %,...

disp(Regression_Monthly_LS)

% Test for Homoscedasticity - Engle's Arch Test  
[h, Pvalue, stat, cValue] = archtest(Regression_Monthly_LS.Residuals.Raw, Lags=[1 3]); 
% h=0 -> H0 not rejected = homo, men p= 0.30. Hac for safety
% med 3 laggar, h=0 -> H0 not rejected = homo, p= 0.72. Hac

% Test for Autocorrelation
[h, Pvalue, stat, cValue] = lbqtest(Regression_Monthly_LS.Residuals.Raw);
% h=0 -> H0 not rejected = ingen autokorrelation, men p=0.50 Hac for safety

% Regression using Hac
[EstCoeffCov, se, coeff] = hac(Regression_Monthly_LS);
t_hac = coeff./se;
p_hac = erfc(0.7071*abs(t_hac));

HAC_table_LS = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'RMW', 'CMA'});

disp(HAC_table_LS)

% Regression using OLS Talwar weighting
Regression_Monthly_Talwar_LS = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + RMW + CMA', ...
    'RobustOpts', 'talwar'); 
disp(Regression_Monthly_Talwar_LS)

% Annualizing
Annualized_alpha_ols_hac = Regression_Monthly_LS.Coefficients.Estimate(1) *12
Annualized_alpha_talwar = Regression_Monthly_Talwar_LS.Coefficients.Estimate(1) *12
%Alpha_M = alpha * 10000          %Pricing error in Basis Points

% Transforming into annual basis points       
Alpha_M_ols_hac_LS = Annualized_alpha_ols_hac * 10000      %Pricing error in Basis Points
Alpha_M_Talwar_LS = Annualized_alpha_talwar * 10000 




%% Part Time Regression Long Short Portfolio Monthly 2003 - 2011
pt_cum_port=cumprod((ls_portfolio(1:108,:))+1);
pt_cum_mark=cumprod(mkt_e(1:108,:)+1);

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
pt_rmw = rmw(1:108,:);
pt_cma = cma(1:108,:);
pt_portfolio_ls = ls_portfolio(1:108,:);

Reg_table=table(pt_mkt_e, pt_smb, pt_hml, pt_rmw, pt_cma, pt_portfolio_ls, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'RMW', 'CMA', 'Portfolio'} );

Regression_Monthly_PT_LS2003 = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + RMW + CMA'); %,...

disp(Regression_Monthly_PT_LS2003)


% Regression using Hac
[EstCoeffCov, se, coeff] = hac(Regression_Monthly_PT_LS2003);
t_hac = coeff./se;
p_hac = erfc(0.7071*abs(t_hac));

HAC_table_PT_LS2003 = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'RMW', 'CMA'});

disp(HAC_table_PT_LS2003)

% Regression using OLS Talwar weighting
Regression_Monthly_Talwar_PT_LS2003 = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + RMW + CMA', ...
    'RobustOpts', 'talwar'); 
disp(Regression_Monthly_Talwar_PT_LS2003)

% Annualizing
Annualized_alpha_ols_hac = Regression_Monthly_PT_LS2003.Coefficients.Estimate(1) *12
Annualized_alpha_talwar = Regression_Monthly_Talwar_PT_LS2003.Coefficients.Estimate(1) *12
%Alpha_M = alpha * 10000          %Pricing error in Basis Points

% Transforming into annual basis points       
Alpha_M_ols_hac_PT_LS2003 = Annualized_alpha_ols_hac * 10000      %Pricing error in Basis Points
Alpha_M_Talwar_PT_LS2003 = Annualized_alpha_talwar * 10000 



%% Part Time Regression Long Short Portfolio Monthly 2012-2021
pt_cum_port=cumprod((ls_portfolio(109:end,:))+1);
pt_cum_mark=cumprod(mkt_e(109:end,:)+1);

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
pt_rmw = rmw(109:end,:);
pt_cma = cma(109:end,:);
pt_portfolio_ls = ls_portfolio(109:end,:);

Reg_table=table(pt_mkt_e, pt_smb, pt_hml, pt_rmw, pt_cma, pt_portfolio_ls, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'RMW', 'CMA', 'Portfolio'} );

Regression_Monthly_PT_LS2021 = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + RMW + CMA'); %,...

disp(Regression_Monthly_PT_LS2021)


% Regression using Hac
[EstCoeffCov, se, coeff] = hac(Regression_Monthly_PT_LS2021);
t_hac = coeff./se;
p_hac = erfc(0.7071*abs(t_hac));

HAC_table_PT_LS2021 = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'RMW', 'CMA'});

disp(HAC_table_PT_LS2021)

% Regression using OLS Talwar weighting
Regression_Monthly_Talwar_PT_LS2021 = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + RMW + CMA', ...
    'RobustOpts', 'talwar'); 
disp(Regression_Monthly_Talwar_PT_LS2021)

% Annualizing
Annualized_alpha_ols_hac = Regression_Monthly_PT_LS2021.Coefficients.Estimate(1) *12
Annualized_alpha_talwar = Regression_Monthly_Talwar_PT_LS2021.Coefficients.Estimate(1) *12
%Alpha_M = alpha * 10000          %Pricing error in Basis Points

% Transforming into annual basis points       
Alpha_M_ols_hac_PT_LS2021 = Annualized_alpha_ols_hac * 10000      %Pricing error in Basis Points
Alpha_M_Talwar_PT_LS2021 = Annualized_alpha_talwar * 10000 




%% Regression Using Data from Fama-French Long Short Portfolio

cum_port=cumprod((ls_portfolio)+1);           
cum_FF=cumprod((FF(:,1))+1);

figure;  hold all;
plot(month, cum_port, 'Linewidth', 2);
plot(month, cum_FF, 'Linewidth', 2);



% Regression using OLS
Reg_table=table(FF_mkt, FF_smb, FF_hml, FF_rmw, FF_cma, ls_portfolio, 'VariableNames', {'Market', ...
    'SMB', 'HML', 'RMW', 'CMA', 'Portfolio'} );    
Regression_FF_ls_Monthly = fitlm(Reg_table,...
    'Portfolio ~ Market + SMB + HML + RMW + CMA'); %,...

disp(Regression_FF_ls_Monthly)



% Regression using Hac
[EstCoeffCov, se, coeff] = hac(Regression_FF_ls_Monthly)
t_hac = coeff./se
p_hac = erfc(0.7071*abs(t_hac))

HAC_table_LS_FF = table(coeff, se, t_hac, p_hac, 'VariableNames', {'Estimates_hac' ...
    'se_hac', 'tStat_hac', 'pValue_hac'}, 'RowNames', {'alpha', 'Market', 'SMB', 'HML',...
    'RMW', 'CMA'});

disp(HAC_table_LS_FF)

% Regression using OLS Talwar weighting
Regression_Monthly_Talwar_LS_FF = fitlm(Reg_table, 'Portfolio ~ Market + SMB + HML + RMW + CMA', ...
    'RobustOpts', 'talwar'); 
disp(Regression_Monthly_Talwar_LS_FF)

% Annualizing
Annualized_alpha_ols_hac = Regression_FF_ls_Monthly.Coefficients.Estimate(1) *12
Annualized_alpha_talwar = Regression_Monthly_Talwar_LS_FF.Coefficients.Estimate(1) *12
%Alpha_M = alpha * 10000          %Pricing error in Basis Points

% Transforming into annual basis points       
Alpha_M_ols_hac_LS_FF = Annualized_alpha_ols_hac * 10000      %Pricing error in Basis Points
Alpha_M_Talwar_LS_FF = Annualized_alpha_talwar * 10000 





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
