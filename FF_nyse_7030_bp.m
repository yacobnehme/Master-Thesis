clear;
clc;
clearvars -except returns cap mcap bm op inv rf Factors ESG NYSE FF;
%%
returns=readmatrix("Returns.xlsx", "Sheet", "Returns");
%returns=readmatrix("LogReturns.xlsx", "Sheet", "LogReturns");
mcap=readmatrix("Market Cap.xlsx", "sheet", "Monthly MC");
cap=readmatrix("Market Cap.xlsx", "Sheet", "Yearly MC");
bm=readmatrix("BM.xlsx", "Sheet", "BM");
op=readmatrix("OP.xlsx", "Sheet", "OP");
inve=readmatrix("INV.xlsx", "Sheet", "INV");
mom=readmatrix("Returns2.xlsx", "Sheet","Mom");
mom=mom(:,2:end);

NYSE=readmatrix("NYSE.xlsx", "Sheet","NYSE");
FF=readmatrix("FF2.xlsx", "sheet", "FF");
%FF=readmatrix("FAMA.xlsx", "sheet", "Blad1");



month=readmatrix("Returns.xlsx", "Sheet", "Months");
month=datetime(month, 'ConvertFrom', 'Excel', 'Format', 'MMM-uuuu');


NYSE=logical(NYSE);
FF=FF/100;
% Transforming into Log returns

%Getting the 6th coloumn from FF i.e getting the risk-free (rf)
rf=FF(:,6);

%% Breakpoints

BP_ME=zeros(size(month));
BP_MARKETCAP=zeros(height(cap),1);
BP_BM=zeros(height(bm),2);
BP_OP=zeros(height(op),2);
BP_INV=zeros(height(inve),2);
BP_MOM=zeros(height(mom),2);
for i=1:height(mcap)
    
    mask=(isnan(mcap(i,:)) | isnan(returns(i,:))... 
         | isnan(bm(ceil(i/12),:)) | isnan(op(ceil(i/12),:)) |isnan(inve(ceil(i/12),:)))==0;

    % mask sorterar ut bolag med fullständing data, framförallt pga att breakpoint för mcap
    % ska vara konstant över alla variabler.

    nyse_mask=NYSE;
    nyse_mask(mask == 0)= 0;

    BP_ME(i,:) = prctile(mcap(i,nyse_mask),50);
    BP_MARKETCAP(ceil(i/12),:) = prctile(cap(ceil(i/12),nyse_mask),50);
    BP_BM(ceil(i/12),1) = prctile(bm(ceil(i/12),nyse_mask),30);
    BP_BM(ceil(i/12),2) = prctile(bm(ceil(i/12),nyse_mask),70);
    BP_OP(ceil(i/12),1) = prctile(op(ceil(i/12),nyse_mask),30);
    BP_OP(ceil(i/12),2) = prctile(op(ceil(i/12),nyse_mask),70);
    BP_INV(ceil(i/12),1) = prctile(inve(ceil(i/12),nyse_mask),30);
    BP_INV(ceil(i/12),2) = prctile(inve(ceil(i/12),nyse_mask),70);
    BP_MOM(i,1) = prctile(mom(i,nyse_mask),30);
    BP_MOM(i,2) = prctile(mom(i,nyse_mask),70);
          
end


%% Market returns
[T,N]=size(returns);
w=zeros(T,N); 
market_ret=zeros(T,1);
for t=1:T

    mask = (isnan(returns(t,:)) | isnan(mcap(t,:)))==0;
    w(t,mask)=mcap(t,mask)/sum(mcap(t,mask));  
    %w(t,mask)=market_cap(t,mask)/sum(market_cap(t,mask));  
    %w(t,mask)=1/sum(mask);

    market_ret(t,1)=sum(w(t,mask).*returns(t,mask));
end

%%  BIG VALUE


portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
w=zeros(T,N);
ret_big_value=zeros(size(market_ret));
portfolio_control=zeros(T,4);
for i=1:height(mcap)
    interval = ceil(i/12)*12-11:ceil(i/12)*12;

    portfolio_mask= cap(ceil(i/12),:) >= BP_MARKETCAP(ceil(i/12),:)...
            & bm(ceil(i/12),:) >= BP_BM(ceil(i/12),2)...
            & ~isnan(cap(ceil(i/12),:))...
            & ~isnan(bm(ceil(i/12),:))...
            & all(~isnan(mcap(interval,:)))...
            & all(~isnan(returns(interval,:)));  
  

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

   % portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     

    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask));
   
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(w(i,:) ~= 0);
    portfolio_control(i,3)=sum(portfolio_returns(i,portfolio_mask) ~= 0);
    portfolio_control(i,4)=sum(portfolio_returns(i,portfolio_mask) == 0);
    portfolio_control(i,5)=sum(isnan(portfolio_returns(i,portfolio_mask)));
    portfolio_control(i,6)=sum(isnan(mcap(i,portfolio_mask)));
    
    
    
    %portfolio_returns(isnan(portfolio_returns)) = 0;
    
    ret_big_value(i,1)=sum(w(i,portfolio_mask).*portfolio_returns(i,portfolio_mask));
end


%% BIG NEUTRAL BM

portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
w=zeros(T,N);
ret_big_neutral_bm=zeros(size(market_ret));
portfolio_control=zeros(T,3);
for i=1:height(mcap)
    
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
   interval = ceil(i/12)*12-11:ceil(i/12)*12;

    portfolio_mask= cap(ceil(i/12),:) >= BP_MARKETCAP(ceil(i/12),:)...
            & bm(ceil(i/12),:) >= BP_BM(ceil(i/12),1) ...
            & bm(ceil(i/12),:) < BP_BM(ceil(i/12),2)...
            & ~isnan(cap(ceil(i/12),:))...
            & ~isnan(bm(ceil(i/12),:))...
            & all(~isnan(mcap(interval,:)))...
            & all(~isnan(returns(interval,:)));  
    

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;
    
    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
    
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 

    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(w(i,:) ~= 0);
    portfolio_control(i,3)=sum(portfolio_returns(i,portfolio_mask) ~= 0);
    portfolio_control(i,4)=sum(portfolio_returns(i,portfolio_mask) == 0);
    portfolio_control(i,5)=sum(isnan(portfolio_returns(i,portfolio_mask)));
    portfolio_control(i,6)=sum(isnan(mcap(i,portfolio_mask)));
    
    %portfolio_returns(isnan(portfolio_returns)) = 0;
   
    ret_big_neutral_bm(i,1)=sum(w(i,portfolio_mask).*portfolio_returns(i,portfolio_mask));
end


%% BIG GROWTH

portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
w=zeros(T,N);
ret_big_growth=zeros(size(market_ret));
for i=1:height(mcap)
   
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
    interval = ceil(i/12)*12-11:ceil(i/12)*12;

    portfolio_mask= cap(ceil(i/12),:) >= BP_MARKETCAP(ceil(i/12),:)...
            & bm(ceil(i/12),:) < BP_BM(ceil(i/12),1)...
            & ~isnan(cap(ceil(i/12),:))...
            & ~isnan(bm(ceil(i/12),:))...
            & all(~isnan(mcap(interval,:)))...
            & all(~isnan(returns(interval,:)));  

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;
    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
    
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
       
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(w(i,:) ~= 0);
    portfolio_control(i,3)=sum(portfolio_returns(i,portfolio_mask) ~= 0);
    portfolio_control(i,4)=sum(portfolio_returns(i,portfolio_mask) == 0);


    %portfolio_returns(isnan(portfolio_returns)) = 0;
    ret_big_growth(i,1)=sum(w(i,portfolio_mask).*portfolio_returns(i,portfolio_mask));
end


%% Small Value

portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
w1=zeros(size(mcap));
w=zeros(T,N);
ret_small_value=zeros(size(market_ret));
for i=1:height(mcap)
    
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
    interval = ceil(i/12)*12-11:ceil(i/12)*12;

    portfolio_mask= cap(ceil(i/12),:) < BP_MARKETCAP(ceil(i/12),:)...
        & bm(ceil(i/12),:) >= BP_BM(ceil(i/12),2)...
            & ~isnan(cap(ceil(i/12),:))...
            & ~isnan(bm(ceil(i/12),:))...
            & all(~isnan(mcap(interval,:)))...
            & all(~isnan(returns(interval,:)));  
    

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask));   
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(w(i,:) ~= 0);
    portfolio_control(i,3)=sum(portfolio_returns(i,portfolio_mask) ~= 0);
    portfolio_control(i,4)=sum(portfolio_returns(i,portfolio_mask) == 0);
    


    %portfolio_returns(isnan(portfolio_returns)) = 0;
    
    ret_small_value(i,1)=sum(w(i,portfolio_mask).*portfolio_returns(i,portfolio_mask));

end

%% Small Neutral BM

portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
w=zeros(T,N);
ret_small_neutral_bm=zeros(size(market_ret));

for i=1:height(mcap)
   
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
    interval = ceil(i/12)*12-11:ceil(i/12)*12;

    
    portfolio_mask= cap(ceil(i/12),:) < BP_MARKETCAP(ceil(i/12),:)...
            & bm(ceil(i/12),:) >= BP_BM(ceil(i/12),1)...
            & bm(ceil(i/12),:) < BP_BM(ceil(i/12),2)...
            & ~isnan(cap(ceil(i/12),:))...
            & ~isnan(bm(ceil(i/12),:))...
            & all(~isnan(mcap(interval,:)))...
            & all(~isnan(returns(interval,:)));  

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(w(i,:) ~= 0);
    portfolio_control(i,3)=sum(portfolio_returns(i,portfolio_mask) ~= 0);
    portfolio_control(i,4)=sum(portfolio_returns(i,portfolio_mask) == 0);


    %portfolio_returns(isnan(portfolio_returns)) = 0;

    ret_small_neutral_bm(i,1)=sum(w(i,portfolio_mask).*portfolio_returns(i,portfolio_mask));
    
end



%% Small Growth

portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
ret_small_growth=zeros(size(market_ret));
w=zeros(T,N);
for i=1:height(mcap)
    
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
    interval = ceil(i/12)*12-11:ceil(i/12)*12;

    portfolio_mask = cap(ceil(i/12),:) < BP_MARKETCAP(ceil(i/12),:)...
            & bm(ceil(i/12),:) < BP_BM(ceil(i/12),1)...
            & ~isnan(cap(ceil(i/12),:))...
            & ~isnan(bm(ceil(i/12),:))...
            & all(~isnan(mcap(interval,:)))...
            & all(~isnan(returns(interval,:)));  

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(w(i,:) ~= 0);
    portfolio_control(i,3)=sum(portfolio_returns(i,portfolio_mask) ~= 0);
    portfolio_control(i,4)=sum(portfolio_returns(i,portfolio_mask) == 0);
    portfolio_control(i,5)=sum(isnan(portfolio_returns(i,portfolio_mask)));
    portfolio_control(i,6)=sum(isnan(mcap(i,portfolio_mask)));

    %portfolio_returns(isnan(portfolio_returns)) = 0;

    ret_small_growth(i,1)=sum(w(i,portfolio_mask).*portfolio_returns(i,portfolio_mask));
end

%% Descriptive SMB (big/small value,neutral growth)
mu_size_bm=[mean(ret_big_value), mean(ret_big_neutral_bm), mean(ret_big_growth);...
    mean(ret_small_value), mean(ret_small_neutral_bm), mean(ret_small_growth)];


%% SMB B/M

smb_bm=(1/3)*(ret_small_value + ret_small_neutral_bm + ret_small_growth)...
    - (1/3)*(ret_big_value + ret_big_neutral_bm + ret_big_growth);

%% Big Robust

portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
ret_big_robust=zeros(size(market_ret));
w=zeros(T,N);
for i=1:height(mcap)
    
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
    interval = ceil(i/12)*12-11:ceil(i/12)*12;
    
    portfolio_mask= cap(ceil(i/12),:) >= BP_MARKETCAP(ceil(i/12),:)...
        & op(ceil(i/12),:) >= BP_OP(ceil(i/12),2)...
            & ~isnan(cap(ceil(i/12),:))...
            & ~isnan(op(ceil(i/12),:))...
            & all(~isnan(mcap(interval,:)))...
            & all(~isnan(returns(interval,:)));  

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(w(i,:) ~= 0);
    portfolio_control(i,3)=sum(portfolio_returns(i,portfolio_mask) ~= 0);
    portfolio_control(i,4)=sum(portfolio_returns(i,portfolio_mask) == 0);


    %portfolio_returns(isnan(portfolio_returns)) = 0;
    ret_big_robust(i,1)=sum(w(i,portfolio_mask).*portfolio_returns(i,portfolio_mask));
end


%% Big Neutral OP

portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
ret_big_neutral_op=zeros(size(market_ret));
w=zeros(T,N);
for i=1:height(mcap)
    
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
    interval = ceil(i/12)*12-11:ceil(i/12)*12;
    
    portfolio_mask = cap(ceil(i/12),:) >= BP_MARKETCAP(ceil(i/12),:)...
            & op(ceil(i/12),:) >= BP_OP(ceil(i/12),1)...
            & op(ceil(i/12),:) < BP_OP(ceil(i/12),2)...
            & ~isnan(cap(ceil(i/12),:))...
            & ~isnan(op(ceil(i/12),:))...
            & all(~isnan(mcap(interval,:)))...
            & all(~isnan(returns(interval,:)));  

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(w(i,:) ~= 0);
    portfolio_control(i,3)=sum(portfolio_returns(i,portfolio_mask) ~= 0);
    portfolio_control(i,4)=sum(portfolio_returns(i,portfolio_mask) == 0);


    %portfolio_returns(isnan(portfolio_returns)) = 0;
    ret_big_neutral_op(i,1)=sum(w(i,portfolio_mask).*portfolio_returns(i,portfolio_mask));
end



%% BIG Weak

portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
ret_big_weak=zeros(size(market_ret));
w=zeros(T,N);
for i=1:height(mcap)
   
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
    interval = ceil(i/12)*12-11:ceil(i/12)*12;


    portfolio_mask = cap(ceil(i/12),:) >= BP_MARKETCAP(ceil(i/12),:)...
            & op(ceil(i/12),:) < BP_OP(ceil(i/12),1)...
            & ~isnan(cap(ceil(i/12),:))...
            & ~isnan(op(ceil(i/12),:))...
            & all(~isnan(mcap(interval,:)))...
            & all(~isnan(returns(interval,:)));  

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(w(i,:) ~= 0);
    portfolio_control(i,3)=sum(portfolio_returns(i,portfolio_mask) ~= 0);
    portfolio_control(i,4)=sum(portfolio_returns(i,portfolio_mask) == 0);


    %portfolio_returns(isnan(portfolio_returns)) = 0;
    ret_big_weak(i,1)=sum(w(i,portfolio_mask).*portfolio_returns(i,portfolio_mask));
end



%% Small Robust
portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
ret_small_robust=zeros(size(market_ret));
w=zeros(T,N);
for i=1:height(mcap)
   
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
    interval = ceil(i/12)*12-11:ceil(i/12)*12;

    portfolio_mask= cap(ceil(i/12),:) < BP_MARKETCAP(ceil(i/12),:) ...
        & op(ceil(i/12),:) >= BP_OP(ceil(i/12),2)...
            & ~isnan(cap(ceil(i/12),:))...
            & ~isnan(op(ceil(i/12),:))...
            & all(~isnan(mcap(interval,:)))...
            & all(~isnan(returns(interval,:)));  

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(w(i,:) ~= 0);
    portfolio_control(i,3)=sum(portfolio_returns(i,portfolio_mask) ~= 0);
    portfolio_control(i,4)=sum(portfolio_returns(i,portfolio_mask) == 0);


    %portfolio_returns(isnan(portfolio_returns)) = 0;
    ret_small_robust(i,1)=sum(w(i,portfolio_mask).*portfolio_returns(i,portfolio_mask));
end



%% Small Neutral OP
portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
ret_small_neutral_op=zeros(size(market_ret));
w=zeros(T,N);
for i=1:height(mcap)
    
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
    interval = ceil(i/12)*12-11:ceil(i/12)*12;

    portfolio_mask = cap(ceil(i/12),:) < BP_MARKETCAP(ceil(i/12),:)...
            & op(ceil(i/12),:) >= BP_OP(ceil(i/12),1)...
            & op(ceil(i/12),:) < BP_OP(ceil(i/12),2)...
            & ~isnan(cap(ceil(i/12),:))...
            & ~isnan(op(ceil(i/12),:))...
            & all(~isnan(mcap(interval,:)))...
            & all(~isnan(returns(interval,:)));  

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(w(i,:) ~= 0);
    portfolio_control(i,3)=sum(portfolio_returns(i,portfolio_mask) ~= 0);
    portfolio_control(i,4)=sum(portfolio_returns(i,portfolio_mask) == 0);


    %portfolio_returns(isnan(portfolio_returns)) = 0;
    ret_small_neutral_op(i,1)=sum(w(i,portfolio_mask).*portfolio_returns(i,portfolio_mask));
end



%% Small Weak
portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
ret_small_weak=zeros(size(market_ret));
w=zeros(T,N);

for i=1:height(mcap)
    
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
    interval = ceil(i/12)*12-11:ceil(i/12)*12;  

    portfolio_mask= cap(ceil(i/12),:) < BP_MARKETCAP(ceil(i/12),:)...
            & op(ceil(i/12),:) < BP_OP(ceil(i/12),1)...
            & ~isnan(cap(ceil(i/12),:))...
            & ~isnan(op(ceil(i/12),:))...
            & all(~isnan(mcap(interval,:)))...
            & all(~isnan(returns(interval,:)));  

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(w(i,:) ~= 0);
    portfolio_control(i,3)=sum(portfolio_returns(i,portfolio_mask) ~= 0);
    portfolio_control(i,4)=sum(portfolio_returns(i,portfolio_mask) == 0);


    %portfolio_returns(isnan(portfolio_returns)) = 0;
    ret_small_weak(i,1)=sum(w(i,portfolio_mask).*portfolio_returns(i,portfolio_mask));
end

%% Descriptive SMB (big/small robust,neutral weak)
mu_size_op=[mean(ret_big_robust), mean(ret_big_neutral_op), mean(ret_big_weak);...
    mean(ret_small_robust), mean(ret_small_neutral_op), mean(ret_small_weak)];

%% SMB OP

smb_op=(1/3)*(ret_small_robust + ret_small_neutral_op + ret_small_weak)...
    - (1/3)*(ret_big_robust + ret_big_neutral_op + ret_big_weak);


%% Big Conservative (Low investment-ratio??)

portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
ret_big_conservative=zeros(size(market_ret));
w=zeros(T,N);
for i=1:height(mcap)
    
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
    interval = ceil(i/12)*12-11:ceil(i/12)*12;
    
    portfolio_mask = cap(ceil(i/12),:) >= BP_MARKETCAP(ceil(i/12),:)...
        & inve(ceil(i/12),:) < BP_INV(ceil(i/12),1)...
            & ~isnan(cap(ceil(i/12),:))...
            & ~isnan(inve(ceil(i/12),:))...
            & all(~isnan(mcap(interval,:)))...
            & all(~isnan(returns(interval,:)));  

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(w(i,:) ~= 0);
    portfolio_control(i,3)=sum(portfolio_returns(i,portfolio_mask) ~= 0);
    portfolio_control(i,4)=sum(portfolio_returns(i,portfolio_mask) == 0);


    %portfolio_returns(isnan(portfolio_returns)) = 0;
    ret_big_conservative(i,1)=sum(w(i,portfolio_mask).*portfolio_returns(i,portfolio_mask));
end


%% Big Neutral Inv

portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
ret_big_neutral_inv=zeros(size(market_ret));
w=zeros(T,N);
for i=1:height(mcap)
    
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
    interval = ceil(i/12)*12-11:ceil(i/12)*12;
    
    portfolio_mask = cap(ceil(i/12),:) >= BP_MARKETCAP(ceil(i/12),:)...
            & inve(ceil(i/12),:) >= BP_INV(ceil(i/12),1)...
            & inve(ceil(i/12),:) < BP_INV(ceil(i/12),2)...
            & ~isnan(cap(ceil(i/12),:))...
            & ~isnan(inve(ceil(i/12),:))...
            & all(~isnan(mcap(interval,:)))...
            & all(~isnan(returns(interval,:)));  

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(w(i,:) ~= 0);
    portfolio_control(i,3)=sum(portfolio_returns(i,portfolio_mask) ~= 0);
    portfolio_control(i,4)=sum(portfolio_returns(i,portfolio_mask) == 0);


    %portfolio_returns(isnan(portfolio_returns)) = 0;
    ret_big_neutral_inv(i,1)=sum(w(i,portfolio_mask).*portfolio_returns(i,portfolio_mask));

end


%% Big Agressive (High investment ratio??)

portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
ret_big_agressive=zeros(size(market_ret));
w=zeros(T,N);
for i=1:height(mcap)
   
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
    interval = ceil(i/12)*12-11:ceil(i/12)*12;
    
    portfolio_mask = cap(ceil(i/12),:) >= BP_MARKETCAP(ceil(i/12),:)...
            & inve(ceil(i/12),:) >= BP_INV(ceil(i/12),2)...
            & ~isnan(cap(ceil(i/12),:))...
            & ~isnan(inve(ceil(i/12),:))...
            & all(~isnan(mcap(interval,:)))...
            & all(~isnan(returns(interval,:)));  

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(w(i,:) ~= 0);
    portfolio_control(i,3)=sum(portfolio_returns(i,portfolio_mask) ~= 0);
    portfolio_control(i,4)=sum(portfolio_returns(i,portfolio_mask) == 0);


    portfolio_returns(isnan(portfolio_returns)) = 0;
    ret_big_agressive(i,1)=sum(w(i,:).*portfolio_returns(i,:));
end


%% Small Conservative (Low Investment ratio??))

portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
ret_small_conservative=zeros(size(market_ret));
w=zeros(T,N);
for i=1:height(mcap)
   
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
    interval = ceil(i/12)*12-11:ceil(i/12)*12;

    portfolio_mask= cap(ceil(i/12),:) < BP_MARKETCAP(ceil(i/12),:)...
        & inve(ceil(i/12),:) < BP_INV(ceil(i/12),1)...
            & ~isnan(cap(ceil(i/12),:))...
            & ~isnan(inve(ceil(i/12),:))...
            & all(~isnan(mcap(interval,:)))...
            & all(~isnan(returns(interval,:)));  

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(w(i,:) ~= 0);
    portfolio_control(i,3)=sum(portfolio_returns(i,portfolio_mask) ~= 0);
    portfolio_control(i,4)=sum(portfolio_returns(i,portfolio_mask) == 0);


    %portfolio_returns(isnan(portfolio_returns)) = 0;
    ret_small_conservative(i,1)=sum(w(i,portfolio_mask).*portfolio_returns(i,portfolio_mask));
end


%% Small Neutral Inv

portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
ret_small_neutral_inv=zeros(size(market_ret));
w=zeros(T,N);
for i=1:height(mcap)
    
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
    interval = ceil(i/12)*12-11:ceil(i/12)*12;        

    portfolio_mask = cap(ceil(i/12),:) < BP_MARKETCAP(ceil(i/12),:)...
            & inve(ceil(i/12),:) >= BP_INV(ceil(i/12),1)...
            & inve(ceil(i/12),:) < BP_INV(ceil(i/12),2)...
            & ~isnan(cap(ceil(i/12),:))...
            & ~isnan(inve(ceil(i/12),:))...
            & all(~isnan(mcap(interval,:)))...
            & all(~isnan(returns(interval,:)));  

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(w(i,:) ~= 0);
    portfolio_control(i,3)=sum(portfolio_returns(i,portfolio_mask) ~= 0);
    portfolio_control(i,4)=sum(portfolio_returns(i,portfolio_mask) == 0);


    %portfolio_returns(isnan(portfolio_returns)) = 0;
    ret_small_neutral_inv(i,1)=sum(w(i,portfolio_mask).*portfolio_returns(i,portfolio_mask));
end


%% Small Agressive (High investment ratio??)

portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
ret_small_agressive=zeros(size(market_ret));
w=zeros(T,N);

for i=1:height(mcap)
    
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
    interval = ceil(i/12)*12-11:ceil(i/12)*12;        

    portfolio_mask= cap(ceil(i/12),:) < BP_MARKETCAP(ceil(i/12),:)...
            & inve(ceil(i/12),:) >= BP_INV(ceil(i/12),2)...
            & ~isnan(cap(ceil(i/12),:))...
            & ~isnan(inve(ceil(i/12),:))...
            & all(~isnan(mcap(interval,:)))...
            & all(~isnan(returns(interval,:)));  
        

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(w(i,:) ~= 0);
    portfolio_control(i,3)=sum(portfolio_returns(i,portfolio_mask) ~= 0);
    portfolio_control(i,4)=sum(portfolio_returns(i,portfolio_mask) == 0);


    %portfolio_returns(isnan(portfolio_returns)) = 0;
    ret_small_agressive(i,1)=sum(w(i,portfolio_mask).*portfolio_returns(i,portfolio_mask));
end

%% Descriptive SMB (big/small conservative,neutral agressive)
mu_size_inv=[mean(ret_big_conservative), mean(ret_big_neutral_inv), mean(ret_big_agressive);...
    mean(ret_small_conservative), mean(ret_small_neutral_inv), mean(ret_small_agressive)];


%% SMB INV

smb_inv=(1/3)*(ret_small_conservative + ret_small_neutral_inv + ret_small_agressive)...
    - (1/3)*(ret_big_conservative + ret_big_neutral_inv + ret_big_agressive);
%% SMB

smb=(1/3)*(smb_bm + smb_op + smb_inv);
mu_smb=mean(smb);

%% HML    
hml=(1/2)*(ret_small_value + ret_big_value) - (1/2)*(ret_small_growth + ret_big_growth);

%% RMW

rmw=(1/2)*(ret_small_robust + ret_big_robust) - (1/2)*(ret_small_weak + ret_big_weak);

%% CMA
cma=((1/2)*(ret_small_conservative + ret_big_conservative))...
    - ((1/2)*(ret_small_agressive + ret_big_agressive));

%% Transforming portfolios into log returns
mkt_e=market_ret-rf;
mkt_e=log(mkt_e+1);
smb=log(smb+1);
hml=log(hml+1);
rmw=log(rmw+1);
cma=log(cma+1);

%% Factor Plots
cum_m=cumsum((mkt_e));
cum_smb=cumsum(smb);
cum_hml=cumsum(hml);
cum_rmw=cumsum(rmw);
cum_cma=cumsum(cma);

figure; hold all;

plot(cum_smb, "linewidth", 2);
plot(cum_hml, "linewidth", 2);
plot(cum_rmw, "linewidth", 2);
plot(cum_cma, "linewidth", 2);
legend("smb","hml", "rmw", "cma", "location", "northwest");

%% Factors market, smb, hml, rmw, mca
Factors=[market_ret, smb, hml, rmw, cma, rf];

FF=log(FF+1);

cum_FF=cumsum(FF(:,1));
figure; hold all;
plot(month, cum_m, "linewidth", 2);
plot(month, cum_FF, "linewidth", 2);
legend("Sample Market Returns", "FF Market Returns", "location", "northwest");
xlabel('Time');
ylabel('Cumulative sum of logarithmic return')

cum_FF=cumsum(FF(:,2));
figure; hold all;
plot(month, cum_smb, "linewidth", 2);
plot(month, cum_FF, "linewidth", 2);
legend("Sample SMB", "FF SMB",'Location', 'northwest');
xlabel('Time');
ylabel('Cumulative sum of logarithmic return')
ylim([-0.1 0.8])

cum_FF=cumsum(FF(:,3));
figure; hold all;
plot(month, cum_hml, "linewidth", 2);
plot(month, cum_FF, "linewidth", 2);
legend("Sample HML", "FF HML",'Location', 'northeast');
xlabel('Time');
ylabel('Cumulative sum of logarithmic return')

cum_FF=cumsum(FF(:,4));
figure; hold all;
plot(month, cum_rmw, "linewidth", 2);
plot(month, cum_FF, "linewidth", 2);
legend("Sample RMW", "FF RMW", 'Location', 'northwest');
ylim([-0.3 0.8])
xlabel('Time');
ylabel('Cumulative sum of logarithmic return')

cum_FF=cumsum(FF(:,5));
figure; hold all;
plot(month, cum_cma, "linewidth", 2);
plot(month, cum_FF, "linewidth", 2);
legend("Sample CMA", "FF CMA", 'Location', 'northwest');
xlabel('Time');
ylabel('Cumulative sum of logarithmic return')

%% Big High Mom

portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
ret_big_hmom=zeros(size(market_ret));
w=zeros(T,N);
for i=1:height(mcap)
    
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
       
    

    portfolio_mask= mcap(i,:) >= BP_ME(i,:)...
        & mom(i,:) >= BP_MOM(i,2)...
        & ~isnan(mcap(i,:))...
            & ~isnan(mom(i,:))...
            & ~isnan(mcap(i,:))...
            & ~isnan(returns(i,:));  
        

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(~isnan(portfolio_returns(i,:)));
    portfolio_control(i,3)=sum(w(i,:) ~= 0);


    portfolio_returns(isnan(portfolio_returns)) = 0;
    ret_big_hmom(i,1)=sum(w(i,:).*portfolio_returns(i,:));
end


%% Big Neutral Mom

portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
ret_big_neutral_mom=zeros(size(market_ret));
w=zeros(T,N);
for i=1:height(mcap)
    
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
        
    
        
    portfolio_mask = mcap(i,:) >= BP_ME(i,:)...
            & mom(i,:) >= BP_MOM(i,1)...
            & mom(i,:) < BP_MOM(i,2)...
            & ~isnan(mcap(i,:))...
            & ~isnan(mom(i,:))...
            & ~isnan(mcap(i,:))...
            & ~isnan(returns(i,:)); 
    
     
        

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(~isnan(portfolio_returns(i,:)));
    portfolio_control(i,3)=sum(w(i,:) ~= 0);


    portfolio_returns(isnan(portfolio_returns)) = 0;
    ret_big_neutral_mom(i,1)=sum(w(i,:).*portfolio_returns(i,:));

end



%% BIG Low Mom

portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
ret_big_lmom=zeros(size(market_ret));
w=zeros(T,N);
for i=1:height(mcap)
   
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
        
       
    portfolio_mask = mcap(i,:) >= BP_ME(i,:)...
            & mom(i,:) < BP_MOM(i,1)...
            & ~isnan(mcap(i,:))...
            & ~isnan(mom(i,:))...
            & ~isnan(mcap(i,:))...
            & ~isnan(returns(i,:));
    
        

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(~isnan(portfolio_returns(i,:)));
    portfolio_control(i,3)=sum(w(i,:) ~= 0);


    portfolio_returns(isnan(portfolio_returns)) = 0;
    ret_big_lmom(i,1)=sum(w(i,:).*portfolio_returns(i,:));
end



%% Small High Mom
portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
ret_small_hmom=zeros(size(market_ret));
w=zeros(T,N);
for i=1:height(mcap)
   
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
        
    
        
    portfolio_mask= mcap(i,:) < BP_ME(i,:)...
        & mom(i,:) >= BP_MOM(i,2)...
        & ~isnan(mcap(i,:))...
        & ~isnan(mom(i,:))...
        & ~isnan(mcap(i,:))...
        & ~isnan(returns(i,:));

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(~isnan(portfolio_returns(i,:)));
    portfolio_control(i,3)=sum(w(i,:) ~= 0);


    portfolio_returns(isnan(portfolio_returns)) = 0;
    ret_small_hmom(i,1)=sum(w(i,:).*portfolio_returns(i,:));
end



%% Small Neutral Mom
portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
ret_small_neutral_mom=zeros(size(market_ret));
w=zeros(T,N);
for i=1:height(mcap)
    
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
        
    portfolio_mask = mcap(i,:) < BP_ME(i,:)...
            & mom(i,:) >= BP_MOM(i,1)...
            & mom(i,:) < BP_MOM(i,2)...
        & ~isnan(mcap(i,:))...
        & ~isnan(mom(i,:))...
        & ~isnan(mcap(i,:))...
        & ~isnan(returns(i,:));
    


    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(~isnan(portfolio_returns(i,:)));
    portfolio_control(i,3)=sum(w(i,:) ~= 0);


    portfolio_returns(isnan(portfolio_returns)) = 0;
    ret_small_neutral_mom(i,1)=sum(w(i,:).*portfolio_returns(i,:));
end



%% Small Low Mom
portfolio_returns=NaN(size(returns));
[T,N]=size(portfolio_returns);
ret_small_lmom=zeros(size(market_ret));
w=zeros(T,N);

for i=1:height(mcap)
    
%     start_month=(i-1)*12+1;
%     end_month=start_month+11;
        
        
    portfolio_mask= mcap(i,:) < BP_ME(i,:)...
        & mom(i,:) < BP_MOM(i,1)...
        & ~isnan(mcap(i,:))...
        & ~isnan(mom(i,:))...
        & ~isnan(mcap(i,:))...
        & ~isnan(returns(i,:));

    

    portfolio_returns(i,:)=returns(i,:).*portfolio_mask;

    portfolio_returns(portfolio_returns == 0 & portfolio_mask == 0) = NaN;
     
    w(i,portfolio_mask)=mcap(i,portfolio_mask)/sum(mcap(i,portfolio_mask)); 
    %w(i,portfolio_mask)=market_cap(i,portfolio_mask)/sum(market_cap(i,portfolio_mask)); 
   
    portfolio_control(i,1)=sum(portfolio_mask);
    portfolio_control(i,2)=sum(~isnan(portfolio_returns(i,:)));
    portfolio_control(i,3)=sum(w(i,:) ~= 0);


    portfolio_returns(isnan(portfolio_returns)) = 0;
    ret_small_lmom(i,1)=sum(w(i,:).*portfolio_returns(i,:));
end

%% Descriptive SMB (big/small High Mom,neutral Mom, Low Mom)
mu_size_mom=[mean(ret_big_hmom), mean(ret_big_neutral_mom), mean(ret_big_lmom);...
    mean(ret_small_hmom), mean(ret_small_neutral_mom), mean(ret_small_lmom)];

%% WML

wml=(1/2)*(ret_small_hmom + ret_big_hmom) - (1/2)*(ret_small_lmom + ret_big_lmom);

wml=log(wml+1);
cum_wml=cumsum(wml);
cum_C=cumsum(Carhart(:,4));
figure; hold all;
plot(month,cum_wml, "linewidth", 2);
plot(month,cum_C, "linewidth", 2);
legend("Sample WML", "Carhart WML",'Location', 'northwest');
xlabel('Time');
ylabel('Cumulative sum of logarithmic return')

