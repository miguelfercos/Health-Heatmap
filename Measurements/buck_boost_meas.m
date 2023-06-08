%Buck
voutbuck=[99970 100492 99007 100455 99171 100067 100419 100700 100424]/1000
ioutbuck=5004 *ones(1,9)/1000
vinbuck=[149239 154088 144075 139090 133961 129063 123883 119047 113969]/1000
iinbuck=[3449 3360 3538 3716 3806 3985 4163 4342 4521]/1000
d1buck=[6774 6595 6953 7311 7491 7850 8200 8566 8925]/100

poutbuck=voutbuck.*ioutbuck;
pinbuck=vinbuck.*iinbuck;
lossesbuck=pinbuck-poutbuck;
figure
plot(d1buck,lossesbuck, '*')
xlabel('d1')
ylabel('losses (W)')
%%
%Boost
voutboost=[99588 99611 99093 99089 100737 99828 97399 97618 99756 99724]/1000
ioutboost=5003 *ones(1,10)/1000
vinboost=[78815 73632 83869 88912 93918 58568 63524 68711 73639 78858]/1000
iinboost=[6536 7026 6111 5738 5513 9089 8036 7397 7027 6537]/1000
d2boost=[2258 2796 1720 1183 824 4402 3692 3154 2796 2258]/100

poutboost=voutboost.*ioutboost;
pinboost=vinboost.*iinboost;
lossesboost=pinboost-poutboost;
figure
plot(d2boost,lossesboost, '*')
xlabel('d2')
ylabel('losses (W)')

%% total

% vout=[99970 100492 99007 100455 99171 100067 100419 100700 100424 99588 99611 99093 99089 100737 99828 97399 97618 99756 99724]/1000
% iout=[5004 5004 5004 5004 5004 5004 5004 5004 5004 5003 5003 5003 5003 5003 5003 5003 5003 5003 5003 ]/1000
% vin=[149239 154088 144075 139090 133961 129063 123883 119047 113969 78815 73632 83869 88912 93918 58568 63524 68711 73639 78858]/1000
% iin=[3449 3360 3538 3716 3806 3985 4163 4342 4521 6536 7026 6111 5738 5513 9089 8036 7397 7027 6537]/1000
% d1=[6774 6595 6953 7311 7491 7850 8200 8566 8925 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 ]/100
% d2=[10000 10000 10000 10000 10000 10000 10000 10000 10000 2258 2796 1720 1183 824 4402 3692 3154 2796 2258]/100
% 
% pout=vout.*iout;
% pin=vin.*iin;
% losses=pin-pout

%%
d1=linspace(0.05,0.95)
d2=1-(d1*80/100)
figure
plot(d1,d2)
hold on
d1=linspace(0.05,0.95)
d2=1-(d1*150/100)
plot(d1,d2)
hold off
%%
CSV=readtable("Bbheatmap.csv",'Delimiter',';');
vout=CSV.vout;
iout=CSV.iout;
vin=CSV.vin;
iin=CSV.iin;%CSV(:,['iin']); outputs a table, not an array
d1=CSV.d1;
d2=CSV.d2;
pin=vin.*iin;%CSV.pin;
pout=iout.*vout;%CSV.iout;
losses=pin-pout;%CSV.losses;


figure
losses_rearr=reshape (losses,[size(unique(d2),1),size(unique(d1),1)])
contourf(unique(d1),unique(d2),losses_rearr-6,75,'LineColor','none')
colormap(hot) %flipud(hot) para representear colores inversos
colorbar;
d1b=linspace(min(unique(d1)),max(unique(d1)));
d2b=1-(d1b*80/100);
hold on
plot(d1b,d2b,'k','Linewidth',2)
d1b=linspace(min(unique(d1)),max(unique(d1)));
d2b=1-(d1b*150/100);
plot(d1b,d2b,'k','Linewidth',2)
hold off
axis([min(unique(d1)) max(unique(d1)) min(unique(d2)),max(unique(d2))])
% f = figure;
%  ax = axes('Parent',f);
%  h = surf(unique(d1),unique(d2),losses_rearr,'Parent',ax);
%  set(h, 'edgecolor','none');
%  %view(ax,[0,90]);
%  colormap(jet);
%  colorbar;
caxis([6.9 30])