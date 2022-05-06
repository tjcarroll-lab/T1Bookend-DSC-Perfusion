function [ sBV, minBV, minvalues, minind ] = waterxsimu( T1blood, T1ev, tblood, fblood, dR1, oBV )
%waterxsimu Summary of this function goes here
%   Detailed explanation goes here
%   T1blood = T1 of blood
%   T1ev = T1 of extravascular region (tissue)
%   tblood = water exchange rate in blood (time it resides before moving)
%   fblood = percent volume of blood
%
% Author: Yong Ik Jeong
% Date: 2017-09-02

if nargin < 6
    oBV = [];
    if nargin < 5
        dR1 = 0.5:0.5:10;
    end
else
    if length(dR1) ~= length(oBV)
        error('dR1 and oBV must be same length');
    end
end

fprintf('Calculating T1 pre...');
aT1evpre = zeros(length(T1blood),length(tblood),length(fblood));
for hh = 1:length(T1blood)
    for ii = 1:length(tblood)
        for jj = 1:length(fblood)
            fev = 1-fblood(jj);
            tev = fev*tblood(ii)/fblood(jj);
            C1 = .5*(1/T1blood(hh) + 1/T1ev + 1/tblood(ii) + 1/tev);
            C2 = .5*sqrt((1/T1blood(hh) - 1/T1ev + 1/tblood(ii) - 1/tev)^2 + 4/(tblood(ii)*tev));
            efblood = .5 - .25*((fev - fblood(jj))*(1/T1blood(hh) - 1/T1ev) + 1/tblood(ii) + 1/tev)/C2;
            efev = 1 - efblood;
            eT1ev = 1/(C1 - C2);
            eT1blood = 1/(C1 + C2);
            
            t = 0:0.1:2;
            S = efblood*(1 - 2*exp(-t/eT1blood)) + efev*(1 - 2*exp(-t/eT1ev));
            
            ft = fittype('1 - 2*exp(-x/a)');
            f0 = fit(t',S',ft,'startpoint',T1ev);
            aT1evpre(hh,ii,jj) = f0.a;
        end
    end
end
fprintf('Done\n');

fprintf('Calculating T1 post...');
%dR1 = 0.5:0.5:10;
T1bloodpost = 1./(dR1 + 1./T1blood);
%T1bloodpost = 200/1000;
aT1evpost = zeros(length(T1blood),length(tblood),length(fblood),length(T1bloodpost));
for hh = 1:length(T1blood)
    T1bloodpost = 1./(dR1 + 1./T1blood(hh));
    for ii = 1:length(tblood)
        for jj = 1:length(fblood)
            for kk = 1:length(T1bloodpost)
                fev = 1-fblood(jj);
                tev = fev*tblood(ii)/fblood(jj);
                C1 = .5*(1/T1bloodpost(kk) + 1/T1ev + 1/tblood(ii) + 1/tev);
                C2 = .5*sqrt((1/T1bloodpost(kk) - 1/T1ev + 1/tblood(ii) - 1/tev)^2 + 4/(tblood(ii)*tev));
                efblood = .5 - .25*((fev - fblood(jj))*(1/T1bloodpost(kk) - 1/ T1ev) + 1/tblood(ii) + 1/tev)/C2;
                efev = 1 - efblood;
                eT1ev = 1/(C1 - C2);
                eT1blood = 1/(C1 + C2);
                
                t = 0:0.1:2;
                S = efblood*(1 - 2*exp(-t/eT1blood)) + efev*(1 - 2*exp(-t/eT1ev));
                
                ft = fittype('1 - 2*exp(-x/a)');
                f0 = fit(t',S',ft,'startpoint',T1ev);
                aT1evpost(hh,ii,jj,kk) = f0.a;
            end
        end
    end
end
fprintf('Done\n');

fprintf('Calculating BV and least square min...');
currmin = Inf;
minind = [];
sBV = zeros(length(T1blood),length(tblood),length(fblood),length(T1bloodpost));
for hh = 1:length(T1blood)
    T1bloodpost = 1./(dR1 + 1./T1blood(hh));
    for ii = 1:length(tblood)
        for jj = 1:length(fblood)
            tmpmin = 0;
            for kk = 1:length(T1bloodpost)
                sBV(hh,ii,jj,kk) = (1/aT1evpost(hh,ii,jj,kk) - 1/aT1evpre(hh,ii,jj))/(1/T1bloodpost(kk) - 1/T1blood(hh));
                if ~isempty(oBV)
                    tmpmin = tmpmin + (sBV(hh,ii,jj,kk) - oBV(kk))^2;
                end
            end
            if ~isempty(oBV)
                if tmpmin < currmin
                    currmin = tmpmin;
                    minind = [hh, ii, jj];
                end
            end
        end
    end
end
if ~isempty(oBV)
    minvalues = [T1blood(minind(1)), tblood(minind(2)), fblood(minind(3))];
    minBV = squeeze(sBV(minind(1),minind(2),minind(3),:));
else
    minvalues = [];
    minBV = [];
end
fprintf('Done\n');

% figure;
% rng(1);
% coloropt = rand(length(tblood),3);
% legendinfo = cell(length(tblood),1);
% ploth = zeros(length(tblood),1);
% for hh = 1:length(T1blood)
%     T1bloodpost = 1./(dR1 + 1./T1blood(hh));
%     for ii = 1:length(tblood)
%         for jj = 1:length(fblood)
%             dR1blood = 1./T1bloodpost - 1./T1blood(hh);
%             tmph = plot(dR1blood,squeeze(sBV(hh,ii,jj,:)),'color',coloropt(ii,:));
%             hold on;
%         end
%         ploth(ii) = tmph;
%         legendinfo{ii} = ['tau=' num2str(tblood(ii)) 'sec'];
%     end
%     % ploth(hh) = tmph;
%     %     legendinfo{hh} = num2str(T1blood(hh));
% end
% legend(ploth,legendinfo);
% hold off;

end
