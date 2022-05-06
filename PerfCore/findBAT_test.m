function [ibat bat brt] = findBAT_test(Signal)

%smoothering
%Signal = ([Signal(1) Signal(1:length(Signal)-1)] + 4*Signal + ...
%          [Signal(2:length(Signal)) Signal(length(Signal))])./6;
      
% predifinition
ibat=0; bat=0; brt=0;

% 20210408 YIJ: Remove signal spikes for better local minimum calculation
%tempbase = mean(Signal(5:10));
tempbasemax = max(Signal(5:10));
tempsig = Signal;
tempsig(Signal > tempbasemax) = tempbasemax;

% find bolus arrival time point
[tf,p] = islocalmin(tempsig);
%[~,nmin] = max(p);
% 20210426 YIJ: if there are multiple big drops, choose earliest one
p(p < 0.75*max(p)) = 0;
nmin = find(p,1);
if isempty(nmin)
    [~,nmin] = max(p);
end
smin = Signal(nmin);
%[smin nmin] = min(Signal);
%[smax nmax] = max(Signal(2:nmin));

% Get rough estimate of arrival time by looking at derivative
%tempdiff = smooth(diff(Signal),.05,'sgolay',1);
tempdiff = diff(Signal);

% tmpfit = fit([1:length(tempdiff)]',tempdiff','smoothingspline','smoothingparam',.9);
% tempdiff = reshape(tmpfit(1:length(tempdiff)),[1 length(tempdiff)]);

tempmin = islocalmin(tempdiff(1:nmin),'minseparation',.10*length(Signal));
if tempmin == 0
    tempbat = 10;
else
    tempmin = find(tempmin);
    tempmin = tempmin(end);
    tempbat = tempmin - (nmin-tempmin);
end

tempmean = mean(Signal(5:tempbat));
%------- Case 14 z 20 x 54 y 51
%disp(tempbat);
%disp(tempmean);
%tempmean = mean(Signal(5:10));
%tempstd = std(Signal(5:10));

count = nmin;
index = 1;

% find Bolus arrival time point
if count > 5
    while index
        mean2min1 = mean(Signal(5:count));  
        std2min1 = std(Signal(5:count));
        nstd2min1 = std2min1/mean2min1;

        mean2min2 = mean(Signal(5:count-1));
        std2min2 = std(Signal(5:count-1));
        nstd2min2 = std2min2/mean2min2;
        
        %fprintf('%f | %f | %f | %f\n', count, nstd2min2/nstd2min1, Signal(count)-smin, ...
        %    std(tempdiff(5:count-2))/std(tempdiff(5:count-1)));
        
%         fprintf('%d | %f (%f) | %f | %f | %f | %f (%f)\n',count,...
%             Signal(count)-smin,0.90*(tempmean - smin),...
%             nstd2min2/nstd2min1,...
%             std(tempdiff(5:count-2))/std(tempdiff(5:count-1)),...
%             std(tempdiff(count-2-2:count-2))/std(tempdiff(count-1-2:count-1)),...
%             mean(tempdiff(count-2:count-1)), 0.10*mean(tempdiff(count:count+1)));
        
        if (0.90*(tempmean - smin) < (Signal(count) - smin)) &...
            (nstd2min2/nstd2min1 > 0.95)
            index=0;
            %disp('option 1');
            %disp(0.9*(tempmean-smin));
        %elseif (0.50*(tempmean - smin) < (Signal(count) - smin)) && std(tempdiff(5:count-2))/std(tempdiff(5:count-1)) > 0.9
            %index=0;
            %disp('option 2');
        elseif (0.80*(tempmean - smin) < (Signal(count) - smin)) && std(tempdiff(count-2-2:count-2))/std(tempdiff(count-1-2:count-1)) > 0.9 && ...
                tempdiff(count-1) > mean(tempdiff(5:tempbat-1)) - 1*std(tempdiff(5:tempbat-1))
            index=0;
            %disp('option 3');
        elseif (0.75*(tempmean - smin) < (Signal(count) - smin)) && tempdiff(count-1) > 0 && tempdiff(count) < 0 && tempdiff(count+1) < 0
            index=0;
            %disp('option 4');
        elseif (0.90*(tempmean - smin) < (Signal(count) - smin)) && ...
                mean(tempdiff(count-2:count-1)) > 0.10*mean(tempdiff(count:count+1))
            index=0;
            %disp('option 5');
        end
        %ADD find local min of derivative from 1 to min of signal
        count = count - 1;
        if count < 5
            break
        end
    end
    bat = count + 1;
else
    %YIJ 20190409: if count is not greater than 5, set bat to 2 instead of
    %zero
    bat = 2;
end

% find starting point
count = 1;
index = 1;
if bat > 5
    while index
        mean2min1 = mean(Signal(count:bat));
        std2min1 = std(Signal(count:bat));
        nstd2min1 = std2min1/mean2min1;

        mean2min2 = mean(Signal(count+1:bat));
        std2min2 = std(Signal(count+1:bat));
        nstd2min2 = std2min2/mean2min2;

        if nstd2min2/nstd2min1 > 0.98
            index=0;
        end
        count = count + 1;
        if count > bat
            ibat = 2;
            break
        else
            ibat = count-1;
        end
    end
else
    %YIJ 20190405: if bat is not greater than 5, set ibat to 1 (instead of
    %zero)
    ibat = 1;
end

%find rec.
diff_sig = diff(Signal);

index = 1;
count=nmin;
if nmin < length(Signal)-1 && nmin > bat %YIJ 20190409 lowest signal point should be later than bat
    while index
        if ibat & bat
            %YIJ 20170620 changed from 0.5 to 0.6 for noisy canine data
            %(20190520) changed back to 0.5
            if ((diff_sig(count)<0) & (0.5*(tempmean-smin)<(Signal(count)-smin)))...
                  | ((Signal(count)-smin)> 0.9*(mean(Signal(ibat:bat)-smin)))  
                index=0;
            end
        end
        count = count + 1;
        if ( count > length(Signal)-1) 
            index=0;
        end
    end
    brt = count -1 ;
    
    if brt == length(Signal)-1
        brt = min(length(Signal)-1, nmin + (nmin - bat) + 2);
    end
else
    %YIJ 20180208: if finding brt fails, set it to same width as bat to
    %peak plus 2 timepoints
    brt = min(length(Signal)-1, nmin + (nmin - bat) + 2); %length(Signal);
end

%fprintf('Test value: %f\n', std(Signal(5:brt))/std(Signal(5:bat)));
