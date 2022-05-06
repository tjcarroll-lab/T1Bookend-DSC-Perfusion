function [ibat bat brt] = findBAT(Signal)

%smoothering
Signal = ([Signal(1) Signal(1:length(Signal)-1)] + 4*Signal + ...
          [Signal(2:length(Signal)) Signal(length(Signal))])./6;
      
% predifinition
ibat=0; bat=0; brt=0;

% find bolus arrival time point
[smin nmin] = min(Signal);
[smax nmax] = max(Signal(2:nmin));

count = nmin;
index = 1;

tempmean = mean(Signal(5:10));
tempstd = std(Signal(5:10));

% find Bolus arrival time point
if count > 5
    while index
        mean2min1 = mean(Signal(1:count));  
        std2min1 = std(Signal(1:count));
        nstd2min1 = std2min1/mean2min1;

        mean2min2 = mean(Signal(1:count-1));
        std2min2 = std(Signal(1:count-1));
        nstd2min2 = std2min2/mean2min2;
        
        if (0.9*(tempmean - smin) < (Signal(count) - smin)) &...
            (nstd2min2/nstd2min1 > 0.95)
            index=0;
        end
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
            if ((diff_sig(count)<0) & (0.6*(tempmean-smin)<(Signal(count)-smin)))...
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
