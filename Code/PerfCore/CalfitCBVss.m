function fitteddata = CalfitCBVss(XCBVss,vecCBVss,cutoffs, PN);

 
 [a xmin] = min(abs(XCBVss - cutoffs(1)));
 [b xmax] = min(abs(XCBVss - cutoffs(2)));

%%%%%%grady Addition, outside this range 0 to 30 is prob error
options = fitoptions('gauss1');
options.Lower = [-Inf 0.01 -Inf];
options.Upper = [Inf 30 Inf];
 
fitteddata = fit(XCBVss(xmin:xmax)',vecCBVss(xmin:xmax)','gauss1', options);
figure;
plot(fitteddata, XCBVss(xmin:xmax),vecCBVss(xmin:xmax))
title(strcat('CBVssFast',num2str(fitteddata.b1), '....', PN)); pause(1);

