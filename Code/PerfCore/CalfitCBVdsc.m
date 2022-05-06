function fitteddata = CalfitCBVdsc(XCBVss,vecCBVss,cutoffs,PN);

 
 [a xmin] = min(abs(XCBVss - cutoffs(1)));
 [b xmax] = min(abs(XCBVss - cutoffs(2)));
 
 %%%%%%grady Addition outisde this range is prob error
options = fitoptions('gauss1');
options.Lower = [-Inf 0.01 -Inf];
%options.Upper = [Inf 30 Inf];
options.Upper = [Inf Inf Inf];

fitteddata = fit(XCBVss(xmin:xmax)',vecCBVss(xmin:xmax)','gauss1', options);
figure;
plot(fitteddata, XCBVss(xmin:xmax),vecCBVss(xmin:xmax))
%title(strcat('CBVdsc',num2str(fitteddata.b1))); pause(1);
title(strcat('CBVdsc',num2str(fitteddata.b1), '....', PN)); pause(1);


