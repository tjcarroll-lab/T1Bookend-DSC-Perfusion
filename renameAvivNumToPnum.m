function renameAvivNumToPnum()

D = dir('*');

for i = 1:length(D)
    
    if i < 10
        patient_number = ['P10' num2str(i)];
    elseif i<100
        patient_number = ['P1' num2str(i)];
    end
    subFold = patient_number;
    movefile(D(i).name, subFold);
    
    outputMatrix(i,1) = D(i);
    outputMatrix(i,2) = patient_number;
    
end

xlswrite('patientAvivNumToPnum.xlsx', outputMatrix);

end