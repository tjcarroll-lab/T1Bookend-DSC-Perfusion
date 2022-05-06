function [] = correctionForT1maps(T1pre, T1post)

TOLERANCE = 1;

for xLocation = 1:size(T1pre, 1)
    counter = 1;
    for yLocation = 1:size(T1pre, 2)
        if (T1pre(xLocation,yLocation) > TOLERANCE)
            counter = counter + 1;
            T1contour(counter) = T1pre(xLocation,yLocation);
            plot(T1contour); pause;
        end
    end
end

end