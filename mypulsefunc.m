function y = mypulsefunc(t)
y = nan(1, size(t, 2));
for i = 1:length(t)
    if t(i) >= 0 && t(i) <= 40
        y(i) = 1;
    else
        y(i) = 0;
    end
end
end