function processBar(total, num, hwait)

if total - num <= 1
    waitbar(num/total, hwait, 'Test to be accompleted');
    pause(0.05);
else
    PerStr = round(num/total*100);
    str = ['Running ', num2str(PerStr), '%'];
    waitbar(num/total, hwait, str);
    pause(0.05);
end

end