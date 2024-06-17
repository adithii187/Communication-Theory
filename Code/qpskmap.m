function symbol = qpskmap(bit1, bit2)
    if (bit1 == 0 && bit2 == 0)
        symbol = 1 + 1i*1;
    elseif (bit1 == 1 && bit2 == 0)
        symbol = -1 + 1i*1;
    elseif (bit1 == 1 && bit2 == 1)
        symbol = -1 - 1i*1; 
    elseif (bit1 == 0 && bit2 == 1)
        symbol = 1 - 1i*1;
    end
end
