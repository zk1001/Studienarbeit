function hrir_choice = hrir_adapt(angle)

% hrir_choice = HRIR_ADAPT(angle) selects the wanted impuls function
% for the simplicity of the main function
    hrirm90 = loadHRIR('Anechoic',80,0,-90,'front');
    hrirm45 = loadHRIR('Anechoic',80,0,-45,'front');
    hrir0 = loadHRIR('Anechoic',80,0,0,'front');
    hrir45 = loadHRIR('Anechoic',80,0,45,'front');
    hrir90 = loadHRIR('Anechoic',80,0,90,'front');
    
    if angle == -90
        hrir_choice = hrirm90;
    elseif angle == -45
        hrir_choice = hrirm45;
    elseif angle == 0
        hrir_choice = hrir0;
    elseif angle == 45
        hrir_choice = hrir45;  
    elseif angle == 90
        hrir_choice = hrir90;  
    end
    
end
