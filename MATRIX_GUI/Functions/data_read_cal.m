% read the cal file
function readme = data_read_cal() 
        fileID = fopen('Memory/cal_data_16ch.txt','r');
        readme = textscan(fileID,'%s','delimiter','\n');
        readme = [readme{:}];
 end
