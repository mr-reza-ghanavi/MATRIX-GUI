% read the text file
function readme = data_read() 
        fileID = fopen('Memory/trial_data.txt','r');
        readme = textscan(fileID,'%s','delimiter','\n');
        readme = [readme{:}];
 end
