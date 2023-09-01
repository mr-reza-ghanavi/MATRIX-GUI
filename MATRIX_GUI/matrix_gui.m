%%%%%                             %%%%%
%%%%%%% Authoring and Copyright %%%%%%%
%%%%%                             %%%%%

% MAtrix_GUI - MATRIX speech intelligibility test in noise MATLAB GUI.
% Author: Reza Ghanavi
% Computing and Audio Research Laboratory - The University of Sydney
% Email: reza.ghanavi@sydney.edu.au
% Website: https://github.com/mr-reza-ghanavi
% Creation: 01-Sep-2023

% Copyright (c) 2023, Reza Ghanavi
% This software is licensed under the GNU General Public License V3.0.
% The software does not come with ANY WARRANTY.

%%%%%                      %%%%%
%%%%%%% Acknowledgement  %%%%%%%
%%%%%                      %%%%%

% The author would like to thank Wei Dai and Yi Shen [1] for providing
% the source code of the updated-maximum-likelihood under the GNU General 
% Public License V3.0 and Dave Fabry for providing Matrix-Speech-Material [2]. 

% References:

% [1] Shen, Y., Dai, W., & Richards, V. M. (2015). A MATLAB toolbox for the efficient
% estimation of the psychometric function using the updated maximum-likelihood adaptive
% procedure. Behavior research methods, 47(1), 13-26.

% [2] Speech material reference : Kelly, H., Lin, G., Sankaran, N., Xia, J., Kalluri, S.,
% & Carlile, S. (2017). Development and evaluation of a mixed gender, multi-talker matrix
% sentence test in Australian English. International journal of audiology, 56(2), 85-91.

%%%%%                                %%%%%
%%%%%%% Application of Matrix Test %%%%%%%
%%%%%                                %%%%%

% The Matrix Test is a valuable tool for comparing different scenarios for the same
% patient/subject, such as aided versus unaided hearing, pre-operative versus 
% post-operative conditions, various hearing devices, or different settings for the
% same hearing device.

% This software uses adaptive speech levels in noise. It means the noise level is
% consistently maintained at an audible level for the subject, with a recommended
% threshold of 65 dB. Initially, the first sentence is presented with a signal-to
% -noise ratio (SNR) of 0 dB.

% This software automatically adjusts the speech level based on the subjects' previous
% responses. If the subject accurately repeats three to five presented words, the speech
% level for the next presentation is reduced. Conversely, if the patient correctly 
% repeats fewer than three words, the speech level for the subsequent presentation is
% increased. 

% Typically, the subject listens to sentences delivered from the frontal loudspeaker,
% accompanied by noise specific to the test. In this software, this noise is only
% introduced during sentence presentation. 

% A 20-item test list typically lasts around 4 minutes. The precision of threshold
% estimation for the 20-trial test lists is usually within a margin of 1 dB. Minor
% variations in Speech Reception Thresholds (SRT) can result in significant differences
% in speech intelligibility due to the steepness of the logistic psychometric function.

%%%%%           %%%%%
%%%%%%% Usage %%%%%%%
%%%%%           %%%%%

% Note: please do not use this software on Cloud. Copy the MATRIX_GUI folder
% on your computer hard drive and open it in Matlab. 

% Before running: right click!! "Add to Path/Selected Folders and Subfolders".

% Syntax:  matrix_gui()

% For best performance, use a multi-channel ASIO audio interface.
% Note: This software only supports audio interfaces up to 16 channels. 

% Before starting, type MATLAB's recognizable audio interface name in the
% "Audio Device" edit field, enter the appropriate number of channels in the
% "Channel No" edit field, and click on an empty area for restarting the application. 

% Use the "Cal Sig Play/Stop" toggle button to play/stop the multi-channel 
% calibration signal via your audio interface. Note: Channel "S-1" always plays
% MATRIX target sentence signal. Other channels play different tracks of speech bubbles.

% Note: You may use your audio material with similar formats to those in the "Target"
% and "Calibration" folders.

% You can enter an appropriate value in the "Calibration [dB SPL]" edit fields
% to adjust the sound pressure level for each channel. You can use a calibrated
% SPL meter to measure the average SPL (use slow averaging) close to the pinna
% of a listener. For standard audiometric research tests, the SPL level for the
% target only and the total noise tracks must be adjusted to about 65 dB SPL at
% the listener's pinna position. Before starting the Matrix test, the signal-to-noise
% ratio (SNR) must be 0 dB SPL.

% Enter your desired trial numbers in the "Trial Numbers" edit field and click on
% an empty area to restart the system. In most audiometric applications, 25 trials
% can provide a good convergence for SRT and slope estimation after training a subject.
% The average of the last ten trials should be considered reliable audiometric results.  

% Add a unique code for your "Subject ID" and a "Test ID" for each test session in
% the main window (no space between characters).

% Use the "Start"  followed by the "Next" button to play randomly generated Matrix
% Base Sentences and record the correct words by clicking them.

% In each test session, the trial data are saved in a temporary file
% ("Results/temp/outdata.mat"). You can use this file in the case of
% any crash in your system.

% After finishing all trials, click "Save data & Close application” in the prompt window. 
%Please note that a review graph prompts and the data is saved in
% ./Results/SubjectID/ SubjectID-TestID.xlsx.

%%%%%                                   %%%%%
%%%%%%% Standard Matrix English Words %%%%%%%
%%%%%                                   %%%%% 

% Matrix English Words (A: Name, B: Verb, C: Quantity, D: Adjective, E: Object)

% A = ["Peter","Kathy","Lucy","Alan","Rachel","Barry","Steven","Thomas","Hannah","Nina"];
% B = ["got","sees","bought","gives","sold","likes","has","kept","wins","wants"];
% C = ["three","nine","five","eight","four","six","two","ten","twelve","some"];
% D = ["large","small","old","dark","thin","green","cheap","pink","red","big"];
% E = ["desks","chairs","shoes","toys","spoons","mugs","ships","rings","tins","beds"];

%%%%%                 %%%%%
%%%%%%% Main Script %%%%%%%
%%%%%                 %%%%%

%% Graphical User Interface (GUI) Function

function matrix_gui()

global app par ml;  
     
par.lock = 0;
par.leak = 0;
par.Fs = 48e3; % sampling frequency

%% User Interface

app = uifigure;
app.Position = [30,150, 1200,430]; % main screen
app.Name = "Matrix Psychometric Test GUI. Author: reza.ghanavi@sydney.edu.au © 2023";
app.Icon = "Icons/hearing.jpg";
app.WindowStyle = 'alwaysontop';
app.Visible = "on";
app.AutoResizeChildren = 'off';
app.Resize = 'off';
app.Color = [0.9,0.9,0.9];
set(figure,'Visible','off')

%% Trial Data Management
   
par.readme = data_read();
write_data(par.readme);

%% Audio Device Name

             par.text=string(par.readme(7));

             ef=uieditfield(app,'Value',par.text);
             ef.Position = [160 370 260 25];
             ef.FontName = "Helvetica";
             ef.FontSize = 16;
             ef.Placeholder = 'Text';
             ef. ValueChangedFcn = ...
                  @(textarea,event) textEnter(textarea);

     function textEnter(textarea)

         name = textarea.Value;
         par.text = cellstr(name);
         par.readme(7)=par.text;
         write_data(par.readme);
         delete(app);
         matrix_gui();

     end
 
device = data_read();
par.devicename = string(device(7));
 
%% Date Function

uilabel(app,'Position',[30 335 130 25],...
    'Text','Trial Date:',FontSize=20,FontName="Helvetica");

format short; date = clock; year = string(date(1));
month = string(date(2)); day = string(date(3));
date = [year,"/",month,"/",day];
par.readme = data_read; date_save = cellstr(join(date));
par.readme(1) = date_save; write_data(par.readme);

uieditfield(app,...
    'Position',[160 300-(1-2)*35 150 25],'Value',[date{:}],FontSize=16);

%% Participant/Trial Information

    field_name = {'First Name:', 'Last Name:',...
        'Sex/Age:','Subject ID:', 'Test ID:'};

    for i=2:6

              uilabel(app,'Position',[30 300-(i-2)*35 130 25],'Text',...
                  string(field_name(i-1)),FontSize=20,FontName="Helvetica");

              par.text=string(par.readme(i));

              editFields(i)=uieditfield(app,'Value',par.text);
              editFields(i).Position = [160 300-(i-2)*35 150 25];
              editFields(i).FontName = "Helvetica";
              editFields(i).FontSize = 16;
              editFields(i).Placeholder = 'Text';
              editFields(i). ValueChangedFcn = ...
                  @(textarea1,event) textEntered1(textarea1,i);

    end

     function textEntered1(textarea1,i)
     
         name = textarea1.Value;
         par.text = cellstr(name);
         par.readme(i)=par.text;
         write_data(par.readme);
     
     end

%% Audio Device Channel Management

readme_cal= data_read_cal();

par.text_cal=string(readme_cal(18));
ch=uieditfield(app,'Value',par.text_cal);
ch.Position = [460 370 50 25];
ch.FontName = "Helvetica";
ch.FontSize = 16;
ch.Placeholder = 'Num';
ch.ValueChangedFcn = @(textch,event) textCH(textch);

function textCH(textch)
num = textch.Value;

if str2double(num) >16

error('Channel numbers more than 16 are not supported.');

end

par.text = cellstr(num);
readme_cal(18)=par.text;
write_data_cal(readme_cal);
delete(app);
matrix_gui();
end

readme_cal= data_read_cal();
par.chnum = str2double(cell2mat(readme_cal(18)));

%% Source Management

readme_cal= data_read_cal();

    for j=1:par.chnum

    
        if j<=8

                  uilabel(app,'Position',[320 300-(j-1)*35 50 25],...
                      'Text',['S-',num2str(j)],FontSize=20,FontName="Helvetica");
                  par.text_cal=string(readme_cal(j));
                  editFields(j)=uieditfield(app,'Value',par.text_cal);
                  editFields(j).Position = [370 300-(j-1)*35 50 25];
                  editFields(j).FontName = "Helvetica";
                  editFields(j).FontSize = 16;
                  editFields(j).Placeholder = 'Num';
                  editFields(j).ValueChangedFcn = ...
                      @(textarea2,event) textEntered2(textarea2,j);

        else
        
            uilabel(app,'Position',[520 300-(j-9)*35 50 25],...
                'Text',['S-',num2str(j)],FontSize=20,FontName="Helvetica");
                  par.text_cal=string(readme_cal(j));
                  editFields(j)=uieditfield(app,'Value',par.text_cal);
                  editFields(j).Position = [570 300-(j-9)*35 50 25];
                  editFields(j).FontName = "Helvetica";
                  editFields(j).FontSize = 16;
                  editFields(j).Placeholder = 'Num';
                  editFields(j).ValueChangedFcn =...
                      @(textarea2,event) textEntered2(textarea2,j);

        end
    
    end

     function textEntered2(textarea2,j)
     
         num = textarea2.Value;
         par.text = cellstr(num);
         readme_cal(j)=par.text;
         write_data_cal(readme_cal);
     
     end

 %% Mute Control 

par.mute_coe = ones(par.chnum,1);

for k=1:par.chnum

    if k<=8

        mut = uicheckbox(app); 
        mut.Position = [423,305-(k-1)*35,70,14]; 
        mut.ValueChangedFcn = @(muty,event) mut_Callback(muty,k);
        mut.Text = "Mute";
        mut.Enable = 'on';
        mut.FontName = "Helvetica";
        mut.FontSize = 20;
        mut.FontColor = [0,0,0];

    else

        mut = uicheckbox(app); 
        mut.Position = [623,305-(k-9)*35,70,14];
        mut.ValueChangedFcn = @(muty,event) mut_Callback(muty,k);
        mut.Text = "Mute";
        mut.Enable = 'on';
        mut.FontName = "Helvetica";
        mut.FontSize = 20;
        mut.FontColor = [0,0,0];

    end
end

function mut_Callback(muty,k)

if muty.Value == 1

    par.mute_coe(k) = 0;

else
    
    par.mute_coe(k) = 1;

end
end

%% Set Trial Numbers 

par.read = data_read_cal();
trial_num=cell2mat(par.read(17));
par.N = str2double(trial_num);

uilabel(app,"Text",'Trial Numbers:',...
    'Position',[30 300-(7-2)*35 130 25],...
    FontSize=20,FontName="Helvetica");

uieditfield(app,'Value',trial_num,...
    'Position',[160 300-(7-2)*35 150 25],...
    'ValueChangedFcn',@(trials,event) trial_fnc(trials),...
    FontSize=16,FontName="Helvetica");

    function trial_fnc(trials)

        trial_num = trials.Value;
        par.N = trial_num;
        par.read(17)={trial_num};
        write_data_cal(par.read);

        delete(app);

        matrix_gui;

    end

%% Start Button

start = uibutton(app); 
start.Position = [700,50,200,200];
start.ButtonPushedFcn = @start_Callback;
start.Tag = "Try1";
start.Text = "Start";
start.FontName = "Helvetica";
start.FontSize = 45;
start.FontColor = [0,0,0];
start.BackgroundColor = [0.8,0.8,0.8];


%% Calibration Button

cal = uibutton(app,"state");
cal.Position = [33 300-(9-2)*35 270 50]; 
cal.ValueChangedFcn = @cal_Callback;
cal.Tag = "cal";
cal.Text = "Cal Sig Play/Stop";
cal.FontWeight = "bold";
cal.FontName = "Yu Gothic";
cal.FontSize = 20;
cal.FontColor = [0,0,0];

function cal_Callback(~,~)

    par.status = cal.Value;

    if par.status == 1

        cal_fnc() % apply callibration coefficients to all channels

        switch_fnc() % play/stop multi-channel calibration signal

    end
end

%% Calibration Function Constructor

    function cal_fnc()
        
        fileID = fopen('Memory/cal_data_16ch.txt','r');
        calin = fscanf(fileID,'%f',[18 1]);
        [X,par.Fs] = audioread("Calibration/s1.wav");
        L = length(X);
        calsig =zeros(L,par.chnum);

        for n=1:par.chnum 

            sig = audioread(['Calibration/s',mat2str(n),'.wav']);
            
            if length(sig)>=L
            sig = sig(1:L);
            else
            sig = [sig;zeros(L-length(sig))];
            end
            calsig(:,n) = 10^(calin(n)/20)*par.mute_coe(n) .* sig; 
        end

    audiowrite("Calibration/main_calsig.wav",calsig,par.Fs);

    end 

drawnow;

%% Calibration Signal Playback Function

    function switch_fnc()

        if ~par.lock
       
    fileReader = dsp.AudioFileReader('Calibration/main_calsig.wav');
    fileInfo = audioinfo('Calibration/main_calsig.wav');
    deviceWriter = audioDeviceWriter('Driver','ASIO',...
        'Device',par.devicename,'SampleRate',fileInfo.SampleRate);
        
        drawnow;

           while ~isDone(fileReader) && par.status==1
               
               audioData = fileReader();
               deviceWriter(audioData);
               par.leak = 1;
              
               if par.status == 0
                   break
               end

               drawnow;
           end

     par.leak=0;

     release(fileReader);
     release(deviceWriter);

        end

     end

%% Labels

uilabel(app,"Text",...
    ['MATRIX Multi-Channel Speech Intelligibility Test in ' ...
    'Noise (Reza Ghanavi © 2023)'],...
    'Position',[500 370 650 25],FontName="Helvetica");

uilabel(app,"Text",'Countdown','Position',[700,335,200,25]);

uilabel(app,...
    'Position',[950,335,200,25],...
    'Text','Base Sentence',FontSize=20,FontName="Helvetica");

uilabel(app,'Position',[30 370 130 25],'Text','Audio Device:'...
    ,FontSize=20,FontName="Helvetica");

uilabel(app,'Position',[320 335 362 25],...
    'Text','Calibration  [dB SPL] ',FontSize=20,FontName="Helvetica");

uilabel(app,'Position',[360 370 100 25],...
    'Text','Channel No: ',FontSize=20,FontName="Helvetica");

 %% label Setup

Label_Set = findall(app,'Type','uilabel');
Red = 0.2; Green = 0.2; Blue = 0.2;
set(Label_Set,'BackgroundColor',[Red Green Blue]);
set(Label_Set,'FontColor','white');
set(Label_Set,'FontSize',14);
set(Label_Set,'HorizontalAlignment','center');

%% Trial Start Function

    function start_Callback(~,~)
%% Trial Counter Parameters 

if ~par.leak

par.counter =1; % counter initialisation
par.j = 1; par.k = 2;
par.out_data = cell(7,2*par.N);
par.expo = zeros(2,2*par.N);
par.NL = 0; % noise level initialisation

%% Logistic Function Parameters

par.ndown = 2;        % up-down sweetpoint selection rule
par.method = 'mean';  % posterior distribution estimation ('mean' and 'mode')                           
par.x0 = 0;           % SNR initialisation
par.x_lim = [-20 20]; % SNR limits

par.alpha = struct(...
    'limits',[-20 20],...     % range. limits=[-20 20]
    'N',100,...               % number of alpha values. N=[100] 
    'scale','lin',...         % Choose between 'lin' and 'log' spacing. select:'lin'
    'dist','norm',...         % prior distribution. Choose between 'norm' and 'flat'. select:'norm'
    'mu',0,...                % mean of the prior distribution.mu=[0]
    'std',5 ...               % standard deviation of the prior distribution. std=[5]  
    );

par.beta = struct(...
    'limits',[0.1 4],...      % range. limits=[0.1 4]
    'N',100,...               % number of beta values. N=[100]
    'scale','lin',...         % Choose between 'lin' and 'log' spacing. select:'lin'
    'dist','norm',...         % prior distribution. Choose between 'norm' and 'flat'. select:'norm'
    'mu',0,...                % mean of the prior distribution.mu=[0]
    'std',1 ...               % standard deviation of the prior distribution. std=[1]  
    );

par.gamma = 0;                % set to [0] for logistic psychometric function estimation

par.lambda = struct(...
    'limits',[0 0.2],...      % range. limits=[0 0.2]
    'N',1,...                 % number of lambda values. set to [1] for not 'lambda' estimation
    'scale','lin',...         % Choose between 'lin' and 'log' spacing. select:'lin'
    'dist','flat',...         % prior distribution. Choose between 'norm' and 'flat'. select:'flat'
    'mu',0,...                % mean of the prior distribution. mu=[0]
    'std',0.1 ...             % standard deviation of the prior distribution. std=[0.1] 
    );

%% Maximum Likelihood Logistic Algorithm

            ml.par = par;
            ml.gamma = par.gamma;
            ml.x = [];
            ml.xnext = par.x0;
            ml.r = [];
            ml.phi = [];
            ml.n = 0;
            ml.psycfun = @(x,a,b,g,l) g+(1-g-l).*(1+exp(-(x-a).*b)).^(-1);
            ml.swpts_idx = [];

            if ml.par.alpha.N > 1

                ml.swpts_idx(end+1) = 2;
            end

            if ml.par.beta.N > 1

                ml.swpts_idx(end+1:end+2) = [1 3];

            end

            if ml.par.lambda.N > 1

                ml.swpts_idx(end+1) = 4;

            end

            ml.swpts_idx = sort(ml.swpts_idx);
            ml.nsteps = length(ml.swpts_idx);
            ml.step_flag = 0;
            ml.track_direction = 0;
            ml.rev_flag = [];
            ml.nrev = 0;

   %% First Sweet Point

            if  ml.par.x0 == ml.par.x_lim(1)

                ml.current_step = 1;

            elseif ml.par.x0 == ml.par.x_lim(2)

                ml.current_step = ml.nsteps;

            else

                ml.current_step = ceil(ml.nsteps/2);
            end  
       
    %% Set the Estimation Space 

            ml.alpha = setParSpace(ml.par.alpha);
            ml.beta = setParSpace(ml.par.beta);
            ml.lambda = setParSpace(ml.par.lambda);
            
            [ml.a,ml.b,ml.l] = meshgrid(ml.alpha,ml.beta,ml.lambda);
            
            A = setPrior(ml.a, ml.par.alpha);
            B = setPrior(ml.b, ml.par.beta);
            L = setPrior(ml.l, ml.par.lambda);

            ml.p = log(prepare_prob(A.*B.*L)); 

%% Go to 'matrix function'

matrix();

end

    end
end

%% Main Matrix Function

function matrix()

global    par 

clf 

%% Matrix English Words (A: Name, B: Verb, C: Quantity, D: Adjective, E: Object)

A = ["Peter","Kathy","Lucy","Alan","Rachel","Barry","Steven","Thomas","Hannah","Nina"];
B = ["got","sees","bought","gives","sold","likes","has","kept","wins","wants"];
C = ["three","nine","five","eight","four","six","two","ten","twelve","some"];
D = ["large","small","old","dark","thin","green","cheap","pink","red","big"];
E = ["desks","chairs","shoes","toys","spoons","mugs","ships","rings","tins","beds"];

%% Matrix Audio File

voice1 = {'a1' 'a2' 'a3' 'a4' 'a5' 'a6' 'a7' 'a8' 'a9' 'a10'};
voice2 = {'b1' 'b2' 'b3' 'b4' 'b5' 'b6' 'b7' 'b8' 'b9' 'b10'};
voice3 = {'c1' 'c2' 'c3' 'c4' 'c5' 'c6' 'c7' 'c8' 'c9' 'c10'};
voice4 = {'d1' 'd2' 'd3' 'd4' 'd5' 'd6' 'd7' 'd8' 'd9' 'd10'};
voice5 = {'e1' 'e2' 'e3' 'e4' 'e5' 'e6' 'e7' 'e8' 'e9' 'e10'};

%% Matrix Signal Constructor 


    for ind = 1:10 % names 
        s1.(voice1{ind}) = audioread(['Target/',voice1{ind},'.wav']);
    end
 
    for ind = 1:10 % verbs
        s2.(voice2{ind}) = audioread(['Target/',voice2{ind},'.wav']);
    end

    for ind = 1:10 % numerals
        s3.(voice3{ind}) = audioread(['Target/',voice3{ind},'.wav']);
    end

    for ind = 1:10 % adjectives
        s4.(voice4{ind}) = audioread(['Target/',voice4{ind},'.wav']);
    end

    for ind = 1:10 % objects
        s5.(voice5{ind}) = audioread(['Target/',voice5{ind},'.wav']);
    end

%% matrix random Target Speech Constructor

randset = randi([1 10],1,5,'single'); 

par.speech = [A(randset(1)),B(randset(2)),...
    C(randset(3)),D(randset(4)),E(randset(5))];

target = [s1.(voice1{randset(1)});s2.(voice2{randset(2)});...
    s3.(voice3{randset(3)});s4.(voice4{randset(4)});s5.(voice5{randset(5)})];

target = [zeros(par.Fs,1);target;zeros(par.Fs/2,1)]; 
L_tar = length(target);

audiowrite("Audio/aud1.wav",target,par.Fs);

sw_L = zeros(par.chnum,1);
sw_L(1) = par.NL;
  
%% Target and Bubble Noise Signals

        fileID = fopen('Memory/cal_data_16ch.txt','r');
        cal_read = fscanf(fileID,'%f',[17 1]);
        main_matx_sig = zeros(L_tar,par.chnum);

    for idx = 1:par.chnum

     [matx_sig,par.Fs] = audioread(['Audio/aud',num2str(idx),'.wav']);

     matx_sig = 10^(sw_L(idx)/20)*10^(cal_read(idx)/20)*par.mute_coe(idx) .* matx_sig;

     main_matx_sig(:,idx) = matx_sig(1:L_tar);

    end

audiowrite("Audio/main_sig.wav",main_matx_sig,par.Fs);


%% DSP Setup

    fileReader = dsp.AudioFileReader('Audio/main_sig.wav');

    fileInfo = audioinfo('Audio/main_sig.wav');

    deviceWriter = audioDeviceWriter('Driver','ASIO','Device',...
        par.devicename,'SampleRate',fileInfo.SampleRate);

%% Text Generator
 
par.val = zeros(5,1);
par.hand=[];
par.valt = 0;

global app;

        for nam = 1:5
            
            box= uicheckbox(app);
            box.Position = [960,315-(nam*60),250,80];
            box.ValueChangedFcn = @(source,event) r_Callback(source,nam);
            box.Text = par.speech(nam);
            box.FontName = "Yu Gothic";
            box.FontSize = 45;
            box.FontColor = [0,0,0];
            par.hand(nam)=box;
  
        end
        
        function r_Callback(source,nam)
           
            par.val(nam) = source.Value;

            if par.val(nam)==1

                source.FontColor = [1,0,0];
                source.FontWeight = 'bold';

            else

                source.FontColor = [0,0,0];
                source.FontWeight = 'normal';

            end
        end
        
%% Next Button

next = uibutton(app); 
next.Position = [700,50,200,200]; 
next.ButtonPushedFcn = @next_Callback;
next.Tag = "next";
next.Text = "Next";
next.FontName = "Helvetica";
next.FontSize = 45;
next.FontColor = [0,0,0];
next.BackgroundColor=[0.8,0.8,0.8];

%% Matrix Multichannel Playback Loop

drawnow;

while ~isDone(fileReader)

    audioData = fileReader();
    deviceWriter(audioData);

    par.lock = 1;

    if isDone(fileReader)

        break;
        
    end

    drawnow;
end

par.lock = 0;

release(fileReader);
release(deviceWriter);

end

%% Matrix Call-Back Function  

function next_Callback(~, ~)

global  par app ml;

if ~par.lock && ~par.leak

delete(par.hand(:));

par.valt = sum(par.val);

par.counter = par.counter +1;

%% Data Sampling

if par.counter <= par.N+1

    for i=1:5

        par.out_data(i,par.j)= cellstr(par.speech(i));
        par.out_data(i,par.k)={par.val(i)*1};

    end
    
    par.out_data(6,par.k)= {round(par.NL,2)};
  
    par.j=par.j+2;
    par.k=par.k+2;

%% Hypothesis Test  

    if     par.valt == 5
        r = 1; 
    elseif par.valt == 4
        r = 0.8;
    elseif par.valt == 3
        r = 0.6;
    elseif par.valt == 2
        r = 0.4;
    elseif par.valt == 1
        r = 0.2;
    elseif par.valt == 0
        r = 0;
    end

%% Maximum Likelihood Iteration Algorithm

            ml.n = ml.n + 1;
            ml.x(end+1,:) = ml.xnext;
            ml.r(end+1,:) = r;

            ml.p = ml.p +...
                log(prepare_prob(ml.psycfun(ml.xnext,ml.a, ml.b, ml.gamma, ml.l)).^r) + ...
                log(prepare_prob(1-ml.psycfun(ml.xnext,ml.a, ml.b, ml.gamma, ml.l)).^(1-r));

            ml.p = ml.p-max(max(max(ml.p)));
            
            if strcmp(ml.par.method, 'mode')

                [~,idx] = max(reshape(ml.p,numel(ml.p),1));
                ml.phi(end+1,:) = [ml.a(idx), ml.b(idx), ml.gamma, ml.l(idx)];

            elseif strcmp(ml.par.method, 'mean')

                pdf_tmp = exp(ml.p);
                pdf_tmp = pdf_tmp/sum(sum(sum(pdf_tmp)));
                % alpha
                alpha_est_tmp = sum(sum(sum(pdf_tmp.*ml.a)));
                % beta
                beta_est_tmp = sum(sum(sum(pdf_tmp.*ml.b)));
                % lambda
                lambda_est_tmp = sum(sum(sum(pdf_tmp.*ml.l)));
                % combine together
                ml.phi(end+1,:) = [alpha_est_tmp, beta_est_tmp, ml.gamma, lambda_est_tmp];

            else

                error('Choose between "mean" and "mode" for "method".');

            end
          
%% Next Signal Strength Estimation (Sweet Points)

             swpt = [logit_sweetpoints(ml.phi(end,:)) ml.par.x_lim(2)];  
             swpt = max(min(swpt(ml.swpts_idx),ml.par.x_lim(2)),ml.par.x_lim(1)); 

%% N-Down, 1-Up Algorithm

            ml.rev_flag(end+1,:) = 0;

            if r>=0.5   

                if ml.step_flag == ml.par.ndown-1 % reached ndown

                    ml.current_step = max(ml.current_step-1,1);
                    newx = swpt(ml.current_step);
                    ml.step_flag = 0;

                    if ml.track_direction == 1

                        ml.rev_flag(end,:) = 1;
                        ml.nrev = ml.nrev + 1;

                    end

                    ml.track_direction = -1;

                elseif ml.step_flag <ml.par.ndown-1 % has not reached ndown

                    newx = swpt(ml.current_step);
                    ml.step_flag = ml.step_flag+1;

                end

            elseif r<0.5    

                ml.current_step = min(ml.current_step+1,ml.nsteps);
                newx = swpt(ml.current_step);
                ml.step_flag = 0;

                if ml.track_direction == -1

                    ml.rev_flag(end,:) = 1;
                    ml.nrev = ml.nrev + 1;

                end

                ml.track_direction = 1;
            end
            
            ml.xnext = newx;
            
            
            ml.PC = sum(ml.r)/length(ml.r); % update percent correct
    
            par.NL = ml.xnext; % update the target signal level
    
%% Counter

par.countdwn = par.N - par.counter+1; 
par.countdwn = num2str(par.countdwn);

uieditfield(app,'Value',par.countdwn,...
'Position',[760 260 80 60],FontSize=60,FontName="Helvetica",...
BackgroundColor=[0.9,0.9,0.9],FontColor=[1,0.2,0],...
FontWeight='bold',Visible='on',Editable='on');


%% Save a Temporary File

 coeff = [(ml.phi(:,1))' ; (ml.phi(:,2))']; % extract SRT and slope

 par.expo(:,par.counter*2)=coeff(:,end);
 par.expo = par.expo(:,1:2*par.N);
 par.outdata = [par.out_data;num2cell(par.expo)];
 out = par.outdata;

 if ~exist('Results/temp', 'dir')
       mkdir('Results/temp');
 end

 save("Results/temp/outdata.mat","out");

%% Return Matrix

matrix();

else
%% Export Data to .xlsx File

selection = uiconfirm(app,'End of Trial!','Declaration','Options',...
    {'Save data & Close application'},'Icon','Icons/success.jpg');

        switch selection
            case 'Save data & Close application'
                
                data = data_read();
                user = data(5);
                test = data(6);

                mkdir('Results/',string(user));

                fname = ['Results','/',string(user),'/',string(user),'-',string(test),'.xlsx'];
                fname = strjoin(fname);
                fname = convertStringsToChars(fname);
                fname = fname(~isspace(fname));

                trial_field = fopen("Memory/trial_data.txt","r");
                trial = textscan(trial_field,'%s','Delimiter','\r');
                trial = [trial{:}];

                fileName_template = 'Template/results_template.xlsx';
                copyfile(fileName_template,fname);

                writecell(trial,fname,'Sheet',1,'Range','B3', 'AutoFitWidth', false);
                writecell(par.outdata,fname,'Sheet',1,'Range','D2','AutoFitWidth', false);

%% Show Results

                figure("Position",[200,300, 600,500],'Name','Results');

                subplot(2,1,1)
                plot(ml.phi(:,1),'k','LineWidth',2);hold;
                line([0 par.N],[0 0],'Color','k','LineStyle','--',...
                    'LineWidth',1);
                T = ['Matrix Sentence Test Results [',user,'-',test,']'];
                T = strjoin(T);
                T = convertStringsToChars(T);
                T = T(~isspace(T));
                title(T,FontSize=12);
                ylabel('SRT (50%) [dB]');
                set(gca,'XTickLabel',[]);
                axis([0 par.N ml.par.alpha.limits+[-3 -7]]);
                grid("on");
                grid("minor");
                
                subplot(2,1,2)
                plot((ml.phi(:,2)./4)*100,'k','LineWidth',2);hold;
                line([0 par.N],[20 20],'Color','k','LineStyle','--',...
                    'LineWidth',1);
                ylabel('Slope [% per-dB]');
                grid("on");
                grid("minor");
                axis([0 par.N [0 40]]);
                set(gca,'XTickLabel',[],'YTick',[0 20 40]);
                xlabel('Trials');


%% Clear Application

delete(app);

        end

end
    
end
end
