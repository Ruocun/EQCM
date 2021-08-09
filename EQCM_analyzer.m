% INTRODUCTION
% 1. Please put the '.mpt' file you want to analyze in the same folder as
% this file.
% 2. Update the mass and molecular weight
% 3. Click 'Run'
% 4. You should expect a text file that contains the analyzed outcome

mass = 1; %mass of the electrode (unit: g) *** unknown mass ***
MW = 183.84+16*3+18; % molecular weight of your sample 

%% Load files
allFiles = dir( '*.mpt' );% get the information all .mpt file in the current folder
allNames = {allFiles(~[allFiles.isdir]).name};% get the names of the .mpt files
Num_El = numel(allNames);% count the number of files in the folder
Ewe = NaN(Num_El,1);% create a matrix to contain working electrode potential
I = NaN(Num_El,1);% create a matrix to contain current
Cath_Cap = NaN(Num_El,1);% create a matrix to contain cathodic capacity
Anod_Cap = NaN(Num_El,1);% create a matrix to contain anodic capacity
sweep_rates = NaN(Num_El,1);% create a matrix to sweep rates
i=1;
%  for i = 1:1:Num_El% a for loop function to extract data from each .mpt file
     % To find out headerlines and sweep rate
    textFileName = char(allNames(i)); % convert the name of the i-th file into a string format readable by Matlab
    fileID = fopen(textFileName,'rt');% open the i-th file for getting headerline and sweep rate infos
    Full_text = textscan(fileID,'%s','Delimiter','\n');% scan the text in the i-th file and save to D
    fclose(fileID);% close the file
    Headerlines_Num = str2num(Full_text{1}{2}(19))*10+str2num(Full_text{1}{2}(20));% Find the number of headerlines in the 2nd row of D{1} which are stored in the 19th and 20th positions and convert the string format to a number format
        % find the sweep rate
    for sweep_rate_finder = 1:Headerlines_Num
        Sweep_Rate_Line = Full_text{1}{sweep_rate_finder};
        Sweep_Rate_Key = 'dE/dt               ';
        Sweep_Rate_Index = strfind(Sweep_Rate_Line, Sweep_Rate_Key);
        if Sweep_Rate_Index == 1
            sweep_rate(i) = sscanf(Sweep_Rate_Line(Sweep_Rate_Index(1) + length(Sweep_Rate_Key):end), '%g', 1);% scan the line that stores the sweep rate data from the position that the Sweep_Rate_Key ends to where the line ends. The number is stored in the i-th position in the sweep_rate matrix.
        end
    end
    % to read the columns
    fileID = fopen(textFileName,'rt');% open the i-th file for getting actual data
    OG_data = textscan(fileID,'%d %d %d %d %d %f %f %f %f %d %f %f %f','Headerlines',Headerlines_Num);% obtain data
    fclose(fileID);% close the file
    
    %% Frequency background subtraction 
    F = double(OG_data{12});% Frequency
    Cyc_No = double(OG_data{10});
    time = double(OG_data{6})/3600;
    Ewe = double(OG_data{8});
    I = double(OG_data{9});
    % find out cycle change index
    Cyc_Index = nan(max(Cyc_No)-min(Cyc_No),1);
    f_cyc_0 = nan(max(Cyc_No)-min(Cyc_No),1);
    time_cyc_0 = nan(max(Cyc_No)-min(Cyc_No),1);
    count = 1;
    for i = 2:numel(Cyc_No)
        if Cyc_No(i) > Cyc_No(i-1)
            Cyc_Index(count) = i;
            f_cyc_0(count) = F(i);
            time_cyc_0(count) = time(i);
            count = count + 1;
        end
    end
    % first point of every cycle
    time_cyc_0 = [time(1);time_cyc_0;time(end)];
    f_cyc_0 = [F(1);f_cyc_0;F(end)];
    % spline background fit
    baseline_func = spline(time_cyc_0',f_cyc_0');
    xx = min(time_cyc_0):0.1:max(time_cyc_0);
    figure(1)
    plot(time_cyc_0,f_cyc_0,'o',xx,ppval(baseline_func,xx),'-')
    hold on
    plot(time, F)
    hold off
    baseline_fit = ppval(baseline_func,time');
    delta_F_clean = F - baseline_fit';
    %% Frequency derivative
    diff_delta_F = diff([eps; delta_F_clean(:)])./diff([eps; time(:)]);
    SGF = 101;% how heavily you are filtering
    filt_diff_delta_F = sgolayfilt(diff_delta_F, 1, SGF);% Savitzky_Golay filtering
    figure(2)
    subplot(1,2,1)
    yyaxis left
    plot(Ewe,-filt_diff_delta_F)
    xlabel('Ewe')
    ylabel('negative frequency change rate')
    yyaxis right
    plot(Ewe,I)
    ylabel('current (mA)')
    subplot(1,2,2)
    yyaxis left
    plot(Ewe,filt_diff_delta_F)
    xlabel('Ewe')
    ylabel('positive frequency change rate')
    yyaxis right
    plot(Ewe,I)
    ylabel('current (mA)')
    
    
    %% Data compilation
    Data = [Ewe, ... % working electrode potential
          I, ... % current
          Cyc_No, ... % cycle number
          time, ... % time in hours
          delta_F_clean,...% frequency with no background
          filt_diff_delta_F]; % frequency change derivatives


    %% Average CV and massograms
    % sort data in cycle based format
    for i = 1:max(Cyc_No)
        Data_sort_cycle{i} = Data(any(Data(:,3)==i,2),:);
    end
    % sort data into cathodic and anodic sweeps (Ewe needs to be unique for
    % the spline to work
    Sweep_turnover = sort([strfind(diff(Data(:,1))'>=0, [0 1]) strfind(diff(Data(:,1))'>=0, [1 0])]);
    Sweep_turn_index = zeros(2*max(Cyc_No),1);
    for i = 1:max(Cyc_No)
            [max_Ewe Sweep_turn_index(2*i)] = max(Data(any(Data(:,3)==i,2),1));
            [min_Ewe Sweep_turn_index(2*i-1)] = min(Data(any(Data(:,3)==i,2),1));
        if i > 1    
            Sweep_turn_index(2*i-1) = Sweep_turn_index(2*i-1)+ Cyc_Index(i-1) - 1;
            Sweep_turn_index(2*i) = Sweep_turn_index(2*i)+ Cyc_Index(i-1) - 1;
        end
    end
    
    for i = 1:length(Sweep_turn_index)
        if i == 1
            Data_sort_sweep{i} = Data(1:Sweep_turn_index(i),:);
        else
            Data_sort_sweep{i} = Data(Sweep_turn_index(i-1)+1:Sweep_turn_index(i),:);
        end
    end
    
    
    % calculate I, F, dF at specific potentials
    Ewe_min = round(min(Ewe),1);
    Ewe_max = round(max(Ewe(Cyc_Index(1):Cyc_Index(2))),1);
    Ewe_art = Ewe_min:0.005:Ewe_max; % set Ewe range with 5 mV spacing, art = artificial? will be confusing :/
    Ewe_art_invert = Ewe_max:0.005:Ewe_min;
    I_Ewe_art_compile = nan(length(Ewe_art),length(Sweep_turn_index)/2-1);
    F_Ewe_art_compile = nan(length(Ewe_art),length(Sweep_turn_index)/2-1);
    dF_Ewe_art_compile = nan(length(Ewe_art),length(Sweep_turn_index)/2-1);
    I_Ewe_art_invert_compile = nan(length(Ewe_art),length(Sweep_turn_index)/2-1);
    F_Ewe_art_invert_compile = nan(length(Ewe_art),length(Sweep_turn_index)/2-1);
    dF_Ewe_art_invert_compile = nan(length(Ewe_art),length(Sweep_turn_index)/2-1);
    
    for i = 2:length(Sweep_turn_index)-1 % ignore the 1st and last sweeps
       Ewe_sort = Data_sort_sweep{i}(:,1);
       I_sort = Data_sort_sweep{i}(:,2);
       F_sort = Data_sort_sweep{i}(:,5);
       dF_sort = Data_sort_sweep{i}(:,6);
       % remove repeated points
       Ewe_repeat = find (diff(Ewe_sort)==0);
       for j = 1:length(Ewe_repeat)
           Ewe_sort(Ewe_repeat(j)-j+1) = [];
           I_sort(Ewe_repeat(j)-j+1) = [];
           F_sort(Ewe_repeat(j)-j+1) = [];
           dF_sort(Ewe_repeat(j)-j+1) = [];
       end
       % spline fit
       I_E_spline = spline(Ewe_sort',I_sort');
       F_E_spline = spline(Ewe_sort',F_sort');
       dF_E_spline = spline(Ewe_sort',dF_sort');
       % calculate fitted values
       if bitget(i,1)==0 %bitget gives 1 if it's odd and 0 if it's even
           I_art = ppval(I_E_spline,Ewe_art);
           F_art = ppval(F_E_spline,Ewe_art);
           dF_art = ppval(dF_E_spline,Ewe_art);
           % compilation
           I_Ewe_art_compile(:,i/2) = I_art'; % anodic cycle
           F_Ewe_art_compile(:,i/2) = F_art';% anodic cycle
           dF_Ewe_art_compile(:,i/2) = dF_art';% anodic cycle
       elseif bitget(i,1)==1
           I_art = ppval(I_E_spline,Ewe_art);
           F_art = ppval(F_E_spline,Ewe_art);
           dF_art = ppval(dF_E_spline,Ewe_art);
           % compilation
           I_Ewe_art_invert_compile(:,(i-1)/2) = I_art';% cathodic cycle
           F_Ewe_art_invert_compile(:,(i-1)/2) = F_art';% cathodic cycle
           dF_Ewe_art_invert_compile(:,(i-1)/2) = dF_art';% cathodic cycle
       end
    end
    
    % calculate averages
    I_anodic_avg = mean(I_Ewe_art_compile(:,2:8),2); % 2nd ~ 8th cycles
    F_anodic_avg = mean(F_Ewe_art_compile(:,2:8),2); % 2nd ~ 8th cycles
    dF_anodic_avg = mean(dF_Ewe_art_compile(:,2:8),2); % 2nd ~ 8th cycles
    I_cathodic_avg = mean(I_Ewe_art_invert_compile(:,1:7),2); % 2nd ~ 8th cycles
    F_cathodic_avg = mean(F_Ewe_art_invert_compile(:,1:7),2); % 2nd ~ 8th cycles
    dF_cathodic_avg = mean(dF_Ewe_art_invert_compile(:,1:7),2); % 2nd ~ 8th cycles
    % calculate standard deviations
    I_anodic_std = std(I_Ewe_art_compile(:,2:8),0,2);
    F_anodic_std = std(F_Ewe_art_compile(:,2:8),0,2);
    dF_anodic_std = std(dF_Ewe_art_compile(:,2:8),0,2);
    I_cathodic_std = std(I_Ewe_art_invert_compile(:,1:7),0,2);
    F_cathodic_std = std(F_Ewe_art_invert_compile(:,1:7),0,2);
    dF_cathodic_std = std(dF_Ewe_art_invert_compile(:,1:7),0,2);
    % put together one cycle
    E_art_all = [flip(Ewe_art');Ewe_art(2:end)'];
    I_avg = [flip(I_cathodic_avg);I_anodic_avg(2:end)];
    F_avg = [flip(F_cathodic_avg);F_anodic_avg(2:end)];
    dF_avg = [flip(dF_cathodic_avg);dF_anodic_avg(2:end)];
    I_std = [flip(I_cathodic_std);I_anodic_std(2:end)];
    F_std = [flip(F_cathodic_std);F_anodic_std(2:end)];
    dF_std = [flip(dF_cathodic_std);dF_anodic_std(2:end)];
    % compilation
    E_I_comp = [E_art_all,I_avg,I_std];
    E_F_comp = [E_art_all,F_avg,F_std];
    E_dF_comp = [E_art_all(2:end-1),dF_avg(2:end-1),dF_std(2:end-1)];%the 1st and last points being artifact
    E_I_F_dF_comp = [E_art_all,I_avg,I_std,F_avg,F_std,[nan;dF_avg(2:end-1);nan],[nan;dF_std(2:end-1);nan]];
    % text output
    Name = char(allNames(1));
    EIFdF_Name = strrep(Name, Name(end-4:end), '_E_I_F_dF.txt');% create new file names
    fileID = fopen(char(EIFdF_Name), 'w');% create a file
    fprintf(fileID,'%10s\t%18s\t%24s\t%20s\t%14s\t%12s\t%25s\r\n','Potential','Current', 'Current', '\g(D)F', '\g(D)F', 'd\g(D)F/dt', 'd\g(D)F/dt');% Name of the columns
    fprintf(fileID,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\r\n','V vs. Li/Li+', 'mA', 'mA', '', '','','');% units of the columns
    fprintf(fileID,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\r\n','', 'avg', 'std', 'avg', 'std','avg','std');% units of the columns
    fprintf(fileID,'%10.5f\t%19.10f\t%25.10f\t%21.10f\t%22.10f\t%20.10f\t%26.10f\r\n',E_I_F_dF_comp');% data exported
    fclose(fileID);% close the file


 %% delta_F vs Q
    % Calculate accumulated charge
    I_turnover = sort([strfind(Data(:,2)'>=0, [0 1]) strfind(Data(:,2)'>=0, [1 0])]); %the index right before turnover
    % separate datasheet into positive current and negative current region
    for i = 1:numel(I_turnover)+1
        if i == 1
           Data_sort_current{i} = Data(1:I_turnover(i),:);
        elseif i > 1 && i < numel(I_turnover)+1
           Data_sort_current{i} = Data(I_turnover(i-1)+1:I_turnover(i),:);
        else
           Data_sort_current{i} = Data(I_turnover(i-1)+1:end,:);
        end
    end
    % check data separation based on negative and positive current
    figure(3)
    for i = 1:numel(I_turnover)+1
        yyaxis left
        plot(Data_sort_current{i}(:,4),Data_sort_current{i}(:,2))
        yyaxis right
        plot(Data_sort_current{i}(:,4),Data_sort_current{i}(:,5))
        hold on
    end
    hold off
    
    % calculate charge
    for i = 1:numel(I_turnover)+1
        for j = 2:length(Data_sort_current{i}(:,2))
        	Q = abs(trapz(Data_sort_current{i}(1:j,4), Data_sort_current{i}(1:j,2)));% integration to find capacity
            Data_sort_current{i}(j,7) = Q;
        end
    end
    
    % create artificial arrays of Q
    Q_max = round(max(Data_sort_current{16}(:,7)),7);% use the anodic sweep of the 8th cycle, round to 7 decimal places
    Q_anodic_art = linspace(0,Q_max,350);
    Q_min = round(max(Data_sort_current{15}(:,7)),7);% use the cathodic sweep of the 8th cycle, round to 7 decimal places
    Q_cathodic_art = linspace(0,Q_min,350); % Q_min is positive because we used absolute value for integration.
    %create empty matrix to host F data
    F_Q_art_compile = nan(length(Q_anodic_art),(numel(I_turnover)+1)/2-3); % 7 columns of data
    F_Q_art_invert_compile = nan(length(Q_cathodic_art),(numel(I_turnover)+1)/2-3);% 7 columns of data
    % spline fit
    for i = 3:numel(I_turnover)-3 % from the pseudo-cathodic (0 current happends at different E) sweep of the 2nd cycle to the pseudo-anodic sweep of the 8th cycle
        Q_sort = Data_sort_current{i}(:,7);
        F_current_sort = Data_sort_current{i}(:,5);
        % spline fit
        F_Q_spline = spline(Q_sort',F_current_sort');
        % calculate fitted values
        if bitget(i,1)==0 %bitget gives 1 if it's odd and 0 if it's even
           F_Q_art = ppval(F_Q_spline, Q_anodic_art);
           % compilation
           F_Q_art_compile(:,i/2-1) = F_Q_art';% anodic cycle
       elseif bitget(i,1)==1
           F_Q_art = ppval(F_Q_spline,Q_cathodic_art);
           % compilation
           F_Q_art_invert_compile(:,(i-1)/2) = F_Q_art';% cathodic cycle
       end
    end
    
    % calculate averages
    F_Q_anodic_avg = mean(F_Q_art_compile(:,:),2); % 2nd ~ 8th cycles
    F_Q_cathodic_avg = mean(F_Q_art_invert_compile(:,:),2); % 2nd ~ 8th cycles
    % calculate standard deviations
    F_Q_anodic_std = std(F_Q_art_compile(:,:),0,2);
    F_Q_cathodic_std = std(F_Q_art_invert_compile(:,:),0,2);
    % put together one cycle
    Q_art_all = [-flip(Q_cathodic_art(5:end)');Q_anodic_art(2:end-1)'];
    F_Q_avg = [flip(F_Q_cathodic_avg(5:end));flip(F_Q_anodic_avg(2:end-1))];% removed the artifacts
    F_Q_std = [flip(F_Q_cathodic_std(5:end));flip(F_Q_anodic_std(2:end-1))];% removed the artifacts
    % compilation
    cath_filler = nan(length(Q_art_all)-length(Q_cathodic_art(1:end-1)),1);
    anod_filler = nan(length(Q_art_all)-length(Q_anodic_art(1:end-1)),1);
    Q_cath = [Q_cathodic_art(1:end-1)';cath_filler];
    F_Q_avg_cath = [F_Q_cathodic_avg(1:end-1);cath_filler];
    F_Q_std_cath = [F_Q_cathodic_std(1:end-1);cath_filler];
    Q_anod = [Q_anodic_art(1:end-1)';anod_filler];
    F_Q_avg_anod = [flip(F_Q_anodic_avg(1:end-1));anod_filler];
    F_Q_std_anod = [flip(F_Q_anodic_std(1:end-1));anod_filler];
    F_Q_comp = [Q_art_all,F_Q_avg,F_Q_std,...
        Q_cath,F_Q_avg_cath,F_Q_std_cath,Q_anod,F_Q_avg_anod,F_Q_std_anod];%assuming charges are all positive
        
    
    % text output
    Name = char(allNames(1));
    QF_Name = strrep(Name, Name(end-4:end), '_Q_F.txt');% create new file names
    fileID = fopen(char(QF_Name), 'w');% create a file
    fprintf(fileID,'%14s\t%12s\t%25s\t%14s\t%12s\t%25s\t%14s\t%12s\t%25s\r\n','Charge','\g(D)F', '\g(D)F','Charge','\g(D)F', '\g(D)F','Charge','\g(D)F', '\g(D)F');% Name of the columns
    fprintf(fileID,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\r\n','mAh', '','','mAh', '','','mAh', '','');% units of the columns
    fprintf(fileID,'%22.10f\t%20.10f\t%26.10f\t%22.10f\t%20.10f\t%26.10f\t%22.10f\t%20.10f\t%26.10f\r\n',F_Q_comp');% data exported
    fclose(fileID);% close the file

