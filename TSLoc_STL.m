 %% Basic rotation game for IPA lab experiments
% this should be a decent starting point for most of our experiments. Lots
% of timing information is being saved which MAY not be needed depending on
% experiment.
% PAB 6/16    
 
%%% Notes from OAK 9/2021
%   I received this file from NAF, who used it to replicate Hyosub Kim's
%   eLife 2019 paper effects (where there is a change in asymptotes when
%   participants switch from receiving a small error clamp + a hit or a
%   miss.
%   Started modifying from there.
 
%% Clear workspace and basic matlab setup
clear; close all; clc;
%Diary is used to keep track of errors or anything displayed in terminal
templogfile = 'temp_log.txt';
%clean old diary file if present
if ~isempty(dir(templogfile))
    delete(templogfile); 
end 
diary(templogfile)
addpath('Functions');
addpath('Sounds')
addpath('TargetFiles')
expcode='TSLoc_STL';

%% Basic PTB stuff
Priority(1);  
% hide mouse
HideCursor; 

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

%% Load sounds
% this is all a little complictated, but it is necessary for accurate audio
% timing. We are opening a single channel and then mixing the sounds
% through. Its reccommended for fast timing, no idea why.
InitializePsychSound(1);  


nrchannels = 2; %make stereo by default for simplicity
std_freq = 48000; %freq of master channel. all sounds converted to std.
% Its recommended to add 15 msecs latency on Windows, to protect against
% shoddy drivers.
% We have bought nice sound cards so not needed
sugLat = [];
% if IsWin
%     sugLat = 0.015;
% end

% Open real default [] soundcard as master device (+8) for playback only (+1), with
% standard low-latency, high timing precision mode, 2 channels, 48kHz:
pamaster = PsychPortAudio('Open', [], 1+8, 1, std_freq, nrchannels, [], sugLat);

% Start master immediately, wait for it to be started. We won't stop the
% master until the end of the session.
PsychPortAudio('Start', pamaster, 0, 0, 1);

% Create to slave audio devices for sound playback (+1), with same
% frequency, channel count et. as master.

%Load sounds and make some adjustments to them
[slow_sound, freq1] = audioread('too_slow.wav'); % too slow voice
if freq1 ~= std_freq % resample to 48khz
    slow_sound = resample(slow_sound, std_freq, freq1);
end
if size(slow_sound,2) == 1 %make stereo if not already
    slow_sound=[slow_sound';slow_sound'];
end
pa_tooslow = PsychPortAudio('OpenSlave', pamaster, 1);
PsychPortAudio('FillBuffer', pa_tooslow, slow_sound);

[miss_sound, freq2] = audioread('incorrect.wav'); % sound when miss the target
if freq2 ~= std_freq % resample to 48khz
    miss_sound = resample(miss_sound, std_freq, freq2);
end
if size(miss_sound,2) == 1 %make stereo if not already
    miss_sound=[miss_sound';miss_sound'];
end
pa_miss = PsychPortAudio('OpenSlave', pamaster, 1);
PsychPortAudio('FillBuffer', pa_miss, miss_sound);

[hit_sound, freq3] = audioread('correct.wav'); % hit
if freq3 ~= std_freq % resample to 48khz
    hit_sound = resample(hit_sound, std_freq, freq3);
end
if size(hit_sound,2) == 1 %make stereo if not already
    hit_sound=[hit_sound';hit_sound'];
end
pa_hit = PsychPortAudio('OpenSlave', pamaster, 1);
PsychPortAudio('FillBuffer', pa_hit, hit_sound);

[start_sound, freq4] = audioread('start.wav'); % game start
if freq4 ~= std_freq % resample to 48khz
    start_sound = resample(start_sound, std_freq, freq4);
end
if size(start_sound,2) == 1 %make stereo if not already
    start_sound=[start_sound';start_sound'];
else
    start_sound=start_sound';
end
pa_start = PsychPortAudio('OpenSlave', pamaster, 1);
PsychPortAudio('FillBuffer', pa_start, start_sound);

[neutral_sound, freq5] = audioread('neutral.wav'); % no performance feedback
if freq5 ~= std_freq % resample to 48khz
    start_sound = resample(dark_sound, std_freq, freq5);
end
if size(neutral_sound,2) == 1 %make stereo if not already
    neutral_sound=[neutral_sound';neutral_sound'];
else
    neutral_sound=neutral_sound';
end
pa_neutral = PsychPortAudio('OpenSlave', pamaster, 1);
PsychPortAudio('FillBuffer', pa_neutral, neutral_sound);

%% create subject's screen
bg_color = [0, 0, 0];
which_screen = 0;%0 is the primary screen
[window_ptr, screen_dimensions] = Screen(which_screen, 'OpenWindow', bg_color);%

Screen('BlendFunction', window_ptr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %suposed to reduce aliasing

%Setup keyboard loop to only check for space and exit key (@)
%Varies at least by operating system, possibly keyboard
RestrictKeysForKbCheck([32 50 160]); %50 = SPACE, 32 = 2/@, 160 = LEFT SHIFT

%% screen/tablet info % WILL VARY BY TABLET OR SCREEN
screenInfo = Screen('Resolution',0);
%using struct "c" to store the constants for easy saving and reading.
c.screenWidth = screenInfo.width; % width in px
c.screenHeight = screenInfo.height; % height in pix
% We need to figure out which display we are using. For reasons unknown the
% planar monitor returns the wrong size - oddly the same numbers as acer
% monitor. Made dir on planar computer called planar
prevdir=pwd;
if (exist('/Planar')>0);
    %we are using Planar monitor
    disp('Using Planar monitor settings');
    c.width_mm = 480; % Planar touch screen
    c.height_mm = 270; % Planar touch screen
    %load calib values
    cd /Planar
    RigCalib;
elseif (exist('/Acer')>0);
    %assuming using acer
    disp('Using Acer monitor settings');
    c.width_mm = 509; % Acer touch screen 
    c.height_mm = 286; % Acer touch screen
    %load calib values
    cd /Acer
    RigCalib;
else
    error('Error: UNKNOWN MONITOR. Missing folder on main drive for either Planar or Acer monitor.');
    c.width_mm = 480; % Planar touch screen
    c.height_mm = 270; % Planar touch screen
end
cd(prevdir);

c.pix_per_mm = (c.screenWidth/c.width_mm+c.screenHeight/c.height_mm)/2; % conversion ratio. not for all monitors...should be updated to have different x,y pix/mm
c.tab_length = 325; % tablet length mm.. %This will need to be updated based on the tablet being used. ATM all same size.
c.tab_height = 203; % tablet height mm..
c.tab_x_max = c.tab_length*100; % max tablet x values
c.tab_y_max = c.tab_height*100; % max tablet y values (wintabmex units are mm/100 or 10 ?m/wintab unit)
%% Standardized visual settings
c.start_x = c.screenWidth/2+c.offset_x; % middle of screen x
c.start_y = c.screenHeight/2+c.offset_y; % middle of screen y
c.tab_center_x = c.tab_x_max/2;
c.tab_center_y = c.tab_y_max/2;
c.cursor_radius = 7; % radius of cursor
 % larger than cursor
%visual start area is 4, but after find start, can be within 1.333 (or 16)
%of start. This reduces false starts in hold period.
c.start_area = 12; % radius of starting area 
c.start_fb_area = 40; %4 pix per mm
c.start_area_hold=c.start_area*(4/3);
green = [0,255,0];        % define green color
red = [200,50,10];        % define red color
black = [0,0,0];          % define black color
grey = [25, 25, 25];   % define grey color
magenta = [255, 0, 255];  % define magenta color
blue = [0,0,255]; % blue
white = [255,255,255]; % white
yellow = [200,200,0]; % yellow
orange = [255,100,0];
target_color = blue; % define target color
cursor_color = white; % cursor color

%% Constants for trial parameters
% These are the fixed task delays
c.hold_time = 0.5; % number of seconds until target appears
c.feedback_time = 0.05; % seconds for fb
c.reach_limit = .30; % max seconds allowed for reach
c.reach_min = .05; % min seconds allowed (not actually active right now)


%% Define instruction screens
IntroText = 'Slice through the target with your hand. \n For now, the cursor will show you where your hand is.';%'The goal of the game is to get your cursor on the target'; 
PreAimText = 'You will no longer see the cursor while you move to the target.';%'You will now receive feedback again'; %1
AimText = 'You will see the cursor again. However, it will not follow your hand movement. \n Reach anywhere except towards the target to see for yourself...'; %2
WashText = 'Do your best to move through the targets with your hand again.\nDo your best to ignore the cursor, which does not show your real hand position!'; %3

%% Enter subject name and info on experimenter screen. Check for trialfile
    ListenChar(2); % supress keystroke output to matlab command window
    FlushEvents('keyDown');
    Screen('Preference', 'TextRenderer', 0);
    Screen('TextSize', window_ptr, 24);
    sub_exp_date = GetEchoString(window_ptr, 'Enter Subject/Experiment/Date:',c.start_x-280,c.start_y-280,[255,255,255],black,0,2,[],[]);
    target_file_name = GetEchoString(window_ptr, 'Enter Target File:',c.start_x-280,c.start_y,[255,255,255],black,0,2,[],[]);
    FlushEvents('keyDown');
    Screen('Flip',window_ptr);
    filename = strcat('TargetFiles/',target_file_name,'.tgt'); % load subject-specific protocol
    loadtrialfile = true;
    %Previous version made trial file, I prefer making them in a batch to
    %counterbalance across participants
while loadtrialfile 
    %check if file exists in current directory
    if ~isempty(dir(filename))
        trial_filetemp = importdata(filename);
        loadtrialfile = false;
    else
        target_file_name = GetEchoString(window_ptr, 'Incorrect Target File Try Again:',c.start_x-280,c.start_y,[255,255,255],black,0,2,[],[]);
        FlushEvents('keyDown');
        Screen('Flip',window_ptr);
        filename = strcat('TargetFiles/',target_file_name,'.tgt'); % load subject-specific protocol
    end
end
%Check if data file already exists and append "_B" to new sub name if it
%does. Otherwise data would be overwritten :(
tempsavename = strcat('Data/',sub_exp_date,'_',target_file_name,'_',expcode,'.mat');
if exist(tempsavename,'file')
    disp('File name already exists, appending "_B" to sub initials')
    sub_exp_date = strcat(sub_exp_date,'_B');
end

%% Get Trial Protocol Info

% determine column # of wanted trialfile columns and get data
targetheaders = strtrim(trial_filetemp.colheaders);
trial_file = trial_filetemp.data;
targetfile.name = target_file_name;
targetfile.endpoint_fb=trial_file(:,strcmp('endpoint_feedback',targetheaders));
targetfile.online_fb=trial_file(:,strcmp('online_feedback',targetheaders));
targetfile.aim_report=trial_file(:,strcmp('aim_report',targetheaders));
targetfile.aim_ring=trial_file(:,strcmp('aim_ring',targetheaders));
targetfile.rotation=trial_file(:,strcmp('rotation',targetheaders));
targetfile.instruction=trial_file(:,strcmp('instruction',targetheaders));
targetfile.target_angle=trial_file(:,strcmp('target_angle',targetheaders));
targetfile.target_dist=trial_file(:,strcmp('target_distance',targetheaders));
targetfile.jump_onset=trial_file(:,strcmp('jump_onset',targetheaders));
targetfile.jump_radius=trial_file(:,strcmp('jump_radius', targetheaders));
targetfile.jump_angle=trial_file(:,strcmp('jump_angle', targetheaders));
targetfile.hollow = trial_file(:,strcmp('hollow', targetheaders));
targetfile.size=trial_file(:,strcmp('size',targetheaders));
targetfile.distractor= trial_file(:,strcmp('distractor', targetheaders));

num_Trials = length(targetfile.target_angle);

% Targets and Reach Distance
% target angles (entry 1) and reach distance (entry 2)
polar_target_loc = [targetfile.target_angle targetfile.target_dist]; % get target locations
% RING RADIUS/REACH DISTANCE
ring_radius = polar_target_loc(:,2)*c.pix_per_mm; % convert to pixels
%To do accurately need to conver in to x and y components and then

% quick loop to make target locations cartesian
target_loc = zeros(num_Trials,2); % init
target_jump_loc = zeros(num_Trials, 2);
jumpinplace = zeros(num_Trials, 1);
for k = 1:num_Trials
    targAng = polar_target_loc(k,1)*pi/180; % chang targ ang to radians
    targRad = ring_radius(k); % reach distance
    target_loc(k,1) = (targRad*cos(targAng))+c.start_x; % target x
    target_loc(k,2) = (targRad*-sin(targAng))+c.start_y; % target y
    if(abs(mod(targetfile.jump_angle(k), 1)) == 0.01)
        jumpinplace(k) = 1;
    end     
    tjAng = (polar_target_loc(k,1)-targetfile.jump_angle(k)) * pi/180;
    tjRad = targetfile.jump_radius(k)*c.pix_per_mm;
    target_jump_loc(k,1) = (tjRad*cos(tjAng))+c.start_x; % target x
    target_jump_loc(k,2) = (tjRad*-sin(tjAng))+c.start_y; % target y
end
% aiming ring visuals
ring_dimensions = [c.start_x-ring_radius, c.start_y-ring_radius, c.start_x+ring_radius, c.start_y+ring_radius];
ring_color = blue;
ring_width = 2;

%% initialize large arrays
% this takes a little while so doing here
% fastest to initialize large arrays
% Should maybe be changed so one long data stream, 
%   but when have done so game locks up...
large_array_size = 5e4; % at the moment anything larger makes game slow
hand_x = nan(num_Trials,large_array_size);
hand_y = nan(num_Trials,large_array_size);
hand_event_idx = nan(num_Trials,large_array_size);
hand_time_tab = nan(num_Trials,large_array_size);
hand_time = nan(num_Trials,large_array_size);
fb_hand_x = nan(num_Trials,large_array_size); %this is redundant with hand_x but transformed to pixels and only last sample per loop
fb_hand_y = nan(num_Trials,large_array_size); %this is redundant with hand_x but transformed to pixels and only last sample per loop
cursor_x_save = nan(num_Trials,large_array_size);
cursor_y_save = nan(num_Trials,large_array_size);
event_idx = nan(num_Trials,large_array_size);
experiment_time = nan(num_Trials,large_array_size);
%% TABLET init - this is slow so doing up here
try
    WinTabMex(0, window_ptr);
catch
    WinTabMex(1); %if get error try closing wintab
    WinTabMex(0, window_ptr);
end

%% DATA INTITIALIZATION
% init small/basic data struct
data = struct('cursor_loc_final',nan(num_Trials,2),'hand_loc_final',nan(num_Trials,2),'aim_raw',nan(num_Trials,2),...
    'target_loc',nan(num_Trials,2),'rt',nan(num_Trials,1),'movement_time',nan(num_Trials,1),...
    'too_slow',nan(num_Trials,1),'target_hit',nan(num_Trials,1),...
    'rotation',nan(num_Trials,1),'instruction',nan(num_Trials,1),'aim_report',nan(num_Trials,1),'aim_ring',nan(num_Trials,1),...
    'online_feedback',nan(num_Trials,1),'endpoint_feedback',nan(num_Trials,1),'target_angle',nan(num_Trials,1),'target_dist',nan(num_Trials,1));
%using separate struct for timers
timers = struct('hold_start',nan(num_Trials,1),'hold_end',nan(num_Trials,1),'aim_onset',nan(num_Trials,1),'aim_report',nan(num_Trials,1),'go',nan(num_Trials,1),'move_start',nan(num_Trials,1),'move_end',nan(num_Trials,1),'FB_onset',nan(num_Trials,1),'FB_offset',nan(num_Trials,1));
%% Timers
t_hold_start = NaN;
t_hold_end = NaN;
t_aim_onset = NaN;
t_aim_report = NaN;
t_go = NaN;
t_move_start = NaN;
t_move_end = NaN;
t_FB_onset = NaN;
t_FB_offset = NaN;

%% Init Events & Time
e_find_start = 1; % first event
e_hold_start = 0;
e_aiming = 0;
e_go_signal = 0;
e_reaching = 0;
e_ring_broken = 0;
e_feedback = 0;
e_data = 0;
e_too_slow = 0;
e_aim_reminder = 0;
aim_x = 0; % init aim
aim_y = 0; % init aim
axisStateX=0;
axisStateY=0;
dist_last=0;
dist_max=0;
x=0;
y=0;
% init loop for tracking missed tablet data
% init time/counter
count = 0; % start indexing counter
handcount = 0;%start at zero and increase with each hand sample
trialnum = 1; % commence trials
trialtoolong=false;
%Tracking false starts
falsestarts=[];
falsestart_counts=0;
refindstart=false;
play_too_slow = false;

%% Show intro text
DrawFormattedText(window_ptr,IntroText,'center','center',white);
Screen('Flip',window_ptr, [], 0);
PsychPortAudio('Start', pa_start, 1, 0, 1); % wake up experimenter!
KbWait(-1);

%% Make sure pen is on tablet and sync time
WinTabMex(2); %Empties the packet queue in preparation for collecting actual data
%needs cleaing up
%start_time1 = GetSecs; 
while 1    
    pkt = WinTabMex(5);
    if ~isempty(pkt)
        start_timea = GetSecs;
        WinTabMex(2); %Clear to sync time
        start_timeb = GetSecs;
        pkt = WinTabMex(5);
        if ~isempty(pkt)
            break
        else
            Screen('DrawText',window_ptr,'Put Pen on Tablet',c.start_x-100,c.start_y-50,red);
            Screen('Flip',window_ptr, [], 0);
        end
    else
        Screen('DrawText',window_ptr,'Put Pen on Tablet',c.start_x-100,c.start_y-50,red);
        Screen('Flip',window_ptr, [], 0);
    end
end
start_time = mean([start_timea start_timeb]);
c.start_time = start_time;
% start time
start_time_tablet = pkt(6); %can do here because already timestamped
tabdata = pkt;
t_experiment_time_last=0; %start_time - start_time
WinTabMex(2); %Clear to start actual data collection
NewTabData = false;
try
    %% MAIN TRIAL LOOP %%
    
    while trialnum <= num_Trials
%         
        
        %Get tablet data
        tabdata = GetTabletData();
        %Update hand indices
        if size(tabdata,2)>1
            newhandi=handcount+1:handcount+size(tabdata,2);%its easier to figure out the new indices here
            count=count+1;
        elseif size(tabdata,2)==1
            newhandi=handcount+1;
            count=count+1;
        end      

        if ~isempty(tabdata)
            % clock loop time
            c.target_radius = targetfile.size(trialnum);
            t_experiment_time = GetSecs-start_time;
            dt=t_experiment_time-t_experiment_time_last;
            NewTabData = true;
            hand_x(trialnum,newhandi)=tabdata(1,:);
            hand_y(trialnum,newhandi)=tabdata(2,:);
            hand_time(trialnum,newhandi) = t_experiment_time_last+((1:length(newhandi))*dt/length(newhandi));
            hand_time_tab(trialnum,newhandi)=(tabdata(6,:)/1000)-start_time_tablet; %Not accurate, has  both offset and drift from comp clock
                        %save last xy in case nonewdata
            if count<2 %because event has not been set yet default to 0
                hand_event_idx(trialnum,newhandi)=0;%take last evemt as no event yet
            else
                hand_event_idx(trialnum,newhandi)=event_idx(trialnum,count-1);%take last evemt as no event yet
            end
            % update large movement data matrix with hand positions
            handcount=handcount+length(newhandi);
            %For screen drawing purposes just take last
            axisStateX = (((tabdata(1,end)-c.tab_center_x)/100)*c.pix_per_mm)+c.start_x;
            axisStateY = (((c.tab_y_max-(tabdata(2,end))-c.tab_center_y)/100)*c.pix_per_mm)+c.start_y;
            tab_count = newhandi;
            [dist_max, dist_last] = RadialDist(hand_x(trialnum,tab_count),hand_y(trialnum,tab_count),c.tab_center_x,c.tab_center_y,c.pix_per_mm);
            %hand_ang = acosd((axisStateX - c.start_x)./dist_max) - targetfile.target_angle(trialnum);
            if (e_find_start == 1) || (abs(targetfile.rotation(trialnum)) == 0)
                cursor_x = axisStateX;
                cursor_y = axisStateY;
            elseif abs(targetfile.rotation(trialnum)) > 0
                rot_matrix = [cos(targetfile.rotation(trialnum)*pi/180) sin(targetfile.rotation(trialnum)*pi/180); -sin(targetfile.rotation(trialnum)*pi/180) cos(targetfile.rotation(trialnum)*pi/180)]; % rot mat
                 hand_r = sqrt((axisStateX-c.start_x)^2+ (axisStateY-c.start_y)^2);
                 targ_ang = targetfile.target_angle(trialnum)*pi/180;
                 rotated = [hand_r*cos(-targ_ang) hand_r*sin(-targ_ang)]*rot_matrix;
                 rotatedX = rotated(1)+c.start_x;
                 rotatedY = rotated(2)+c.start_y;
                 cursor_x = rotatedX;
                 cursor_y = rotatedY;
%                cursor_x = axisStateX;
%               cursor_y = axisStateY;
            end
            fb_hand_x(trialnum,count) = axisStateX; % index hand position x
            fb_hand_y(trialnum,count) = axisStateY; % index hand position y
            cursor_x_save(trialnum,count) = cursor_x; % index hand position x
            cursor_y_save(trialnum,count) = cursor_y; % index hand position y
            experiment_time(trialnum,count)=t_experiment_time;
            time_dt(trialnum,count)=dt;%remove
            t_experiment_time_last = t_experiment_time;
        elseif NewTabData
            NewTabData = false; %clear as no current data
        end
        
        
        %% FIND START
        if e_find_start == 1
            if NewTabData
                event_idx(trialnum,count) = 1;
                Screen('DrawText',window_ptr,num2str(trialnum),30,25,white);
                if e_aim_reminder 
                     Screen('DrawText',window_ptr,'Remember to report aim',c.start_x-200,c.start_y-50,red);
                end
                
                if dist_max <= c.start_fb_area% && targetfile.online_fb(trialnum) == 1%Give online feedback
                     cursor = [cursor_x-c.cursor_radius, cursor_y-c.cursor_radius, cursor_x+c.cursor_radius, cursor_y+c.cursor_radius];
                     Screen(window_ptr, 'FillOval', cursor_color, cursor);
                     %Show start
                     Screen(window_ptr, 'FrameOval', white, [c.start_x-c.start_area,c.start_y-c.start_area,c.start_x+c.start_area,c.start_y+c.start_area]);
                else
                     % find start helping circle - we don't want to show cursor FB here
                     find_circle = [c.start_x-dist_last,c.start_y-dist_last,c.start_x+dist_last,c.start_y+dist_last];
                     if refindstart %if they false started make find center circle red
                          Screen('FrameOval',window_ptr, red, find_circle);
                     else
                           Screen('FrameOval',window_ptr, white, find_circle);
                     end
                end
             Screen('Flip',window_ptr,[],[],2);%dont wait for sync
                % is subject in start area?
             if dist_last < c.start_area
                  e_find_start = 0;
                  e_hold_start = 1;
                  t_hold_start = GetSecs-start_time;
                  if dist_max>dist_last
                      dist_max=dist_last;
                  end
            end
         end
        end
        
        %% HOLD START
        if e_hold_start == 1
            if e_aim_reminder 
                 Screen('DrawText',window_ptr,'Remember to report aim',c.start_x-200,c.start_y-50,red);
            end
%             if targetfile.online_fb(trialnum) == 1%Give online feedback
%                     cursor = [cursor_x-c.cursor_radius, cursor_y-c.cursor_radius, cursor_x+c.cursor_radius, cursor_y+c.cursor_radius];
%                     Screen(window_ptr, 'FillOval', cursor_color, cursor);
%                     %Show start
%                     Screen(window_ptr, 'FrameOval', white, [c.start_x-c.start_area,c.start_y-c.start_area,c.start_x+c.start_area,c.start_y+c.start_area]);
%             else
                % show centered white oval
                Screen(window_ptr, 'FillOval', white, [c.start_x-c.start_area,c.start_y-c.start_area,c.start_x+c.start_area,c.start_y+c.start_area]);
              %end
            Screen('Flip',window_ptr,[],[],2);%dont wait for synch            
            % holding at start?
            if NewTabData
                event_idx(trialnum,count) = 2;
                if dist_max > c.start_area_hold % they left the start position
                    e_find_start = 1;
                    e_hold_start = 0;
                    %Make sure have actually held for half of hold time
                    %Otherwise gets triggered while finding start
                    leavetime = GetSecs-start_time;
                    if leavetime-t_hold_start >= (c.hold_time/2)
                        refindstart=true;
                        falsestart_counts=falsestart_counts+1;
                        falsestarts(falsestart_counts).mi = trialnum;
                        falsestarts(falsestart_counts).exp_time = leavetime;
                        falsestarts(falsestart_counts).event = event_idx(trialnum,count);
                    end
                end
            end
            % we want timming outside the new data statement so that timing
            % is more precise
            % did subject hold long enough?
            hold_count = GetSecs-start_time;
            if hold_count-t_hold_start >= c.hold_time;
                t_hold_end = GetSecs-start_time;
                t_real_hold_time=t_hold_start-t_hold_start;
                %Show start
%                 if targetfile.online_fb(trialnum) == 1%Give online feedback
%                     cursor = [cursor_x-c.cursor_radius, cursor_y-c.cursor_radius, cursor_x+c.cursor_radius, cursor_y+c.cursor_radius];
%                     Screen(window_ptr, 'FillOval', cursor_color, cursor);
%                     %Show start
%                     Screen(window_ptr, 'FrameOval', white, [c.start_x-c.start_area,c.start_y-c.start_area,c.start_x+c.start_area,c.start_y+c.start_area]);
%                 else
                    % show centered white oval
                 Screen(window_ptr, 'FillOval', white, [c.start_x-c.start_area,c.start_y-c.start_area,c.start_x+c.start_area,c.start_y+c.start_area]);
%                 end                % draw aim ring target
                if targetfile.aim_report(trialnum) == 1 && targetfile.aim_ring(trialnum) == 1
                    Screen('FrameOval', window_ptr, ring_color, ring_dimensions(trialnum,:), ring_width);       
                end
                target = [target_loc(trialnum,1)-c.target_radius, target_loc(trialnum,2)-c.target_radius, target_loc(trialnum,1)+c.target_radius, target_loc(trialnum,2)+c.target_radius];
                targ_jump = [target_jump_loc(trialnum,1)-c.target_radius, target_jump_loc(trialnum,2)-c.target_radius, target_jump_loc(trialnum,1)+c.target_radius, target_jump_loc(trialnum,2)+c.target_radius];
                if(targetfile.hollow(trialnum) == 1)
                    Screen(window_ptr, 'FrameOval', blue, target, 2);
                else
                    Screen(window_ptr, 'FillOval', blue, target);
                end
                
                if e_aim_reminder 
                    Screen('DrawText',window_ptr,'Remember to report aim',c.start_x-200,c.start_y-50,red);
                end
                
                [~, t_aim_onset, ~, ~, ~]=Screen('Flip',window_ptr);%dont wait for synch, keep on screan
                t_aim_onset=t_aim_onset-start_time;
                % reset mouse for new aim
                if targetfile.aim_report(trialnum) == 1
                    SetMouse(0,0,window_ptr);
                end
                
                e_hold_start = 0;
                e_aiming = 1; % time to aim!
            end
        end
        
        
        %% AIM SELCTION / SHOW STIMULI
        if e_aiming == 1
            if NewTabData
                event_idx(trialnum,count) = 3;
                % subject starts reaching before an "acceptable" aim is chosen
                if dist_max > c.start_area_hold
                    refindstart=true;
                    falsestart_counts=falsestart_counts+1;
                    falsestarts(falsestart_counts).mi = trialnum;
                    falsestarts(falsestart_counts).exp_time = GetSecs-start_time;
                    falsestarts(falsestart_counts).event = event_idx(trialnum,count);
                    if targetfile.aim_report(trialnum) == 1
                        Screen('DrawText',window_ptr,'Remember to report aim',c.start_x-200,c.start_y-50,red);
                        Screen('Flip',window_ptr, [], 1);
                        e_aim_reminder = 1;
                    end
                    e_aiming = 0;
                    e_find_start = 1;
                end
            end
            
            % CALCULATE AIM CHOICE
            % NOTE: unfortunately, program will read tablet as mouse input as well, so you have to use
            % a hack to ignore any entries that are consonant with their start point tablet position.
            [mouse_x,mouse_y] = GetMouse;
            mouse_dist_fromStart = sqrt(((mouse_x-c.start_x)^2)+((mouse_y-c.start_y)^2));
            
            if targetfile.aim_report(trialnum) == 1 && mouse_dist_fromStart > (c.start_area*2) % only record/animate mouse position if it's not in/near the start area
                [aim_x,aim_y] = GetMouse; % log
                aim_ang = atan2((aim_y-c.start_y),(aim_x-c.start_x)); % angle of aim
                aim_distance = sqrt(((aim_x-c.start_x)^2)+((aim_y-c.start_y)^2)); % distance of aim from center point
                
                if aim_distance < (ring_radius(trialnum)*2) && aim_distance > (ring_radius(trialnum)/2) % if aim is acceptably close to the ring, show it and move on
                    t_aim_report=GetSecs-start_time; %could move up in loop for slightly more accurate tim
                    aim_marker = [aim_x-c.cursor_radius, aim_y-c.cursor_radius, aim_x+c.cursor_radius, aim_y+c.cursor_radius];
                    if targetfile.online_fb(trialnum) == 1%Give online feedback
                        cursor = [cursor_x-c.cursor_radius, cursor_y-c.cursor_radius, cursor_x+c.cursor_radius, cursor_y+c.cursor_radius];
                        Screen(window_ptr, 'FillOval', cursor_color, cursor);
                        %Show start
                        Screen(window_ptr, 'FrameOval', white, [c.start_x-c.start_area,c.start_y-c.start_area,c.start_x+c.start_area,c.start_y+c.start_area]);
                    else
                        % show centered white oval
                        Screen(window_ptr, 'FillOval', white, [c.start_x-c.start_area,c.start_y-c.start_area,c.start_x+c.start_area,c.start_y+c.start_area]);
                    end
                    % target
                    if(targetfile.hollow(trialnum) == 1)
                        Screen(window_ptr, 'FrameOval', blue, target, 2);
                    else
                        Screen(window_ptr, 'FillOval', blue, target);
                    end
                    % draw aim ring target
%                     if targetfile.aim_ring(trialnum) == 1
%                         Screen('FrameOval', window_ptr, ring_color, ring_dimensions(trialnum,:), ring_width);       
%                     end
                    % No Cross Hair aim mark
%                     % draw aim
%                     Screen('DrawLine',window_ptr, white, aim_x, aim_y-8,aim_x,aim_y+8,4);
%                     Screen('DrawLine',window_ptr, white, aim_x-8, aim_y,aim_x+8,aim_y,4);
                    
                    Screen('Flip',window_ptr,[],1,2);%dont wait for synch, keep on screan
                    e_aim_reminder = 0;
                    e_aiming = 0;
                    e_go_signal = 1;
                    firstloop = true;
                end
                
            elseif targetfile.aim_report(trialnum) == 0
                aim_ang = NaN;
                aim_raw = [NaN NaN];
                e_aiming = 0;
                e_go_signal = 1;
                firstloop = true;
                t_aim_report = NaN;
            end
            
        end
        
                
        %% "GO" CUE
        if e_go_signal == 1
           event_idx(trialnum,count) = 5;
            if firstloop == true
                            %Screen('Flip',window_ptr);
    %             % draw aim ring target
    %             if targetfile.aim_ring(trialnum) == 1
    %                 Screen('FrameOval', window_ptr, ring_color, ring_dimensions(trialnum,:), ring_width);       
    %             end
                    if targetfile.online_fb(trialnum) == 1%Give online feedback
                        cursor = [cursor_x-c.cursor_radius, cursor_y-c.cursor_radius, cursor_x+c.cursor_radius, cursor_y+c.cursor_radius];
                        Screen(window_ptr, 'FillOval', cursor_color, cursor);
                        %Show start
                        Screen(window_ptr, 'FillOval', white, [c.start_x-c.start_area,c.start_y-c.start_area,c.start_x+c.start_area,c.start_y+c.start_area]);
                    else
                        % show centered white oval
                        Screen(window_ptr, 'FillOval', white, [c.start_x-c.start_area,c.start_y-c.start_area,c.start_x+c.start_area,c.start_y+c.start_area]);
                    end
                    if dist_last >= targetfile.jump_onset*c.pix_per_mm
                        if(jumpinplace(trialnum) == 0)
                            if(targetfile.distractor(trialnum) == 1)
                                Screen(window_ptr, 'FillOval', target_color, target);
                                Screen(window_ptr, 'FillOval', target_color, targ_jump);
                            else
                                if(targetfile.hollow(trialnum) == 1)
                                    Screen(window_ptr, 'FrameOval', target_color, targ_jump, 2);
                                else
                                    Screen(window_ptr, 'FillOval', target_color, targ_jump);
                                end
                            end
                        else
                            [~, onsetTime, ~, ~, ~]=Screen('Flip', window_ptr);
                            Screen(window_ptr, 'FillOval', cursor_color, cursor);
                            WaitSecs('UntilTime', onsetTime + 0.02);
                            jumpinplace(trialnum) = 0;
                        end
                    else
                        if(targetfile.hollow(trialnum) == 1)
                            Screen(window_ptr, 'FrameOval', target_color, target, 2);
                        else
                            Screen(window_ptr, 'FillOval', target_color, target);
                        end
                    end
    %             if targetfile.cross_hairs(trialnum) == 1
    % %             % draw aim
    %                 Screen('DrawLine',window_ptr, white, aim_x, aim_y-8,aim_x,aim_y+8,4);
    %                 Screen('DrawLine',window_ptr, white, aim_x-8, aim_y,aim_x+8,aim_y,4);
    %             end
                [~, t_go, ~, ~, ~]=Screen('Flip',window_ptr);%wait for sync
                t_go=t_go-start_time;
                if dist_last > c.start_area_hold
                    e_go_signal = 0;
                    e_reaching = 1;
                    %update to corrected hand time
                    t_move_start = DistTime(hand_x(trialnum,tab_count),hand_y(trialnum,tab_count),hand_time(trialnum,tab_count),c.tab_center_x,c.tab_center_y,c.pix_per_mm,c.start_area);
                    RT = t_move_start - t_go; % log subject's RT
                end
                firstloop = false;
            else
                 %Screen('Flip',window_ptr);
    %             % draw aim ring target
    %             if targetfile.aim_ring(trialnum) == 1
    %                 Screen('FrameOval', window_ptr, ring_color, ring_dimensions(trialnum,:), ring_width);       
    %             end
                    if targetfile.online_fb(trialnum) == 1%Give online feedback
                        cursor = [cursor_x-c.cursor_radius, cursor_y-c.cursor_radius, cursor_x+c.cursor_radius, cursor_y+c.cursor_radius];
                        Screen(window_ptr, 'FillOval', cursor_color, cursor);
                        %Show start
                        Screen(window_ptr, 'FrameOval', white, [c.start_x-c.start_area,c.start_y-c.start_area,c.start_x+c.start_area,c.start_y+c.start_area]);
                    else
                        % show centered white oval
                        Screen(window_ptr, 'FillOval', white, [c.start_x-c.start_area,c.start_y-c.start_area,c.start_x+c.start_area,c.start_y+c.start_area]);
                    end
                    if dist_last >= targetfile.jump_onset*c.pix_per_mm
                        if(jumpinplace(trialnum) == 0)
                            if(targetfile.distractor(trialnum) == 1)
                                Screen(window_ptr, 'FillOval', target_color, target);
                                Screen(window_ptr, 'FillOval', target_color, targ_jump);
                            else
                                if(targetfile.hollow(trialnum) == 1)
                                    Screen(window_ptr, 'FrameOval', target_color, targ_jump, 2);
                                else
                                    Screen(window_ptr, 'FillOval', target_color, targ_jump);
                                end
                            end
                        else
                            [~, onsetTime, ~, ~, ~]=Screen('Flip', window_ptr);
                            Screen(window_ptr, 'FillOval', cursor_color, cursor);
                            WaitSecs('UntilTime', onsetTime + 0.02);
                            jumpinplace(trialnum) = 0;
                        end
                    else
                        if(targetfile.hollow(trialnum)==1)
                            Screen(window_ptr, 'FrameOval', target_color, target, 2);
                        else
                            Screen(window_ptr, 'FillOval', target_color, target);
                        end
                    end
    %             if targetfile.cross_hairs(trialnum) == 1
    % %             % draw aim
    %                 Screen('DrawLine',window_ptr, white, aim_x, aim_y-8,aim_x,aim_y+8,4);
    %                 Screen('DrawLine',window_ptr, white, aim_x-8, aim_y,aim_x+8,aim_y,4);
    %             end
                Screen('Flip',window_ptr,[],[],2);%dont wait for sync
                if dist_last > c.start_area_hold
                    t_flip_start = t_experiment_time;
                    e_go_signal = 0;
                    e_reaching = 1;
                    %update to corrected hand time
                    t_move_start = DistTime(hand_x(trialnum,tab_count),hand_y(trialnum,tab_count),hand_time(trialnum,tab_count),c.tab_center_x,c.tab_center_y,c.pix_per_mm,c.start_area );
                    RT = t_move_start - t_go; % log subject's RT
                end
            end
         end
        
        
        %% REACHING
        if e_reaching == 1
            if NewTabData
                event_idx(trialnum,count) = 6;
    %             % draw aim ring target
    %             if targetfile.aim_ring(trialnum) == 1
    %                 Screen('FrameOval', window_ptr, ring_color, ring_dimensions(trialnum,:), ring_width);       
    %             end
                %draw start
                Screen(window_ptr, 'FrameOval', white, [c.start_x-c.start_area,c.start_y-c.start_area,c.start_x+c.start_area,c.start_y+c.start_area]);

                % target
                if dist_last >= targetfile.jump_onset*c.pix_per_mm
                        if(jumpinplace(trialnum) == 0)
                            if(targetfile.distractor(trialnum) == 1)
                                Screen(window_ptr, 'FillOval', target_color, target);
                                Screen(window_ptr, 'FillOval', target_color, targ_jump);
                            else
                                if(targetfile.hollow(trialnum) == 1)
                                    Screen(window_ptr, 'FrameOval', target_color, targ_jump, 2);
                                else
                                    Screen(window_ptr, 'FillOval', target_color, targ_jump);
                                end
                            end
                        else
                            [~, onsetTime, ~, ~, ~]=Screen('Flip', window_ptr);
                            Screen(window_ptr, 'FillOval', cursor_color, cursor);
                            WaitSecs('UntilTime', onsetTime + 0.02);
                            jumpinplace(trialnum) = 0;
                        end
                else
                    if(targetfile.hollow(trialnum)==1)
                       Screen(window_ptr, 'FrameOval', target_color, target, 2);
                    else
                       Screen(window_ptr, 'FillOval', target_color, target);
                    end
                end
    %             if targetfile.cross_hairs(trialnum) == 1
    % %             % draw aim
    %                 Screen('DrawLine',window_ptr, white, aim_x, aim_y-8,aim_x,aim_y+8,4);
    %                 Screen('DrawLine',window_ptr, white, aim_x-8, aim_y,aim_x+8,aim_y,4);
    %             end

                % online FB?
                if targetfile.online_fb(trialnum) == 1 && dist_last <= ring_radius(trialnum) % online FB within ring
                    cursor = [cursor_x-c.cursor_radius, cursor_y-c.cursor_radius, cursor_x+c.cursor_radius, cursor_y+c.cursor_radius];
                    Screen(window_ptr, 'FillOval', cursor_color, cursor);
                elseif targetfile.online_fb(trialnum) == 1 && dist_last > ring_radius(trialnum)
                    % end cursor angle
                    cursor_ang = atan2((cursor_y-c.start_y),(cursor_x-c.start_x)); %angle of cursor
                    fb_x = (ring_radius(trialnum)*cos(cursor_ang))+c.start_x;
                    fb_y =  (ring_radius(trialnum)*sin(cursor_ang))+c.start_y;
                    cursor = [fb_x-c.cursor_radius, fb_y-c.cursor_radius, fb_x+c.cursor_radius, fb_y+c.cursor_radius];
                    Screen(window_ptr, 'FillOval', cursor_color, cursor);
                end

                Screen('Flip',window_ptr,[],[],2);%dont wait for synch

                % Broken ring?
                if dist_max >= ring_radius(trialnum) % reach complete
                    e_reaching = 0;
                    e_feedback = 1;
                    t_move_end = DistTime(hand_x(trialnum,tab_count),hand_y(trialnum,tab_count),hand_time(trialnum,tab_count),c.tab_center_x,c.tab_center_y,c.pix_per_mm,ring_radius(trialnum));
                    movement_time = t_move_end - t_move_start;

                    %% FEEDBACK!

                    % feedback info
                    % end cursor angle
                    cursor_ang = atan2((cursor_y-c.start_y),(cursor_x-c.start_x)); %angle of cursor
                    fb_x = (ring_radius(trialnum)*cos(cursor_ang))+c.start_x;
                    fb_y =  (ring_radius(trialnum)*sin(cursor_ang))+c.start_y;
                    % end hand angle
                    reach_ang = atan2((axisStateY-c.start_y),(axisStateX-c.start_x)); % angle of hand reach
                    hand_ep_x = (ring_radius(trialnum)*cos(reach_ang))+c.start_x;
                    hand_ep_y = (ring_radius(trialnum)*sin(reach_ang))+c.start_y;

                    %% FB visuals

    %                 % draw aim ring target
    %                 if targetfile.aim_ring(trialnum) == 1
    %                     Screen('FrameOval', window_ptr, ring_color, ring_dimensions(trialnum,:), ring_width);       
    %                 end
                    %draw start
                    Screen(window_ptr, 'FrameOval', white, [c.start_x-c.start_area,c.start_y-c.start_area,c.start_x+c.start_area,c.start_y+c.start_area]);

                    % Sounds
                    % "too slow" sound
                    if movement_time >= 2*c.reach_limit % over 600 ms
                        e_too_slow = 1;
                        play_too_slow = true;
                    end
                    % with a low-latency sound card we could beat the
                    % visual FB. Move later if so.
                    if targetfile.endpoint_fb(trialnum) == 0 && targetfile.online_fb(trialnum) == 0; % no FB trials
                        fb_color = blue;
                        %PsychPortAudio('Start', pa_neutral, 1, 0, 1); % neutral
                        hit = 0;
                    elseif targetfile.rotation(trialnum)~=0
                        fb_color = blue;
                        %PsychPortAudio('Start', pa_neutral, 1, 0, 1); % neutral
                        hit = 1; % this is not right - they may have missed here depending on the rotation and target size
                    elseif sqrt((fb_x-target_loc(trialnum,1))^2 + (fb_y-target_loc(trialnum,2))^2) <= c.target_radius+1
                        WaitSecs(0.3);
                        %PsychPortAudio('Start', pa_hit, 1, 0, 1); % success!
                        fb_color = blue;
                        hit = 1;
                    else
                        WaitSecs(0.3);
                        %PsychPortAudio('Start', pa_miss, 1, 0, 1); % missed
                        fb_color = blue;
                        hit = 0;
                    end
                    sound_start = GetSecs;

                    % Target
                        if(targetfile.hollow(trialnum) == 1)
                            Screen(window_ptr, 'FrameOval', fb_color, targ_jump, 2);
                        else
                            if(targetfile.distractor(trialnum) == 1)
                                Screen(window_ptr, 'FillOval', target_color, targ_jump);
                            else
                                Screen(window_ptr, 'FillOval', fb_color, targ_jump);
                            end
                        end
    %                 if targetfile.cross_hairs(trialnum) == 1
    %     %             % draw aim
    %                     Screen('DrawLine',window_ptr, white, aim_x, aim_y-8,aim_x,aim_y+8,4);
    %                     Screen('DrawLine',window_ptr, white, aim_x-8, aim_y,aim_x+8,aim_y,4);
    %                 end

                    % ep cursor fb
                    if targetfile.endpoint_fb(trialnum) == 1
                        fb_cursor = [fb_x-c.cursor_radius, fb_y-c.cursor_radius, fb_x+c.cursor_radius, fb_y+c.cursor_radius];
                        Screen(window_ptr, 'FillOval', cursor_color, fb_cursor);
                    end
                    [~, t_FB_onset, ~, ~, ~]=Screen('Flip',window_ptr,[],1);% wait for sync, keep on screen
                    t_FB_onset=t_FB_onset-start_time; %IF USING ONLINE FEEDBACK AND NOT EP THIS WILL NOT BE ACCURATE
                    fb_start = GetSecs; % start fb clock

                end
            end
        end
        
        %% FEEDBACK loop to continue visuals over length of fb time
        if e_feedback == 1;
            
            event_idx(trialnum,count) = 7;
            fb_hold = GetSecs;
            
           if play_too_slow && e_too_slow %play too slow 800ms after fb sound
%                 WaitSecs(.8);
%                 WinTabMex(3);
%                 Screen('Flip',window_ptr,[],[],2);
%                 DrawFormattedText(window_ptr,'Too slow.','center','center',[255 165 0]); 
%                 Screen('Flip',window_ptr, [], 1);
                PsychPortAudio('Start', pa_tooslow, 1, 0, 1); % "too slow"
                play_too_slow = false;
                %pause(4) % time out
                WinTabMex(2);
            end
            
            if fb_hold-fb_start > c.feedback_time
                e_feedback = 0;
                e_data = 1;
                Screen('Flip',window_ptr,[],[],2);%dont wait for synch
                find_circle = [c.start_x-dist_last,c.start_y-dist_last,c.start_x+dist_last,c.start_y+dist_last];
                Screen('FrameOval',window_ptr, white, find_circle);
                Screen('Flip',window_ptr,[],[],2);%dont wait for synch
                [~, t_FB_offset, ~, ~, ~]=Screen('Flip',window_ptr,[],1);% wait for sync, keep on screen
                t_FB_offset=t_FB_offset-start_time;
            end
        end
        
        
        %% COLLECT DATA AND PRINT ANY INSTRUCTIONS
        if e_data == 1
            
            % re-collect experiment time
            experiment_time(trialnum,count) = GetSecs-start_time;
            event_idx(trialnum,count) = 8;           
            %% STORE FAST DATA
            % basic trial data
            data.cursor_loc_final(trialnum,:) = [fb_x fb_y];
            data.hand_loc_final(trialnum,:) = [hand_ep_x hand_ep_y];
            data.aim_raw(trialnum,:) = [aim_x aim_y];
            data.target_loc(trialnum,:) = target_loc(trialnum,:);
            data.rt(trialnum) = RT;
            data.movement_time(trialnum) = movement_time;
            data.too_slow(trialnum) = e_too_slow;
            data.target_hit(trialnum) = hit;
            data.rotation(trialnum) = targetfile.rotation(trialnum);
            data.instruction(trialnum) = targetfile.instruction(trialnum);
            data.aim_report(trialnum) = targetfile.aim_report(trialnum);
            data.aim_ring(trialnum) = targetfile.aim_ring(trialnum);
            data.online_feedback(trialnum) = targetfile.online_fb(trialnum);
            data.endpoint_feedback(trialnum) = targetfile.endpoint_fb(trialnum);
            data.target_angle(trialnum) = targetfile.target_angle(trialnum);
            %store timers
            timers.hold_start(trialnum) = t_hold_start;
            timers.hold_end(trialnum) = t_hold_end;
            timers.aim_onset(trialnum) = t_aim_onset;
            timers.aim_report(trialnum) = t_aim_report;
            timers.go(trialnum) = t_go;
            timers.move_start(trialnum) = t_move_start;
            timers.move_end(trialnum) = t_move_end;
            timers.FB_onset(trialnum) = t_FB_onset;
            timers.FB_offset(trialnum) = t_FB_offset;
            %% Check for wait screen and show instructions
            if targetfile.instruction(trialnum) == 1 % use for short text
                WinTabMex(3); % Stop/Pause data acquisition.
                DrawFormattedText(window_ptr,PreAimText,'center','center',white); 
                Screen('Flip',window_ptr, [], 1);
                PsychPortAudio('Start', pa_start, 1, 0, 1); % wake up experimenter!
                KbWait(-1);
                WinTabMex(2); % Start data acquisition/clear queue.
            end
            if targetfile.instruction(trialnum) == 2 % use for longer text
                WinTabMex(3); % Stop/Pause data acquisition.
                DrawFormattedText(window_ptr,AimText,'center','center',white); 
                Screen('Flip',window_ptr, [], 1);
                PsychPortAudio('Start', pa_start, 1, 0, 1); % wake up experimenter!
                KbWait(-1);
                WinTabMex(2); % Start data acquisition/clear queue.
            end
            if targetfile.instruction(trialnum) == 3 % use for longer text
                WinTabMex(3); % Stop/Pause data acquisition.
                DrawFormattedText(window_ptr,WashText,'center','center',white); 
                Screen('Flip',window_ptr, [], 1);
                PsychPortAudio('Start', pa_start, 1, 0, 1); % wake up experimenter!
                KbWait(-1);
                WinTabMex(2); % Start data acquisition/clear queue.
            end
            
            % reset trial variables
            fb_x = 0;
            fb_y = 0;
            hand_ep_x = 0;
            hand_ep_y = 0;
            e_feedback = 0;
            e_data = 0;
            e_too_slow = 0;
            aim_x = 0;
            aim_y = 0;
            trialtoolong=false;
            refindstart=false;
            trialnum = trialnum + 1;
            count = 0;
            handcount = 0;
            e_find_start = 1;
            
        end
        
        % make sure large matrices aren't re-computed past column limit (will cause error)
        % this means trial was way too long anyway, so should likely be discarded
        if count > large_array_size && trialtoolong==false
            msg=strcat('WARNING: Trial',{' '},num2str(trialnum),' too long, losing data');
            disp(msg);
            count = large_array_size;
            trialtoolong=true;
        elseif count > large_array_size
            count = large_array_size;
        end
        
    end
    
    %% TASK COMPLETED! %%
    
    WaitSecs(1);
    
    % close sound port (important for stability of Matlab)
    PsychPortAudio('Stop', pamaster);
    PsychPortAudio('Stop', pa_tooslow);
    PsychPortAudio('Stop', pa_miss);
    PsychPortAudio('Stop', pa_hit);
    PsychPortAudio('Stop', pa_neutral);
    PsychPortAudio('Stop', pa_start);
    % Close audio device, shutdown driver:
    PsychPortAudio('Close');

    % end messege %
    Screen('DrawText', window_ptr, 'Great Job!', c.start_x-150,c.start_y-150,green);
    Screen('DrawText', window_ptr, '[saving...]', c.start_x-150,c.start_y-20,blue);
    Screen('Flip',window_ptr, [], 1);
    
    %Check for Data directory and create if missing
    try
        cd Data
    catch
        mkdir('Data');
        cd Data
    end   
     %% Put all data into one struct %%
    subdata.sub_exp_date=sub_exp_date;
    subdata.targetfilename = target_file_name;
    subdata.trialdata=data;
    subdata.timers=timers;
    subdata.constants = c;
    subdata.targetfile = targetfile;
    subdata.cursor = struct('event_idx',event_idx,'fb_hand_x',fb_hand_x,'fb_hand_y',fb_hand_y,'cursor_x_save',cursor_x_save,'cursor_y_save',cursor_y_save,'experiment_time',experiment_time);
    subdata.tablet = struct('hand_event_idx',hand_event_idx,'hand_x',hand_x,'hand_y',hand_y,'hand_time',hand_time,'hand_time_tab',hand_time_tab);
    subdata.falsestarts = falsestarts;
    subsavename = strcat(sub_exp_date,'_',target_file_name,'_',expcode);
    %last thing is to stop diary and copy to datafile
    diary off
    finallogfile = strcat('../',templogfile);
    fid=fopen(finallogfile);
    subdata.logfile=textscan(fid,'%s');
    fclose(fid);
    delete(finallogfile);
    save(subsavename,'subdata','-v7.3')
    %% Exit %%
    cd ..
    exittime = GetSecs;
    KbWait(-1,[],exittime+120); %might need to specify device/kb
    ListenChar(0); % re-activate matlab listening to keyboard (if not working, try CTRL-C)
    sca
    Priority(0);
    WinTabMex(3); % Stop/Pause data acquisition.
    WinTabMex(1);

catch err
    % close sound port (important for stability of Matlab)
    PsychPortAudio('Stop', pamaster);
    PsychPortAudio('Stop', pa_tooslow);
    PsychPortAudio('Stop', pa_miss);
    PsychPortAudio('Stop', pa_hit);
    PsychPortAudio('Stop', pa_neutral);
    PsychPortAudio('Stop', pa_start);
    % Close audio device, shutdown driver:
    PsychPortAudio('Close');
    %% NOTE: can also put a save here if you want to save in the event of an error
    % save basic data into matlab file (small file)
      save(sub_exp_date,'data'); 
    % SAVE BIG MATRICES %
    cursordata = strcat(sub_exp_date,'CURSOR');
    save(cursordata,'event_idx','fb_hand_x','fb_hand_y','cursor_x_save','cursor_y_save','experiment_time');
    tabsave = strcat(sub_exp_date,'TAB');
    save(tabsave,'hand_event_idx','hand_x','hand_y','hand_time');
    WinTabMex(3); % Stop/Pause data acquisition.
    WinTabMex(1);
    ShowCursor;
    disp(err.message);
    disp(err.stack);
    diary off
    movefile(templogfile,strcat(sub_exp_date,'error.txt'));
    sca
    Priority(0);
    ListenChar(0); % re-activate matlab listening to keyboard (if not working, try CTRL-C)
end