function [] = run_FShift_HighFreqFlick(sub,flag_training, flag_isolum, flag_block)
% run_FShift_HighFreqFlick(sub,flag_training, flag_isolum, flag_block) also known as HiFli
%   runs experiment SSVEP_FShift_HighFreqFlick
%       sub:            participant number
%       flag_training:  1 = do training
%       flag_isolum:    1 = do isoluminance adjustment
%       flag_block:     1 = start with block 1
%           e.g. run_FShift_PerIrr(1,1, 0, 1)
% 
% current version includes two irrelevant colors in periphery
%
% log
%   - inlay flicker between background color and white or RDK color and white?
%   - fixed position or vary across participants?


% Christopher Gundlach, Maria Dotzer,  Leipzig, 2025,2024,2023,2021,2020
if nargin < 4
    help run_FShift_PerIrr
    return
end

%% parameters
% sub = 1; flag_training = 1; flag_block = 1; flag_isolum = 1;
% design
p.sub                   = sub;                  % subject number
p.flag_block            = flag_block;           % block number to start
p.flag_training         = flag_training;        % do training

p.ITI                   = [1000 1000];          % inter trial interval in ms
p.targ_respwin          = [200 1000];           % time window for responses in ms

% screen
p.scr_num               = 1;                    % screen number
p.scr_res               = [1920 1080];          % resolution
p.scr_refrate           = 480;                  % refresh rate in Hz (e.g. 85)
p.scr_color             = [0.2 0.2 0.2 1];      % default: [0.05 0.05 0.05 1]; ; color of screen [R G B Alpha]
p.scr_imgmultipl        = 4;

% some isoluminace parameters
p.isol.TrlAdj           = 5;                    % number of trials used for isoluminance adjustment
p.isol.MaxStd           = 10;                   % standard deviation tolerated
p.isol.run              = false;                % isoluminance run?
% p.isol.override         = [0.0980 0.0392 0 1; 0 0.0596 0.1490 1;0 0.0745 0 1]; % these are the ones used for p.isol.bckgr = p.scr_color(1:3);
p.isol.override         = [0.3529 0.1412 0 1; 0 0.2745 0.6863 1]; % these are the ones used for p.isol.bckgr = p.scr_color(1:3);

% p.isol.bckgr            = p.scr_color(1:3)+0.2;          % isoluminant to background or different color?
p.isol.bckgr            = p.scr_color;          % isoluminant to background or different color?


% stimplan
p.stim.condition        = [1 2 3 4];
p.stim.RDKcenter        = [1 2; 1 2; 5 6; 5 6]; % defines which RDK are shown in center [1 2 = convent. RDKs; 5 6 = inlet flicker]
p.stim.RDKperi          =  [3 4; 3 4; ...   % defines which RDK are shown in the periphery [3 4 = convent. RDKs; 7 8 = inlet flicker]
                            7 8; 7 8];
p.stim.RDK2attend       = [1 2 5 6];        % defines which RDK to attend in which condition
p.stim.eventnum_e       = [0 0 0 0 1 2];    % ratio of eventnumbers for experiment
p.stim.eventnum_t       = [0 1 2];          % ratio of eventnumbers for training
p.stim.con_repeats      = [25 25 25 25];    % trial number/repeats for each eventnum and condition
p.stim.con_repeats_t    = [1];              % trial number/repeats for each eventnum and condition
p.stim.trialnum_t       = 20;               % trial number in training
p.stim.time_postcue     = 2;                % post.cue time in s
p.stim.time_precue      = [1.5 2];          % precue time in s; [upper lower] for randomization
p.stim.event.type       = 2;                % types of events (1 = targets only, 2 = targets + distrators)
p.stim.event.length     = 0.3;              % lengt of events in s
p.stim.event.min_onset  = 0.2;              % min post-cue time before event onset in s
p.stim.event.min_offset = 0;                % min offset from target end to end of trial in s
p.stim.event.min_dist   = 0.8;              % min time between events in s
p.stim.blocknum         = 20;               % number of blocks
p.stim.ITI              = [1 1];            % ITI range in seconds
p.stim.frames_postcue   = p.stim.time_postcue*p.scr_refrate;


% introduce RDK structure
RDK.RDK(1).id           = 1;                            % id of RDK
RDK.RDK(1).size         = [154 308];                    % width and height of RDK in pixel; only even values [38 = 9.6°]
RDK.RDK(1).centershift  = [0 0];                        % position of RDK center; x and y deviation from center in pixel
RDK.RDK(1).col          = [1 0.4 0 1; p.scr_color(1:3) 0];% "on" and "off" color
RDK.RDK(1).freq         = 0;                            % flicker frequency, frequency of a full "on"-"off"-cycle
RDK.RDK(1).mov_freq     = 120;                          % Defines how frequently the dot position is updated; 0 will adjust the update-frequency to your flicker frequency (i.e. dot position will be updated with every "on"-and every "off"-frame); 120 will update the position for every frame for 120Hz or for every 1. quadrant for 480Hz 
RDK.RDK(1).num          = 85;                           % number of dots % 85
RDK.RDK(1).mov_speed    = 1;                            % movement speed in pixel
RDK.RDK(1).mov_dir      = [0 1; 0 -1; -1 0; 1 0];       % movement direction  [0 1; 0 -1; -1 0; 1 0] = up, down, left, right
RDK.RDK(1).dot_size     = 10;                           % size of dots
RDK.RDK(1).shape        = 1;                            % 1 = square RDK; 0 = ellipse/circle RDK;
RDK.RDK(1).in_flag      = 0;                            % display inlet?
RDK.RDK(1).in_col       = [];                           % inlet color for on off
RDK.RDK(1).in_freq      = 0;                            % flicker frequency, frequency of a full "on"-"off"-cycle
RDK.RDK(1).in_dot_size  = 6;                            % size of dots
RDK.RDK(1).flicker_type = 'SSVEP';                      % have only inlet or SSVEP flicker


p.stim.pos_shift        = [-255 0; 255 0];              % position shift in pixel for stimuli in periphery [255 = 7.8°] either left or right
% p.stim.freqs            = {[26 29];[17 20 23]};         % frequencies of {[center1 center2];[peri1 peri2 peri3]}
p.stim.freqs            = {[20 23];[14 17]};            % conventional frequencies of {[center1 center2];[peri1 peri2]}
p.stim.in_freqs         = {[68 71];[62 65]};            % high inlet frequencies of {[center1 center2];[peri1 peri2]}
p.stim.colors           = ...                           % "on" and "off" color
    {[1 0.4 0 1; p.scr_color(1:3) 0];...
    [0 0.4 1 1; p.scr_color(1:3) 0]};
    % plot_colorwheel([1 0.4 0; 0 0.4 1; 0 1 0; 1 0 1],'ColorSpace','propixxrgb','LAB_L',50,'NumSegments',60,'AlphaColWheel',1,'LumBackground',100)
    % plot_colorwheel([1 0.4 0; 0 0.4 1],'ColorSpace','propixxrgb','LAB_L',50,'NumSegments',60,'AlphaColWheel',1,'LumBackground',100)
% p.stim.in_colors           = ...                           % "on" and "off" color for inlet (from background to to white)
%     [1 1 1 1; p.scr_color(1:3) 1];
p.stim.in_colors           = ...                           % "on" and "off" color for inlet (from background to to white) same brightnes as background on average
    [0.4 0.4 0.4 1; 0 0 0  1];
p.stim.color_names      = {'redish';'blue'};
 
RDK.event.type          = 'globalmotion';       % event type global motion
RDK.event.duration      = p.stim.event.length;  % time of coherent motion
RDK.event.coherence     = .4;                   % percentage of coherently moving dots 0.4 [changed from 0.4 to 0.3 to 0.4]
RDK.event.direction     = RDK.RDK(1).mov_dir;   % movement directions for events

% fixation cross
p.crs.color             = [0.4 0.4 0.4 1];      % color of fixation cross
p.crs.size              = 12;                   % size of fixation
p.crs.width             = 2;                    % width of fixation cross
p.crs.cutout            = 0;                    % 1 = no dots close to fixation cross

% trigger
p.trig.rec_start        = 253;                  % trigger to start recording
p.trig.rec_stop         = 254;                  % trigger to stop recording
p.trig.tr_start         = 77;                   % trial start; main experiment
p.trig.tr_stop          = 88;                   % trial end; main experiment
p.trig.tr_con_type      = [1 2 3 4]*10;         % condition type
p.trig.type             = [1 2; 5 7];           % [first: target, distractor; second: target, distractor]
p.trig.button           = 150;                   % button press
p.trig.event_type       = [201 202];              % target, distractor

% possible condition triggers:
% {[1 101 201 111 121 211 221]; [2 102 202 112 122 212 222]; [3 103 203 113 123 213 223]; ...
% [4 104 204 114 124 214 224]; [5 105 205 115 125 215 225]; [6 106 206 116 126 216 226]}

% logfiles
p.log.path              = '/home/stimulation120/matlab/user/christopher/stim_FShift_HighFli/logfiles/';
p.log.exp_name          = 'SSVEP_FShift_HighFreqFlickr';
p.log.add               = '_a';


%% check for logfile being present
filecheck=dir(sprintf('%sVP%02.0f_timing*',p.log.path,p.sub));
if ~isempty(filecheck)
    reply = input(sprintf('\nVP%02.0f existiert bereits. Datei überschreiben? [j/n]... ',p.sub),'s');
    if strcmp(reply,'j')
        p.filename = sprintf('VP%02.0f_timing',p.sub);
    else
        [temp name_ind]=max(cellfun(@(x) numel(x), {filecheck.name}));
        p.filename = sprintf('%s%s',filecheck(name_ind).name(1:end-4),p.log.add);
    end
else
    p.filename = sprintf('VP%02.0f_timing',p.sub);
end

t.isol = {};
% routine to check for older isoluminance adjustments
for i_file = 1:numel(filecheck)
    t.in = load(fullfile(filecheck(i_file).folder,filecheck(i_file).name));
    t.datenum{i_file} = filecheck(i_file).datenum;
    t.isol{i_file} = t.in.p.isol;
    
end



%% Screen init
ps.input = struct('ScrNum',p.scr_num,'RefRate',p.scr_refrate,'PRPXres',p.scr_res,'BckGrCol',p.scr_color,'PRPXmode',2);
[~, ps.screensize, ps.xCenter, ps.yCenter, ps.window, ps.framerate, ps.RespDev, ps.keymap] = PTExpInit_GLSL(ps.input,1);

% some initial calculations
% for debugging: ps.xCenter = round(p.scr_res(1)/2); ps.yCenter = round(p.scr_res(2)/2)
% fixation cross
ps.center = [ps.xCenter ps.yCenter];
p.crs.half = p.crs.size/2;
p.crs.bars = [-p.crs.half p.crs.half 0 0; 0 0 -p.crs.half p.crs.half];

% shift into 4 quadrants (running with 480 Hz)
ps.shift = [-ps.xCenter/2, -ps.yCenter/2; ps.xCenter/2, -ps.yCenter/2;... % shifts to four quadrants: upper left, upper right, lower left, lower right
    -ps.xCenter/2, ps.yCenter/2; ps.xCenter/2, ps.yCenter/2];

p.crs.lines = [];
for i_quad=1:p.scr_imgmultipl
    p.crs.lines = cat(2, p.crs.lines, [p.crs.bars(1,:)+ps.shift(i_quad,1); p.crs.bars(2,:)+ps.shift(i_quad,2)]); %array with start and end points for the fixation cross lines, for all four quadrants
end

%% keyboard and ports setup ???
KbName('UnifyKeyNames')
Buttons = [KbName('ESCAPE') KbName('Q') KbName('SPACE') KbName('j') KbName('n') KbName('1!') KbName('2@') KbName('3#')];
RestrictKeysForKbCheck(Buttons);
key.keymap=false(1,256);
key.keymap(Buttons) = true;
key.keymap_ind = find(key.keymap);
[key.ESC, key.SECRET, key.SPACE, key.YES, key.NO] = deal(...
    Buttons(1),Buttons(2),Buttons(3),Buttons(4),Buttons(5));

%% start experiment
% initialize randomization of stimulation frequencies and RDK colors
% inititalize RDKs [RDK1 and RDK2 task relevant at center;  RDK3 RDK4 RDK5 not and in periphery]
rng(1,'v4')
% quasi randomize position
t.pos = [];
for i_rep = 1:100
    t.pos = cat(1,t.pos,p.stim.pos_shift(randsample(1:2,2),:));
end
% keep position fixed across participants?
t.pos = repmat(p.stim.pos_shift(1,:),100,1);


% rand('state',p.sub)
rng(p.sub,'v4')

% assign colors
RDK.RDK(1).col_init = RDK.RDK(1).col;
RDK.RDK(2:8) = deal(RDK.RDK(1));
t.col = repmat(p.stim.colors,4,1);
[RDK.RDK(:).col_init] = deal(t.col{:});

% assign flicker type
[RDK.RDK(5:8).in_flag] = deal(1); % which onses are unconventional RDKs?
[RDK.RDK(5:8).flicker_type] = deal('SSVEP_inlet'); % which onses are unconventional RDKs?

% assign ID
t.ids = num2cell(1:numel(RDK.RDK));
[RDK.RDK(:).id] = t.ids{:};

% assign inlet color
[RDK.RDK(5:8).in_col] = deal(p.stim.in_colors); 


% randomize frequencies for conventional RDKs
t.val = num2cell(p.stim.freqs{1}(randperm(2)));
[RDK.RDK(1:2).freq] = t.val{:};
t.val = num2cell(p.stim.freqs{2}(randperm(2)));
[RDK.RDK(3:4).freq] = t.val{:};

% randomize frequencies for inlet flickering RDKs
t.val = num2cell(p.stim.in_freqs{1}(randperm(2)));
[RDK.RDK(5:6).in_freq] = t.val{:};
t.val = num2cell(p.stim.in_freqs{2}(randperm(2)));
[RDK.RDK(7:8).in_freq] = t.val{:};

% assign colors
[RDK.RDK(:).col] = deal(RDK.RDK(:).col_init);
t.col = repmat(p.stim.color_names,4,1);
[RDK.RDK(:).colnames] = deal(t.col{:});
p.isol.override = p.isol.override;
p.isol.init_cols = cell2mat(cellfun(@(x) x(1,:),{RDK.RDK(:).col},'UniformOutput',false)');



% % randomize colors? no
% % [RDK.RDK([1 2 5 6]).col] = deal(p.stim.colors{randperm(4)});
% % [RDK.RDK(3:4).col] = deal(RDK.RDK(1:2).col);
% p.colrandidx = randperm(3); p.colrandidx = [p.colrandidx(1:2) p.colrandidx];
% [RDK.RDK(:).col] = deal(p.stim.colors{p.colrandidx});
% [RDK.RDK(:).colnames] = deal(p.stim.color_names{p.colrandidx});
% p.isol.override = p.isol.override(p.colrandidx,:);
% p.isol.init_cols = cell2mat(cellfun(@(x) x(1,:),{RDK.RDK(:).col},'UniformOutput',false)');

% random position shift in periphery
[RDK.RDK([3:4 7:8]).centershift] = deal(t.pos(p.sub,:));

% initialize blank variables
timing = []; button_presses = []; resp = []; randmat = [];

%% initial training
if p.flag_training
    fprintf(1,'\nTraing starten mit q')
    inp.prompt_check = 0;
    while inp.prompt_check == 0             % loop to check for correct input
        [key.keyisdown,key.secs,key.keycode] = KbCheck;
        if key.keycode(key.SECRET)==1
            flag_trainend = 0; inp.prompt_check = 1;
        end
        Screen('Flip', ps.window, 0);
    end
    
    
    i_bl = 1;
    flag_trainend = 0;
    while flag_trainend == 0 % do training until ended
        %rand('state',p.sub*i_bl) % determine randstate
        rng(p.sub*i_bl,'v4')
        randmat.training{i_bl} = rand_FShift_HighFreqFlick(p, RDK,  1);
        [timing.training{i_bl},button_presses.training{i_bl},resp.training{i_bl}] = ...
            pres_FShift_HighFreqFlick(p, ps, key, RDK, randmat.training{i_bl}, i_bl,1);
        save(sprintf('%s%s',p.log.path,p.filename),'timing','button_presses','resp','randmat','p', 'RDK')
        pres_feedback(resp.training{i_bl},p,ps, key,RDK)
               
        % loop for training to be repeated
        fprintf(1,'\nTraing wiederholen? (j/n)')
        inp.prompt_check = 0;
        while inp.prompt_check == 0             % loop to check for correct input
            [key.keyisdown,key.secs,key.keycode] = KbCheck; 
            if key.keycode(key.YES)==1
                i_bl = i_bl + 1; flag_trainend = 0; inp.prompt_check = 1;
            elseif key.keycode(key.NO)==1
                flag_trainend = 1; inp.prompt_check = 1;
            end
            Screen('Flip', ps.window, 0);
        end  
        
    end
end

%% then isoluminance adjustment
% do the heterochromatic flicker photometry
ttt=WaitSecs(0.7);
if flag_isolum == 1
%     
%     PsychDefaultSetup(2);
%     Datapixx('Open');
%     Datapixx('SetPropixxDlpSequenceProgram', 0);
%     Datapixx('RegWrRd');
     
    
    
    % % start isoluminance script only RGB output (no alpha) [for all colors]
    % [Col2Use] = PRPX_IsolCol_480_adj(...
    %     [p.isol.bckgr(1:3); p.isol.init_cols(:,1:3)],...
    %     p.isol.TrlAdj,...
    %     p.isol.MaxStd,...
    %     cellfun(@(x) x(1), {RDK.RDK.centershift})',...
    %     RDK.RDK(1).size);
    % start isoluminance script only RGB output (no alpha) [for only a few colors]
    [Col2Use] = PRPX_IsolCol_480_adj(...
        [p.isol.bckgr(1:3); p.isol.init_cols(1:2,1:3)],...
        p.isol.TrlAdj,...
        p.isol.MaxStd,...
        cellfun(@(x) x(1), {RDK.RDK(1:2).centershift})',...
        [1 1].*max(RDK.RDK(1).size));
    
    for i_RDK = 1:numel(RDK.RDK)
        RDK.RDK(i_RDK).col(1,:) = [Col2Use(2+mod(i_RDK+1,2),:) 1];
    end
    % index function execution
    p.isol.run = sprintf('originally run: %s',datestr(now));
    p.isol.coladj = [Col2Use(2:end,:) ones(size(Col2Use,1)-1,1)];
    save(sprintf('%s%s',p.log.path,p.filename),'timing','button_presses','resp','randmat','p', 'RDK')
    
    fprintf('\nadjusted colors:\n')
    for i_col = 1:size(p.isol.coladj,1)
        fprintf('RDK%1.0f [%1.4f %1.4f %1.4f %1.4f]\n', i_col,p.isol.coladj(i_col,:))
    end
    
    Screen('CloseAll')
    Datapixx('SetPropixxDlpSequenceProgram', 0);
    Datapixx('RegWrRd');
    Datapixx('close');
else
    % select colors differently
    fprintf(1,'\nKeine Isoluminanzeinstellung. Wie soll verfahren werden?')
    % specify options
    % option1: use default values
    isol.opt(1).available = true;
    t.cols = cell2mat({RDK.RDK(1:2).col}');
    isol.opt(1).colors = t.cols(1:2:end,:);
    isol.opt(1).text = sprintf('default: %s',sprintf('[%1.2f %1.2f %1.2f] ',isol.opt(1).colors(:,1:3)'));
    % option2: use isoluminance values of previously saved dataset
    if ~isempty(t.isol) % file loaded 
        [t.t t.idx] = max(cell2mat(t.datenum));
        if any(strcmp(fieldnames(t.isol{t.idx}),'coladj')) % and adjusted colors exist?s
            isol.opt(2).available = true;
            isol.opt(2).colors = t.isol{t.idx}.coladj(1:end,:);
            isol.opt(2).text = sprintf('aus gespeicherter Datei: %s',sprintf('[%1.2f %1.2f %1.2f] ',isol.opt(2).colors(:,1:3)'));
        else
            isol.opt(2).available = false;
            isol.opt(2).colors = [];
            isol.opt(2).text = [];
        end
    else
        isol.opt(2).available = false;
        isol.opt(2).colors = [];
        isol.opt(2).text = [];
    end
    % option3: use manual override
    if ~isempty(p.isol.override)
        isol.opt(3).available = true;
        isol.opt(3).colors = p.isol.override;
        isol.opt(3).text = sprintf('manuell definiert in p.isol override: %s',sprintf('[%1.2f %1.2f %1.2f] ',isol.opt(3).colors(:,1:3)'));
    else
        isol.opt(3).available = false;
        isol.opt(3).colors = [];
        isol.opt(3).text = [];
    end
    % check for buttons
    IsoButtons = Buttons(6:8);
    isol.prompt.idx = find([isol.opt(:).available]);
    t.prompt = [];
    for i_prompt = 1:numel(isol.prompt.idx)
        t.prompt = [t.prompt sprintf('\n(%1.0f) %s',i_prompt,isol.opt(isol.prompt.idx(i_prompt)).text)];
    end
    
    % display options
    fprintf('%s',t.prompt)
    inp.prompt_check = 0;
    while inp.prompt_check == 0             % loop to check for correct input
        [key.keyisdown,key.secs,key.keycode] = KbCheck;
        if any(key.keycode)
            inp.prompt_check = 1;
        end
        Screen('Flip', ps.window, 0);
    end
    Col2Use = isol.opt(isol.prompt.idx(key.keycode(IsoButtons(1:numel(isol.prompt.idx)))==1)).colors;
    % use selected colors
    for i_RDK = 1:numel(RDK.RDK)
        RDK.RDK(i_RDK).col(1,:) = Col2Use(1+mod(i_RDK+1,2),:);
    end
    % index function execution
    switch isol.prompt.idx(key.keycode(IsoButtons(1:numel(isol.prompt.idx)))==1)
        case 1
            p.isol.run = sprintf('default at %s',datestr(now));
        case 2
            p.isol.run = sprintf('reloaded at %s from %s',datestr(now),datestr(t.datenum{t.idx}));
        case 3
            p.isol.run = sprintf('override at %s',datestr(now));
    end
    p.isol.coladj = Col2Use;
%     save(sprintf('%s%s',p.log.path,p.filename),'timing','button_presses','resp','randmat','p', 'RDK')
    
    fprintf('\nselected colors:\n')
    for i_col = 1:size(p.isol.coladj,1)
        fprintf('RDK%1.0f [%1.4f %1.4f %1.4f %1.4f]\n', i_col,p.isol.coladj(i_col,:))
    end
end

%% redo initialization
ps.input = struct('ScrNum',p.scr_num,'RefRate',p.scr_refrate,'PRPXres',p.scr_res,'BckGrCol',p.scr_color,'PRPXmode',2);
[~, ps.screensize, ps.xCenter, ps.yCenter, ps.window, ps.framerate, ps.RespDev, ps.keymap] = PTExpInit_GLSL(ps.input,1);

% some initial calculations
% fixation cross
ps.center = [ps.xCenter ps.yCenter];
p.crs.half = p.crs.size/2;
p.crs.bars = [-p.crs.half p.crs.half 0 0; 0 0 -p.crs.half p.crs.half];

% shift into 4 quadrants (running with 480 Hz)
ps.shift = [-ps.xCenter/2, -ps.yCenter/2; ps.xCenter/2, -ps.yCenter/2;... % shifts to four quadrants: upper left, upper right, lower left, lower right
    -ps.xCenter/2, ps.yCenter/2; ps.xCenter/2, ps.yCenter/2];

p.crs.lines = [];
for i_quad=1:p.scr_imgmultipl
    p.crs.lines = cat(2, p.crs.lines, [p.crs.bars(1,:)+ps.shift(i_quad,1); p.crs.bars(2,:)+ps.shift(i_quad,2)]); %array with start and end points for the fixation cross lines, for all four quadrants
end

% keyboard setup
KbName('UnifyKeyNames')
Buttons = [KbName('ESCAPE') KbName('Q') KbName('SPACE') KbName('j') KbName('n') KbName('1!') KbName('2@') KbName('3#')];
RestrictKeysForKbCheck(Buttons);
key.keymap=false(1,256);
key.keymap(Buttons) = true;
key.keymap_ind = find(key.keymap);
[key.ESC, key.SECRET, key.SPACE, key.YES, key.NO] = deal(...
    Buttons(1),Buttons(2),Buttons(3),Buttons(4),Buttons(5));


%% do training again?
% loop for training to be repeated
fprintf(1,'\nTraing starten (j/n)')
inp.prompt_check = 0;
while inp.prompt_check == 0             % loop to check for correct input
    [key.keyisdown,key.secs,key.keycode] = KbCheck;
    if key.keycode(key.YES)==1
        flag_trainend = 0; inp.prompt_check = 1;
    elseif key.keycode(key.NO)==1
        flag_trainend = 1; inp.prompt_check = 1;
    end
    Screen('Flip', ps.window, 0);
end

if ~exist('i_bl'); i_bl = 1; end
while flag_trainend == 0 % do training until ended
    %rand('state',p.sub*i_bl) % determine randstate
    rng(p.sub*i_bl,'v4')
    randmat.training{i_bl} = rand_FShift_HighFreqFlick(p, RDK,  1);
    [timing.training{i_bl},button_presses.training{i_bl},resp.training{i_bl}] = ...
        pres_FShift_HighFreqFlick(p, ps, key, RDK, randmat.training{i_bl}, i_bl,1);
    save(sprintf('%s%s',p.log.path,p.filename),'timing','button_presses','resp','randmat','p', 'RDK')
    pres_feedback(resp.training{i_bl},p,ps, key,RDK)
    
    % loop for training to be repeated
    fprintf(1,'\nTraing wiederholen? (j/n)')
    inp.prompt_check = 0;
    while inp.prompt_check == 0             % loop to check for correct input
        [key.keyisdown,key.secs,key.keycode] = KbCheck;
        if key.keycode(key.YES)==1
            i_bl = i_bl + 1; flag_trainend = 0; inp.prompt_check = 1;
        elseif key.keycode(key.NO)==1
            flag_trainend = 1; inp.prompt_check = 1;
        end
        Screen('Flip', ps.window, 0);
    end
    
end


%% present each block
% randomization
% rand('state',p.sub);                         % determine randstate
rng(p.sub,'v4')
randmat.experiment = rand_FShift_HighFreqFlick(p, RDK,  0);    % randomization
for i_bl = p.flag_block:p.stim.blocknum
    % start experiment
    [timing.experiment{i_bl},button_presses.experiment{i_bl},resp.experiment{i_bl}] = ...
        pres_FShift_HighFreqFlick(p, ps, key, RDK, randmat.experiment, i_bl,0);
    % save logfiles
    save(sprintf('%s%s',p.log.path,p.filename),'timing','button_presses','resp','randmat','p', 'RDK')
          
    pres_feedback(resp.experiment{i_bl},p,ps, key, RDK)    
end

fprintf(1,'\n\nENDE\n')

%Close everything
Datapixx('SetPropixxDlpSequenceProgram', 0);
Datapixx('RegWrRd');
Datapixx('close');
ppdev_mex('Close', 1);
ListenChar(0);
sca;


end

