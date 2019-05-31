function DtSweeps()
    modelTypes = {'HRC Serbe','HRC rect Behnia T4','HRC rect Behnia T5'};
    
    % simulation length in milliseconds
    T = 4e3;
    
    % spatial extent of stimulus in degrees
    X = 35;
    
    % number of times to run the simulation
    repeats = 80;
    
    % how fast each frame updates (Hz)
    updateRate = 60;
    
    % size of pixels in space (degrees)
    pixelWidth = 5;
    
    % determines which spatial filters to use in make_ternary_stim
    % 1 for gaussian, 0 for delta function
    spatialFilterType = 1;

    % number of frames between correlations
    % we will initially simulate 13 different delays
    dts = 0:12;
    % Stimuli tensor has dimensions (time,leftVsRightInput,dts,parity)
    % number of dts is length(dts)+1 because the last dt will be
    % uncorrelated stimuli, or infinite dt
    stimuli = zeros(T*repeats,2,length(dts)+1,2);
    
    % generate a stimulus for every dt
    for ii = 1:length(dts)
        dt = dts(ii);
        parity = 1; % 1 is for phi
        stimuli(:,:,ii,1) = makeTernaryStim(T*repeats,X,updateRate,pixelWidth,dt,parity,spatialFilterType);
        parity = -1; % -1 is for reverse phi
        stimuli(:,:,ii,2) = makeTernaryStim(T*repeats,X,updateRate,pixelWidth,dt,parity,spatialFilterType);
    end
    
    % Last dt is infinite (uncorrelated stimuli)
    uncDt = -1;
    parity = 1;
    stimuli(:,:,end,1) = makeTernaryStim(T*repeats,X,updateRate,pixelWidth,uncDt,parity,spatialFilterType);
    % "negative parity" uncorrelated is the same as uncorrelated
    stimuli(:,:,end,2) = stimuli(:,:,end,1);
    
    % loop through each model type
    for ii = 1:length(modelTypes)
        thisModelType = modelTypes{ii};
        
        % turn spaces into underscores to use them as function handels
        thisModelType(thisModelType == ' ') = '_';
        modelHandle = str2func(thisModelType);
        
        % perform simulation
        output = getModelResponse(stimuli,repeats,modelHandle);
        % output has dimensions [dt,parity,direction]
        
        % generate plots of the dt sweeps
        figure();
        X = (1000/60)*dts';
        for parity = [1,2]
            subplot(1,2,parity);
            
            % Normalize by the uncorrelated response
            plot(X,squeeze(output(1:end-1,parity,:)/output(end,1,1)));
            ylim([0 2]);
            xlabel('interval (ms)');
            ylabel('response (a.u.)');
            if parity == 1
                title([modelTypes{ii} ' Positive Correlations']);
            else
                title([modelTypes{ii} ' Negative Correlations']);
            end
            hold on;
            h = plot([0 (1000/60)*dts(end)],[1 1]);
            set(h,'Color','black');
            legend({'PD','ND','unc'})
        end
    end
end

% s1 and s2 are the stimuli for the two seperate inputs
function out = HRC_Serbe(s1,s2)

    % change to luminance, not contrast
    s1 = (s1+1)/2;
    s2 = (s2+1)/2;

    % a butterworth filter is a first order lowpass, whose time constant is
    % defined as (1/(2*pi*tau))/(sampleRate/2) = 1/(pi*tau*sampleRate)
    % tau is in milliseconds and sample rate is in kHz. Sample rate here is
    % 1 kHz
    % Tm1 LP
    [BL1,AL1] = butter(1,1/pi/230,'low');
    % Tm2 LP
    [BL2,AL2] = butter(1,1/pi/100,'low');
    % Tm4 LP
    [BL4,AL4] = butter(1,1/pi/200,'low');
    % Tm9 LP
    [BL9,AL9] = butter(1,1/pi/630,'low');

    % Tm1 HP
    [BH1,AH1] = butter(1,1/pi/1230,'high');
    % Tm2 HP
    [BH2,AH2] = butter(1,1/pi/360 ,'high');
    % Tm4 HP
    [BH4,AH4] = butter(1,1/pi/250 ,'high');
    % Tm9 HP
    %[BH9,AH9] = butter(1,1/pi/inf ,'high');

    tm1 = @(x) filter(BL1,AL1,subplus(-filter(BH1,AH1,x)));
    tm2 = @(x) filter(BL2,AL2,subplus(-filter(BH2,AH2,x)));
    tm4 = @(x) filter(BL4,AL4,subplus(-filter(BH4,AH4,x)));
    tm9 = @(x) filter(BL9,AL9,subplus(               -x ));


    out = ((tm1(s1).*tm2(s2))...
          +(tm1(s1).*tm4(s2))...
          +(tm9(s1).*tm1(s2))...
          +(tm4(s1).*tm2(s2))...
          +(tm9(s1).*tm2(s2))...
          +(tm9(s1).*tm4(s2)))/6;
end


function out = HRC_rect_Behnia_T4(s1,s2)

    q = load('data/BehniaData.mat');

    f1 = q.Mi1;
    f2 = q.Tm3;

    f1 = f1./sum(f1.^2);
    f2 = f2./sum(f2.^2);

    out = subplus(filter(f1,1,s1)).*subplus(filter(f2,1,s2));
end

function out = HRC_rect_Behnia_T5(s1,s2)

    q = load('data/BehniaData.mat');

    f1 = q.Tm1;
    f2 = q.Tm2;

    f1 = f1./sum(f1.^2);
    f2 = f2./sum(f2.^2);

    out = subplus(filter(f1,1,s1)).*subplus(filter(f2,1,s2));
end