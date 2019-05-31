function ToyModels()
    modelTypes = {'HRC','HRC-rect','HRC-LN','BL'};

    % simulation length in milliseconds
    T = 4e3;
    % spatial extent of stimulus in degrees
    X = 35;
    % number of times to repeat the experiment, decreases noise
    repeats = 8;
    % frequency of stimulus change (Hz)
    updateRate = 60;
    % width of stimulus pixels (degrees)
    pixelWidth = 5;

    % Temporal delay of the correlations in the stimulus (frames)
    dt = 1;

    for filterType = {'simple','realistic'}
        % Generate Stimuli
        % Blur stimuli in case of 'Realistic' filter computation
        if strcmp(filterType,'simple')
            spatialFilterType = 0;
        else
            spatialFilterType = 1;
        end

        % Stimuli tensor has dimensions (time,leftOrRightInput,dt,parity)
        stimuli = zeros(T*repeats,2,2,2);
        parity = 1;
        stimuli(:,:,1,1) = makeTernaryStim(T*repeats,X,updateRate,pixelWidth,dt,parity,spatialFilterType);
        parity = -1;
        stimuli(:,:,1,2) = makeTernaryStim(T*repeats,X,updateRate,pixelWidth,dt,parity,spatialFilterType);
        
        % a dt of -1 signified "uncorrelated" or "infinite delay"
        % we'll put the uncorrelated stimulus as the last dt
        uncDt = -1;
        parity = 1;
        stimuli(:,:,end,1) = makeTernaryStim(T*repeats,X,updateRate,pixelWidth,uncDt,parity,spatialFilterType);
        stimuli(:,:,end,2) = stimuli(:,:,end,1);

        % Compute output
        results = zeros(length(modelTypes),5);
        for ii = 1:length(modelTypes)
            modelType = modelTypes{ii};
            % replace dashes with underscores so we can turn the name into
            % a function call
            modelType(modelType == '-') = '_';
            modelHandle = str2func(modelType);
            
            % Create a model function with the temporal filters specified
            [delayFilter,nonDelayFilter] = getFilters(filterType);
            specificModelHandle = @(s1,s2) modelHandle(delayFilter,nonDelayFilter,s1,s2);
            
            % compute the results of the models
            output = getModelResponse(stimuli,repeats,specificModelHandle);
            % output format: dt,parity,direction
            % results format: model,[parity,direction;rand]
            results(ii,1:4) = output(1,:);
            results(ii,5) = output(2,1,1);
        end
        
        % opponent, subtract one direction of motion from the opposite
        % direction
        results(:,6) = results(:,1) - results(:,3);
        results(:,7) = results(:,2) - results(:,4);
        
        % Plot results
        figure();
        numModels = length(modelTypes);
        barPos = [1:5 7 8];
        for model=1:numModels;
            subplot(numModels,1,model);
            bar(barPos,results(model,:));
            set(gca,'XTickLabel',{'PD+','PD-','ND+','ND-','rand','Net+','Net-'});
            % Plot dotted line at rand
            hold on;
            h = plot([0 6],repelem(results(model,5),2));
            set(h,'LineStyle','--','Color','black')
            % Add title
            titleTxt = [modelTypes{model} ' with ' filterType{1} ' filters'];
            title(titleTxt);
        end
    end
end

function [delayFilter, nonDelayFilter] = getFilters(filterType)
    filterLength = 200;
    t = 0:filterLength;
    tau = round(1000/60);
    
    if strcmp(filterType,'realistic')
        delayFilter = t.*exp(-t/tau);
        delayFilter = delayFilter'/sum(delayFilter);
        
        nonDelayFilter = zeros(filterLength,1);
        nonDelayFilter(1) = 1;
    else
        delayFilter = zeros(filterLength,1);
        delayFilter(1+tau) = 1;
        
        nonDelayFilter = zeros(filterLength,1);
        nonDelayFilter(1) = 1;
    end
end

function out = HRC(delayFilter,nonDelayFilter,s1,s2)
    out = filter(delayFilter,1,s1).*filter(nonDelayFilter,1,s2);
end

function out = HRC_rect(delayFilter,nonDelayFilter,s1,s2)
    out = subplus(filter(delayFilter,1,s1)).*subplus(filter(nonDelayFilter,1,s2));
end

function out = HRC_LN(delayFilter,nonDelayFilter,s1,s2)
    out = subplus(filter(delayFilter,1,s1)+filter(nonDelayFilter,1,s2));
end

function out = BL(delayFilter,nonDelayFilter,s1,s2)
    a = filter(nonDelayFilter,1,s1);
    b = filter(delayFilter,1,s2);
    out = subplus(a-b);
end