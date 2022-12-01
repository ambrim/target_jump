%% generate trial file
baselineTrialTypes = [1, 1, 0, 360, 4, 0;... % endpoint_feedback, online_feedback, rotation, target_angle, size
                      1, 1, 0, 360, 24, 0];
                  
nofbBaselineTrialTypes = [0, 0, 0, 360, 4, 0;... % endpoint_feedback, online_feedback, rotation, target_angle, size
                      0, 0, 0, 360, 24, 0;...
                      0, 0, 0, 360, 4, 0;...
                      0, 0, 0, 360, 24, 0];
                  
clampTutorialTrialTypes = [1, 1, 1.75, 360, 4, 1;... % endpoint_feedback, online_feedback, rotation, target_angle, size
                      1, 1, -1.75, 360, 4, 2;...
                      1, 1, 1.75, 360, 24, 1;... % endpoint_feedback, online_feedback, rotation, target_angle, size
                      1, 1, -1.75, 360, 24, 2];
                  
testBlockTrialTypes = [1, 1, 1.75, 360, 4, 1;... % endpoint_feedback, online_feedback, rotation, target_angle, size
                      1, 1, -1.75, 360, 4, 2;...
                      1, 1, 1.75, 360, 4, 2;...
                      1, 1, -1.75, 360, 4, 1;...
                      1, 1, 1.75, 360, 24, 1;...
                      1, 1, -1.75, 360, 24, 2;...
                      1, 1, 1.75, 360, 24, 2;...
                      1, 1, -1.75, 360, 24, 1];

                  
%targetAngles = [80, 200, 320];
targetSizes = [4 ,24];                 
numBaselineBlocks = 3;
numNoFBBaselineBlocks = 5;
numTestBlocks = 32;

baseline_length = 2;
noFBbaseline_length = 4;

trialContents = nan((baseline_length*numBaselineBlocks) + (...
    noFBbaseline_length * numNoFBBaselineBlocks) + (...
    length(testBlockTrialTypes) * numTestBlocks * 3) + 4, 6);
thisTrialNum = 1;
oneback_rwd_outer = NaN;
twoback_rwd_outer = NaN;
oneback_size_outer = NaN;
twoback_size_outer = NaN;
% baseline with feedback
for i = 1:numBaselineBlocks
    shuffleGood = 0;
    while shuffleGood == 0        
        shuffleGood = 1;
        
        almostOrder = randperm(baseline_length);
        tempTrialTypes = baselineTrialTypes(almostOrder,:);
        tempTrialTypes = [tempTrialTypes,...
            zeros(baseline_length,1)];
        
        trialRows = [];
        targetOrder = [];
        for r = 1:2
            order = randperm(2);
            theseSizes = targetSizes(order); 
            for t = 1:2 %t now stands for target SIZE
                [rows, cols] = size(tempTrialTypes);
                found = 0;
                for p = 1:rows
                    if (found == 0) && (tempTrialTypes(p,7) == 0) ...
                            && (tempTrialTypes(p,5) == theseSizes(t))
                        thisRow = tempTrialTypes(p,1:6);
                        tempTrialTypes(p,7) = 1;
                        found = 1;
                    end
                end
                trialRows = [trialRows; thisRow];
            end
        end
        
        
        oneback_rwd = oneback_rwd_outer;
        twoback_rwd = twoback_rwd_outer;
        oneback_size = oneback_size_outer;
        twoback_size = twoback_size_outer;
        [rows, cols] = size(trialRows);
        for j = 1:rows
            if shuffleGood == 1
%                 if (oneback_size == trialRows(j, 5)) && ...
%                         (oneback_size == twoback_size)
%                     shuffleGood = 0;
%                 end
                twoback_rwd = oneback_rwd;
                oneback_rwd = trialRows(j, 6);
                twoback_size = oneback_size;
                oneback_size = trialRows(j, 5);
            end
        end
    end
    oneback_rwd_outer = oneback_rwd;
    twoback_rwd_outer = twoback_rwd;
    oneback_size_outer = oneback_size;
    twoback_size_outer = twoback_size;
    thisBlock = trialRows;
    
    iteration = 1;
    for t = thisTrialNum:(thisTrialNum + baseline_length - 1)
        trialContents(t,:) = thisBlock(iteration, :);
        iteration = iteration + 1;
    end
    thisTrialNum = (thisTrialNum + baseline_length);
end
lastBLTrial = thisTrialNum - 1;

% baseline without feedback
for i = 1:numNoFBBaselineBlocks
    shuffleGood = 0;
    while shuffleGood == 0        
        shuffleGood = 1;
        
        almostOrder = randperm(noFBbaseline_length);
        tempTrialTypes = nofbBaselineTrialTypes(almostOrder,:);
        tempTrialTypes = [tempTrialTypes,...
            zeros(noFBbaseline_length,1)];
        
        trialRows = [];
        targetOrder = [];
        for r = 1:2
            order = randperm(2);
            theseSizes = targetSizes(order);
            for t = 1:2 %t now stands for target SIZE
                [rows, cols] = size(tempTrialTypes);
                found = 0;
                for p = 1:rows
                    if (found == 0) && (tempTrialTypes(p,7) == 0) ...
                            && (tempTrialTypes(p,5) == theseSizes(t))
                        thisRow = tempTrialTypes(p,1:6);
                        tempTrialTypes(p,7) = 1;
                        found = 1;
                    end
                end
                trialRows = [trialRows; thisRow];
            end
        end
        
        
        oneback_rwd = oneback_rwd_outer;
        twoback_rwd = twoback_rwd_outer;
        oneback_size = oneback_size_outer;
        twoback_size = twoback_size_outer;
        [rows, cols] = size(trialRows);
        for j = 1:rows
            if shuffleGood == 1
                if (oneback_size == trialRows(j, 5)) && ...
                        (oneback_size == twoback_size)
                    shuffleGood = 0;
                end
                twoback_rwd = oneback_rwd;
                oneback_rwd = trialRows(j, 6);
                twoback_size = oneback_size;
                oneback_size = trialRows(j, 5);
            end
        end
    end
    oneback_rwd_outer = oneback_rwd;
    twoback_rwd_outer = twoback_rwd;
    oneback_size_outer = oneback_size;
    twoback_size_outer = twoback_size;
    thisBlock = trialRows;
    
    iteration = 1;
    for t = thisTrialNum:(thisTrialNum + noFBbaseline_length - 1)
        trialContents(t,:) = thisBlock(iteration, :);
        iteration = iteration + 1;
    end
    thisTrialNum = (thisTrialNum + noFBbaseline_length);
end
lastNoFBBLTrial = thisTrialNum - 1;

% tutorial trials (x2)
order = randperm(4);
thisBlock = clampTutorialTrialTypes(order, :);
iteration = 1;
for t = thisTrialNum:(thisTrialNum + 3)
    trialContents(t,:) = thisBlock(iteration, :);
    iteration = iteration + 1;
end
thisTrialNum = (thisTrialNum + 4);
lastTutorialTrial = thisTrialNum - 1;


% test block
oneback_rotn_outer = NaN;
twoback_rotn_outer = NaN;
for i = 1:numTestBlocks
    shuffleGood = 0;
    while shuffleGood == 0        
        shuffleGood = 1;
        
        almostOrder = randperm(length(testBlockTrialTypes));
        tempTrialTypes = testBlockTrialTypes(almostOrder,:);
        tempTrialTypes = [tempTrialTypes,...
            zeros(length(tempTrialTypes),1)];
        
        trialRows = [];
        targetOrder = [];
        for r = 1:4 %each size now has four trial types associated to it rot x rwd
            order = randperm(2);
            theseSizes = targetSizes(order);
            for t = 1:2 %t now stands for target SIZE
                [rows, cols] = size(tempTrialTypes);
                found = 0;
                for p = 1:rows
                    if (found == 0) && (tempTrialTypes(p,7) == 0) ...
                            && (tempTrialTypes(p,5) == theseSizes(t))
                        thisRow = tempTrialTypes(p,1:6);
                        tempTrialTypes(p,7) = 1;
                        found = 1;
                    end
                end
                trialRows = [trialRows; thisRow];
            end
        end
        
        
        oneback_rotn = oneback_rotn_outer;
        twoback_rotn = twoback_rotn_outer;
        oneback_rwd = oneback_rwd_outer;
        twoback_rwd = twoback_rwd_outer;
        [rows, cols] = size(trialRows);
        for j = 1:rows
            if shuffleGood == 1
                if (oneback_rotn == trialRows(j, 3)) && ...
                        (oneback_rotn == twoback_rotn)
                    shuffleGood = 0;
                end
                if (oneback_rwd == trialRows(j, 6)) && ...
                        (oneback_rwd == twoback_rwd)
                    shuffleGood = 0;
                end
                twoback_rotn = oneback_rotn;
                oneback_rotn = trialRows(j, 3);
                
                twoback_rwd = oneback_rwd;
                oneback_rwd = trialRows(j, 6);
            end
        end
    end
    oneback_rotn_outer = oneback_rotn;
    twoback_rotn_outer = twoback_rotn;
    oneback_rwd_outer = oneback_rwd;
    twoback_rwd_outer = twoback_rwd;
    thisBlock = trialRows;
    
    iteration = 1;
    for t = thisTrialNum:3:(thisTrialNum + length(testBlockTrialTypes)*3 -1)
        trialContents(t,1) = 0;
        trialContents(t,2) = 0;
        trialContents(t,3) = 0;
        trialContents(t,4) = thisBlock(iteration, 4);
        trialContents(t,5) = thisBlock(iteration, 5);
        trialContents(t,6) = thisBlock(iteration, 6);
        trialContents(t+1,:) = thisBlock(iteration, :);
        trialContents(t+2,1) = 0;
        trialContents(t+2,2) = 0;
        trialContents(t+2,3) = 0;
        trialContents(t+2,4) = thisBlock(iteration, 4);
        trialContents(t+2,5) = thisBlock(iteration, 5);
        trialContents(t+2,6) = thisBlock(iteration, 6);
        iteration = iteration + 1;
    end
    thisTrialNum = (thisTrialNum + length(testBlockTrialTypes)*3);
end


headers = {'trials', 'endpoint_feedback', 'online_feedback',...
    'aim_report', 'aim_ring', 'rotation', 'instruction',...
    'target_angle', 'target_distance', 'size', 'aim_at_targ', 'reward'};
[rows, cols] = size(trialContents);
trials = nan(rows,1);
endpoint_feedback = nan(rows,1);
online_feedback = nan(rows,1);
aim_report = zeros(rows,1);
aim_ring = zeros(rows,1);
rotation = nan(rows,1);
instruction = nan(rows,1);
target_angle = nan(rows,1);
target_distance = nan(rows,1);
size = nan(rows,1);
aim_at_targ = ones(rows,1);
reward = nan(rows, 1);
for r = 1:rows
    trials(r,1) = r;
    target_distance(r,1) = 80;
    if (r>lastNoFBBLTrial) && (r<=lastTutorialTrial)
        aim_at_targ(r,1) = 0;
    end
    if r==lastBLTrial
        instruction(r,1) = 1;
    elseif r == lastNoFBBLTrial
        instruction(r,1) = 2;
    elseif r == lastTutorialTrial
        instruction(r,1) = 3;
    else
        instruction(r,1) = 0;
    end
    
    endpoint_feedback(r,1) = trialContents(r,1);
    online_feedback(r,1) = trialContents(r,2);
    rotation(r,1) = trialContents(r,3);
    target_angle(r,1) = trialContents(r,4);
    size(r,1) = trialContents(r,5);
    reward(r, 1) = trialContents(r, 6);
end
trialMatrix = [trials, endpoint_feedback, online_feedback,...
    aim_report, aim_ring, rotation, instruction,...
    target_angle, target_distance, size, aim_at_targ, reward];

saveme = [headers; num2cell(trialMatrix)];

cd('C:\Users\labadmin\Desktop\example Olivia STL code\TargetFiles')
writetable(cell2table(num2cell(trialMatrix),...
    'VariableNames', headers), 'STL_Money_T016.txt', 'Delimiter', '\t')