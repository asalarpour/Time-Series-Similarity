function TS_computeDistanceforDatasetAllInMex(datasetNum, distNum, divideNum)


% ######################################################################## 
% set addresses
folder = TS_getAddress();

% Import dataset
[tsSet, tsSpec] = TS_selectDataset(datasetNum, folder); %#ok<ASGLU>
    
% local variable
distCell = cell(0);
reverseStr = ''; %#ok<NASGU>

% ######################################################################## 

switch distNum

    case 1  % Calculate DTW_FULL 
        
        winRatio = -1;
        tic
        [dtwFullDist, dtwFullLen, dtwFullLenMax] = ...
            TS_dtwDistanceMexAll(tsSet, winRatio);
        dtwFullTotalTime = toc;
        
        if ( any(isinf(dtwFullDist)))
            error('there is error in dtwfull');
        end

        distCell{1}.distMat = squareform(dtwFullDist);
        distCell{1}.dtwLen = squareform(dtwFullLen);
        distCell{1}.lenMax = squareform(dtwFullLenMax);
        distCell{1}.totalTime = dtwFullTotalTime;
        distCell{1}.winRatio = winRatio;
        
        nameSave = [folder.distance 'Dataset_' num2str(datasetNum) '_' num2str(distNum) '_DTWFull_mex.mat'];
        save(nameSave, 'distCell', '-v7.3')
        
% ########################################################################
        
    case 2  % Calculate DTW_CONSTRAINED 
        winInterval = linspace(1, 30, 30) / 100;
        for winIdx = 1: length(winInterval)
            
            winRatio = winInterval(winIdx);
            tic
            [dtwConDist, dtwConLen, dtwConLenMax] = ...
                TS_dtwDistanceMexAll(tsSet, winRatio);
            dtwConTotalTime = toc;
            
            if ( any(isinf(dtwConDist)))
                error('there is error in dtwcon');
            end
            
            distCell{1}.distMat = squareform(dtwConDist);
            distCell{1}.dtwLen = squareform(dtwConLen);
            distCell{1}.lenMax = squareform(dtwConLenMax);
            distCell{1}.totalTime = dtwConTotalTime;
            distCell{1}.winRatio = winRatio;
            distCell{1}.winIdx = winIdx;

            nameSave = [folder.distance 'Dataset_' num2str(datasetNum)...
                '_' num2str(distNum) '_DTWCon_winRatio_' num2str(winIdx) '_mex.mat'];
            save(nameSave, 'distCell', '-v7.3')
        end

% ########################################################################
                
    case 3  % Calculate LCSS_EPS 

        tmpPt = [];
        for k = 1: size(tsSet, 1)
            tmpPt = [tmpPt; tsSet(k).ts]; %#ok<AGROW>
        end
        tmpF = mean(tmpPt) + std(tmpPt);
        tmpE = mean(tmpPt) - std(tmpPt);
        
        epsMax = sqrt(sum((tmpF- tmpE) .^ 2));
        epsMin = 0.02 * epsMax;
        epsInterv = linspace(epsMin, epsMax, divideNum);
        
        for epsIdx = 1: length(epsInterv)
            
            eps = epsInterv(epsIdx);
            tic
            lcssDist = TS_lcssDistanceMexAll(tsSet, eps);
            lcssTotalTime = toc;
            
            if ( any(isinf(lcssDist)))
                error('there is error in lcss');
            end
            
            distCell{1}.distMat = squareform(lcssDist);
            distCell{1}.totalTime = lcssTotalTime;
            distCell{1}.eps = eps;
            distCell{1}.epsIdx = epsIdx;

            nameSave = [folder.distance 'Dataset_' num2str(datasetNum)...
                '_' num2str(distNum) '_LCSS_eps_' num2str(epsIdx) '_mex.mat'];
            save(nameSave, 'distCell', '-v7.3')
        end
        
% ########################################################################

    case 4  % Calculate EDR_EPS 
        
        tmpPt = [];
        for k = 1: size(tsSet, 1)
            tmpPt = [tmpPt; tsSet(k).ts]; %#ok<AGROW>
        end
        tmpF = mean(tmpPt) + std(tmpPt);
        tmpE = mean(tmpPt) - std(tmpPt);
        
        epsMax = sqrt(sum((tmpF- tmpE) .^ 2));
        epsMin = 0.02 * epsMax;
        epsInterv = linspace(epsMin, epsMax, divideNum);
        
        for epsIdx = 1: length(epsInterv)
            
            eps = epsInterv(epsIdx);
            tic
            [edrDist, edrLenMax] = TS_edrDistanceMexAll(tsSet, eps);
            edrTotalTime = toc;
            
            if ( any(isinf(edrDist)))
                error('there is error in edr');
            end
            
            distCell{1}.distMat = squareform(edrDist);
            distCell{1}.lenMax = squareform(edrLenMax);
            distCell{1}.totalTime = edrTotalTime;
            distCell{1}.eps = eps;
            distCell{1}.epsIdx = epsIdx;

            nameSave = [folder.distance 'Dataset_' num2str(datasetNum)...
                '_' num2str(distNum) '_EDR_eps_' num2str(epsIdx) '_mex.mat'];
            save(nameSave, 'distCell', '-v7.3')
        end
        
        
% ########################################################################
        
    case 5  % Calculate ERP_GAP 
    
        tmpPt = [];
        for k = 1: size(tsSet, 1)
            tmpPt = [tmpPt; tsSet(k).ts]; %#ok<AGROW>
        end
        tmp1 = mean(tmpPt) + 3 * std(tmpPt);
        tmp2 = mean(tmpPt) - 3 * std(tmpPt);
        gapInterv = [zeros(size(tmp1))' tmp1' tmp2'] ;
        
        for gapIdx = 1: size(gapInterv, 2)
            
            gap = gapInterv(:, gapIdx)';
            tic
            [erpDist, erpLenMax] = TS_erpDistanceMexAll(tsSet, gap);
            erpTotalTime = toc;
            
            if ( any(isinf(erpDist)))
                error('there is error in erp');
            end
            
            distCell{1}.distMat = squareform(erpDist);
            distCell{1}.lenMax = squareform(erpLenMax);
            distCell{1}.totalTime = erpTotalTime;
            distCell{1}.gap = gap;
            distCell{1}.gapIdx = gapIdx;

            nameSave = [folder.distance 'Dataset_' num2str(datasetNum)...
                '_' num2str(distNum) '_ERP_gap_' num2str(gapIdx) '_mex.mat'];
            save(nameSave, 'distCell', '-v7.3')
        end
% ########################################################################

    case 6  % Calculate SWALE_EPS_GAP_reward 
        
        tmpPt = [];
        for k = 1: size(tsSet, 1)
            tmpPt = [tmpPt; tsSet(k).ts]; %#ok<AGROW>
        end
        tmpF = mean(tmpPt) + std(tmpPt);
        tmpE = mean(tmpPt) - std(tmpPt);
        
        epsMax = sqrt(sum((tmpF- tmpE) .^ 2));
        epsMin = 0.02 * epsMax;
        epsInterv = linspace(epsMin, epsMax, divideNum);
        reward = 50 * max(max(tmpPt));
        
        gapMax = reward;
        gapMin = 0;
        gapInterv = linspace(gapMin, gapMax, divideNum);
    

        
        for epsIdx = 1: length(epsInterv)
            for gapIdx = 1: length(gapInterv)
                
                eps = epsInterv(epsIdx);
                gapC = gapInterv(gapIdx);
                tic
                [swaleDist, swaleLenMax] = TS_swaleDistanceMexAll(tsSet, eps, gapC, reward); 
                swaleTotalTime = toc;
                
                if ( any(isinf(swaleDist)))
                    error('there is error in swale');
                end
                distCell{1}.distMat = squareform(swaleDist);
                distCell{1}.lenMax = squareform(swaleLenMax);
                distCell{1}.totalTime = swaleTotalTime;
                distCell{1}.gapC = gapC;
                distCell{1}.gapIdx = gapIdx;
                distCell{1}.eps = eps;
                distCell{1}.epsIdx = epsIdx;
                distCell{1}.reward = reward;
                nameSave = [folder.distance 'Dataset_' num2str(datasetNum)...
                    '_' num2str(distNum) '_SWALE_eps_' num2str(epsIdx) '_gapC_' num2str(gapIdx) '_mex.mat'];
                save(nameSave, 'distCell', '-v7.3')

            end
        end
        
% ########################################################################

    case 7 % Calculate TWED_LAMBDA_NU
        
        tmpPt = [];
        for k = 1: size(tsSet, 1)
            tmpPt = [tmpPt; tsSet(k).ts]; %#ok<AGROW>
        end
        
        lambdaMax = sum(max(tmpPt));
        lambdaMin = 0;
        lambdaInterv = linspace(lambdaMin, lambdaMax, divideNum);
        
        nuInterv = 10 .^ -(0:divideNum - 1);
        
        for lambdaIdx = length(lambdaInterv):-1:1
            for nuIdx = 1: length(nuInterv)
                
                nu = nuInterv(nuIdx);
                lambda = lambdaInterv(lambdaIdx);
                
                tic
                [twedDist, twedLenMax] = TS_twedDistanceMexAll(tsSet, lambda, nu); 
                twedTotalTime = toc;
                
                if ( any(isinf(twedDist)))
                    error('there is error in twed');
                end
                distCell{1}.distMat = squareform(twedDist);
                distCell{1}.lenMax = squareform(twedLenMax);
                distCell{1}.totalTime = twedTotalTime;
                distCell{1}.nu = nu;
                distCell{1}.nuIdx = nuIdx;
                distCell{1}.lambda = lambda;
                distCell{1}.lambdaIdx = lambdaIdx;
                
                nameSave = [folder.distance 'Dataset_' num2str(datasetNum)...
                    '_' num2str(distNum) '_TWED_nu_' num2str(nuIdx) '_lambda_' num2str(lambdaIdx) '_mex.mat'];
                save(nameSave, 'distCell', '-v7.3')

            end
        end
        
% ########################################################################

    case 8 % Calculate MSM_COST 
        
        costInterv = 10 .^ (-2:5);

        for costIdx = 1: length(costInterv)
                
            cost = costInterv(costIdx);
            
            tic
            [msmDist, msmLenMax] = TS_msmDistanceMexAll(tsSet, cost);
            msmTotalTime = toc;
            
            if ( any(isinf(msmDist)))
                error('there is error in msm');
            end
            distCell{1}.distMat = squareform(msmDist);
            distCell{1}.lenMax = squareform(msmLenMax);
            distCell{1}.totalTime = msmTotalTime;
            distCell{1}.cost = cost;
            distCell{1}.costIdx = costIdx;
            
            nameSave = [folder.distance 'Dataset_' num2str(datasetNum)...
                '_' num2str(distNum) '_MSM_cost_' num2str(costIdx) '_mex.mat'];
            save(nameSave, 'distCell', '-v7.3')

        end
        
        
% ########################################################################        

    case 9 % Calculate WDTW_G 
        
        gInterv = [0 0.05 0.25 0.5 1];        
        
        for gIdx = 1: length(gInterv)
                
            g = gInterv(gIdx);
            
            tic
            [wdtwDist, wdtwLen, wdtwLenMax] = ...
                TS_weightedDtwDistanceMexAll(tsSet, g);
            wdtwTotalTime = toc;
            
            if ( any(isinf(wdtwDist)))
                error('there is error in wdtw');
            end
            
            distCell{1}.distMat = squareform(wdtwDist);
            distCell{1}.wdtwLen = squareform(wdtwLen);
            distCell{1}.lenMax = squareform(wdtwLenMax);
            distCell{1}.totalTime = wdtwTotalTime;
            distCell{1}.g = g;
            distCell{1}.gIdx = gIdx;
            
            nameSave = [folder.distance 'Dataset_' num2str(datasetNum)...
                '_' num2str(distNum) '_WDTW_g_' num2str(gIdx) '_mex.mat'];
            save(nameSave, 'distCell', '-v7.3')

        end
        
% ########################################################################          

    case 10  % Calculate Hausdorff 
                
        tic
        hausDist = TS_hausdorffDistanceMexAll(tsSet);
        hausTotalTime = toc;
        
        if ( any(isinf(hausDist)))
            error('there is error in hausdorff');
        end
        distCell{1}.distMat = squareform(hausDist);
        distCell{1}.totalTime = hausTotalTime;

            
        nameSave = [folder.distance 'Dataset_' num2str(datasetNum)...
            '_' num2str(distNum) '_Hausdorff_mex.mat'];
        save(nameSave, 'distCell', '-v7.3')
         
% ########################################################################

    case 11  % Calculate Frechet 
                
        tic
        frechetDist = TS_frechetDistanceMexAll(tsSet);
        frechetTotalTime = toc;
        
        if ( any(isinf(frechetDist)))
            error('there is error in frechet');
        end
        distCell{1}.distMat = squareform(frechetDist);
        distCell{1}.totalTime = frechetTotalTime;

            
        nameSave = [folder.distance 'Dataset_' num2str(datasetNum)...
            '_' num2str(distNum) '_Frechet_mex.mat'];
        save(nameSave, 'distCell', '-v7.3')
         
% ########################################################################

    case 12  % Calculate discreteFrechet 
                
        tic
        disFrechetDist = TS_discreteFrechetDistanceMexAll(tsSet);
        disFrechetTotalTime = toc;
        
        if ( any(isinf(disFrechetDist)))
            error('there is error in discrete frechet');
        end
        distCell{1}.distMat = squareform(disFrechetDist);
        distCell{1}.totalTime = disFrechetTotalTime;

            
        nameSave = [folder.distance 'Dataset_' num2str(datasetNum)...
            '_' num2str(distNum) '_disFrechet_mex.mat'];
        save(nameSave, 'distCell', '-v7.3')
         
% ########################################################################

    case 13  % Calculate SSPD 
                
        tic
        sspdDist = TS_sspdDistanceMexAll(tsSet);
        sspdTotalTime = toc;
        
        if ( any(isinf(sspdDist)))
            error('there is error in sspd');
        end
        distCell{1}.distMat = squareform(sspdDist);
        distCell{1}.totalTime = sspdTotalTime;

            
        nameSave = [folder.distance 'Dataset_' num2str(datasetNum)...
            '_' num2str(distNum) '_SSPD_mex.mat'];
        save(nameSave, 'distCell', '-v7.3')
         
% ########################################################################

    case 14 % Calculate TQuEST_THRESH 
        
        tmpPt = [];
        for k = 1: size(tsSet, 1)
            tmpPt = [tmpPt; tsSet(k).ts]; %#ok<AGROW>
        end
        tmpMax = max(tmpPt);
        
        threshMax = TS_euclideanDistance(tmpMax,[0,0]);
        threshMin = 0;
        threshInterv = linspace(threshMin, threshMax, divideNum);       
                
        for threshIdx = 1: length(threshInterv)
            
            thresh = threshInterv(threshIdx);
            tic
            tquestDist = TS_tquestDistanceMexAll(tsSet, thresh);
            tquestTotalTime = toc;
            
            if ( any(isinf(squareform(tquestDist))) || any(isnan(squareform(tquestDist))))
                error('there is error in tquest');
            end
            
            distCell{1}.distMat = squareform(tquestDist);
            distCell{1}.totalTime = tquestTotalTime;
            distCell{1}.thresh = thresh;
            distCell{1}.threshIdx = threshIdx;

            nameSave = [folder.distance 'Dataset_' num2str(datasetNum)...
                '_' num2str(distNum) '_TQuEST_thresh_' num2str(threshIdx) '_mex.mat'];
            save(nameSave, 'distCell', '-v7.3')
        end
        
% ########################################################################

    case 15 % Calculate CIDTW 
        
        tic
        [cidtwDist, cidtwLen, cidtwLenMax] = TS_complexityInvariantDtwDistanceMexAll(tsSet);
        cidtwTotalTime = toc;
        
        if ( any(isinf(cidtwDist)))
            error('there is error in cidtw');
        end
        distCell{1}.distMat = squareform(cidtwDist);
        distCell{1}.totalTime = cidtwTotalTime;
        distCell{1}.cidtwLen = squareform(cidtwLen);
        distCell{1}.lenMax = squareform(cidtwLenMax);


            
        nameSave = [folder.distance 'Dataset_' num2str(datasetNum)...
            '_' num2str(distNum) '_CIDTW_mex.mat'];
        save(nameSave, 'distCell', '-v7.3')
         
% ########################################################################        
        
    case 16 % Calculate DDTW_ALPHA 
        tsSetDiff = tsSet;
        for i = 1: length(tsSet)
            tsSetDiff(i).ts = diff(tsSet(i).ts, 1, 1);
        end
            
        alphaInterv = linspace(1, pi / 2, divideNum);
                
        for alphaIdx = 1: length(alphaInterv)
            
            alpha = alphaInterv(alphaIdx);
            tic
            dtw1 = TS_dtwDistanceMexAll(tsSet);
            dtw2 = TS_dtwDistanceMexAll(tsSetDiff);
            a = cos(alpha);
            b = sin(alpha);

            ddtwDist = a * dtw1 + b * dtw2;
            ddtwTotalTime = toc;
            
            if ( any(isinf(ddtwDist)))
                error('there is error in ddtw');
            end
            
            distCell{1}.distMat = squareform(ddtwDist);
            distCell{1}.totalTime = ddtwTotalTime;
            distCell{1}.alpha = alpha;
            distCell{1}.alphaIdx = alphaIdx;

            nameSave = [folder.distance 'Dataset_' num2str(datasetNum)...
                '_' num2str(distNum) '_DDTW_alpha_' num2str(alphaIdx) '_mex.mat'];
            save(nameSave, 'distCell', '-v7.3')
        end
                
% ########################################################################        
 
end


