%% EEGLAB Preproc - UGO 
% 1. Load & Channel Fix
% 2. Rimozione DEFINITIVA Canali Extra (ECG, EGG, VEOG...)
% 3. Filtro 0.5-45 Hz
% 4. ASR (Soglia 100 + Rejection ON + Window OFF)
% 5. ICA (Picard) + ICLabel
% 6. Interpolazione (SOLO canali cerebrali mancanti, NO physio)
% 7. Reportistica

clear; close all; clc;
eeglab_root = fileparts(which('eeglab'));
% if isempty(eeglab_root), error('EEGLAB non trovato!'); end
% addpath(eeglab_root); 
eeglab nogui;

% --- CONFIGURAZIONE ---
datapath = 'C:\Users\ugo_p\Desktop\PAC_tryouts\Data_Pill\';
output_dir = 'C:\Users\ugo_p\Desktop\PAC_tryouts\Data_Pill\EEG_preproc_ASR_100\';
report_file = fullfile(output_dir, 'cleaning_report_log_ugo.txt');

if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% Inizializza report
if ~exist(report_file, 'file')
    fid = fopen(report_file, 'w');
    fprintf(fid, 'Subject\tCondition\tOriginal_Sec\tCleaned_Sec\tRemoved_Sec\tKept_Percent\n');
    fclose(fid);
end

% LISTA SOGGETTI COMPLETA
% subject_list = {'P001','P002', 'P003', 'P004', 'P005', 'P006', 'P007', 'P008', 'P009', 'P010', ...
%                 'P011', 'P012', 'P013', 'P014', 'P015', 'P017', 'P018', ...
%                 'P019', 'P021', 'P022', 'P023', 'P024',...
%                 'P025', 'P026', 'P027', 'P028', 'P029', 'P030', 'P031', ...
%                 'P032','P033', 'P034', 'P035', 'P036', 'P037', 'P038', ...
%                 'P039','P040', 'P0166', 'P0200'};

subject_list = {'P018'};

chanlocs_file = fullfile(eeglab_root, 'plugins', 'dipfit', 'standard_BEM', 'elec', 'standard_1005.elc');

% PAROLE CHIAVE PER RIMOZIONE (Incluso VEOG per beccare VEOGexFT10)
physio_keywords = {'egg', 'ecg', 'ekg', 'heog', 'veog', 'eog', 'emg', 'resp', 'gsr'};

% --- PARAMETRI ---
param_filter_locut = 0.2;  % cut off at 0.1
param_filter_hicut = 43;   % cut off at 48


param_asr_burst    = 20;  % Soglia ASR
param_heart_thresh = 0.80; % ICLabel Heart
param_eye_thresh   = 0.80; % ICLabel Eye

% =========================================================================
% LOOP SOGGETTI
% =========================================================================
for s_idx = 1:length(subject_list)
    subject = subject_list{s_idx};
    fprintf('\n\n>>> PROCESSING %s <<<\n', subject);
    
    subj_folder = fullfile(datapath, subject, 'BASELINE');
%     if ~exist(subj_folder, 'dir'), fprintf('Folder not found. Skip.\n'); continue; end
    if ~exist(subj_folder, 'dir'), fprintf('Folder not found. Skip.\n'); end
    
    files = dir(fullfile(subj_folder, '*.vhdr'));
    
    for f = 1:length(files)
        file_name = files(f).name;
        
        if contains(file_name, 'PRECOL', 'IgnoreCase', true), cond = 'pre';
        elseif contains(file_name, 'OA_ST', 'IgnoreCase', true), cond = 'post';
        else, continue; end
        
        fprintf('   --> File: %s (%s)\n', file_name, cond);
        
        try
            %% 1. LOAD & LABEL FIX
            EEG = pop_loadbv(subj_folder, file_name);
            
            fcz_idx = find(contains({EEG.chanlocs.labels}, 'FCZexTP9', 'IgnoreCase', true));
            if ~isempty(fcz_idx), EEG = pop_chanedit(EEG, 'changefield', {fcz_idx, 'labels', 'FCz'}); end
            EEG = pop_chanedit(EEG, 'lookup', chanlocs_file);
            
            % SALVIAMO LA LISTA ORIGINALE PER DEFINIRE I "BRAIN CHANNELS"
            orig_chanlocs = EEG.chanlocs;
            
            %% 2. RIMOZIONE FISICA CANALI EXTRA
            % Rimuove tutto ciò che contiene 'ecg', 'veog' (quindi anche VEOGexFT10), etc.
            keep_mask = true(1, EEG.nbchan);
            removed_list = {};
            for k=1:EEG.nbchan
                if contains(EEG.chanlocs(k).labels, physio_keywords, 'IgnoreCase', true)
                    keep_mask(k) = false;
                    removed_list{end+1} = EEG.chanlocs(k).labels;
                end
            end
            
            if ~isempty(removed_list)
                EEG = pop_select(EEG, 'channel', find(keep_mask));
                fprintf('       [CLEAN] Rimossi per sempre: %s\n', strjoin(removed_list, ', '));
            end
            
            % Definiamo ORA quali sono i canali cerebrali target per l'interpolazione futura.
            % Escludiamo dalla lista originale quelli che abbiamo appena deciso di rimuovere.
            brain_mask_for_interp = ~contains({orig_chanlocs.labels}, physio_keywords, 'IgnoreCase', true);
            final_brain_chanlocs = orig_chanlocs(brain_mask_for_interp);
            
            %% 3. FILTERING 
%             EEG = pop_eegfiltnew(EEG, 'locutoff', param_filter_locut, 'hicutoff', param_filter_hicut);
            EEG = pop_eegfiltnew(EEG, 'locutoff', param_filter_locut); 
            EEG = pop_eegfiltnew(EEG, 'hicutoff', param_filter_hicut); 
            
            %% 4. channel rejection just based on correlation with neighbours
            pnts_pre = EEG.pnts;
            sec_pre = pnts_pre / EEG.srate;
            
%             fprintf('       Running ASR (Threshold=%d, Rejection=ON)...\n', param_asr_burst);
            
            EEG_cleanraw = pop_clean_rawdata(EEG, ...
                'FlatlineCriterion', 'off', ...      % OFF
                'ChannelCriterion', 0.8, ...         % Remove channels correlating <80% with neighbours
                'LineNoiseCriterion', 'off', ...     % OFF
                'Highpass', 'off', ...               % OFF
                'BurstCriterion', param_asr_burst, ...         
                'BurstRejection', 'on', ...         
                'Distance', 'Euclidian', ...
                'WindowCriterion', 'off');           % OFF

            % Reportistica
            pnts_post = EEG_cleanraw.pnts;
            sec_post = pnts_post / EEG_cleanraw.srate;
            perc_kept = (pnts_post / pnts_pre) * 100;
            sec_removed = sec_pre - sec_post;
            
            fprintf('       DATA KEPT: %.2f%%\n', perc_kept);
            
            fid = fopen(report_file, 'a');
            fprintf(fid, '%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f%%\n', ...
                    subject, cond, sec_pre, sec_post, sec_removed, perc_kept);
            fclose(fid);
            
            %% 5. ICA (PICARD)
%             EEG = pop_reref(EEG, []);
%             data_rank = rank(double(EEG.data'));
%             
%             fprintf('       Running Picard ICA (Rank %d)...\n', data_rank);
%             plugin_askinstall('picard', 'picard', 1); % install Picard plugin [UGO]
            % High-pass filter (FIR filter with 1.5 Hz cutoff)
            % Note: EEGLAB's pop_eegfiltnew uses cutoff at -6 dB point.
            EEG = pop_eegfiltnew(EEG_cleanraw, 'locutoff', 2.5); 
        
%             EEG = eeg_picard(EEG, 'mode', 'standard', 'pca', data_rank, 'verbose', false);
%             EEG =  pop_runica(EEG, 'extended', 1, 'interupt', 'on');


            data_rank = rank(double(EEG.data'));
            EEG = pop_runica(EEG, ...
            'extended', 1, ...
            'pca', data_rank, ...
            'interupt', 'on');
            
            %% Transfer ICA weights to EEG_cleanraw (unfiltered dataset)
            EEG_cleanraw.icaweights = EEG.icaweights;
            EEG_cleanraw.icasphere  = EEG.icasphere;
            EEG_cleanraw.icawinv    = EEG.icawinv;
            EEG_cleanraw.icaact     = [];
            EEG_cleanraw.icachansind = EEG.icachansind;
            
            %% Recompute ICA activations on cleanraw data
            EEG_cleanraw = eeg_checkset(EEG_cleanraw, 'ica');

            %% 6. ICLABEL
            EEG = pop_iclabel(EEG_cleanraw, 'default');
            flag_thresholds = [NaN NaN; ...             % Brain
                               NaN NaN; ...              % Muscle
                               param_eye_thresh 1; ...  % Eye
                               param_heart_thresh 1; ...% Heart
                               NaN NaN; ...             % Line Noise
                               NaN NaN; ...              % Channel Noise
                               NaN NaN];                % Other
            EEG = pop_icflag(EEG, flag_thresholds);
            EEG = pop_subcomp(EEG, [], 0);
            
            %% 7. INTERPOLAZIONE SICURA (SOLO BRAIN)
            % Qui usiamo 'final_brain_chanlocs' calcolato allo step 2.
            % Questo garantisce che ECG, VEOG, EGG NON vengano mai reinseriti.
            % Vengono interpolati solo i canali cerebrali (es. Fz, C3) che ASR (Flatline) potrebbe aver rimosso.
            
            fprintf('       Interpolazione canali mancanti (Solo Cerebrali)...\n');
            EEG = pop_interp(EEG, final_brain_chanlocs);
            
            % NOW we rereference
            EEG = pop_reref(EEG, []); % Reref finale
            
            %% 8. SAVE
            out_name = sprintf('%s_%s_cleaned_ugo.set', subject, cond);
            EEG = pop_saveset(EEG, 'filename', out_name, 'filepath', output_dir);
            
        catch ME
            fprintf('       [ERROR] %s: %s\n', file_name, ME.message);
            fid = fopen(report_file, 'a');
            fprintf(fid, '%s\t%s\tERROR\tERROR\t%s\t0%%\n', subject, cond, ME.message);
            fclose(fid);
        end
    end
end
fprintf('\nDONE. Report salvato in: %s\n', report_file);

