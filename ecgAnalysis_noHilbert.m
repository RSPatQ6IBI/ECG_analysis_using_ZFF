path_zff = '/idiap/user/rprasad/Projects/ECG_analysis_matlab_codes/FeatFrmZFRnZTW/';
path_code = '/idiap/user/rprasad/Projects/ECG_analysis_matlab_codes/';
path_data = '/idiap/temp/rprasad/datasets_local/ECG_datasets/mit-bih-master/sig/';
path_ann = '/idiap/temp/rprasad/datasets_local/ECG_datasets/mit-bih-master/atrs/';
dir_data = dir(path_data);
addpath( path_code , path_data , path_ann , path_zff );
fs = 360;
%%
se_arr = []
for count = 10:50 %:length(dir_data)
    fileName = dir_data(count).name();
    %     [fileName num2str(count) ]
    x = readmatrix([path_data fileName]);
    ecg = x(:,2)/max(abs(x(:,2)));
    [s1,s2] = strtok(fileName,'.');
    a = readtable([path_ann s1 'atr' s2]);
    a = a{:,'Var2'};
    
    % to check the plots
    flag_hil = 0;
    flag_tp_fp = 1 ;
    flag_plot = 1; plot_err = 1;
    tol_win = 40;
     
    % 2 for false positives in GT and derived ;
    % 1 for true positive in both
    [res ,ax ] = display_results(ecg, fs , a, flag_hil, flag_tp_fp , flag_plot);
    se_val = compute_acc(res, a, tol_win, plot_err, ecg , ax);
    se_arr = [se_arr; se_val];
end
disp('<<<<<__________________>>>>>>>>')
disp(se_arr)

%%% ZFF computation
function[gci,es] = derive_zff(x,fs)
[~,gci,es,~] = zfsig_ECG(x,fs,650);
end

function[gci , es , d_ecg ] = qrs_with_hilbert(ecg , fs , flag_hil)
d_ecg = diff(ecg);
if(flag_hil)
    d_ecg = abs(hilbert(ecg.*[0;d_ecg]));
else
    d_ecg = abs(ecg.*[0;d_ecg]);
end
[gci,es] = derive_zff(d_ecg,fs);

end

function[ ret_arr ] = validate_ann(gt_ann , gci , es , ecg , d_ecg , flag_tp_fp )
runMean_factor = 10;
% for derived locations
gci_fp = gci(es < RunMean( es , runMean_factor )/5 );
gci_tp = gci(es > RunMean( es , runMean_factor )/5 );

% for ground truth locations
gt_ener_ecg =  d_ecg(gt_ann);
fp_gt = gt_ann(gt_ener_ecg < RunMean( gt_ener_ecg , runMean_factor )/5);
tp_gt = gt_ann(gt_ener_ecg > RunMean( gt_ener_ecg , runMean_factor )/5);

if(flag_tp_fp == 1)
    line_tp_gt = zeros(length(ecg) , 1);
    line_tp_gt(tp_gt) = 1;
    line_tp_der = zeros(length(ecg) , 1);
    line_tp_der(gci_tp) = 1;
    ret_arr = [line_tp_gt line_tp_der];
elseif( flag_tp_fp == 2 )
    line_fp_gt = zeros(length(ecg) , 1);
    line_fp_gt(fp_gt) = 1;
    line_fp_der = zeros(length(ecg) , 1);
    line_fp_der(gci_fp) = 1;
    ret_arr = [line_fp_gt line_fp_der];
end

end

function[ret_arr, ax] = display_results(ecg, fs , ann , flag_hil, flag_tp_fp , flag_plot)

[ gci , es , d_ecg ] = qrs_with_hilbert(ecg , fs , flag_hil );
[ ret_arr ] = validate_ann( ann , gci , es , ecg , d_ecg , flag_tp_fp );
line_fp_gt = ret_arr(:,1);
line_fp_der = ret_arr(:,2);
line_baseline = zeros(length(ecg),1); line_baseline(ann)=1;
if(flag_plot)
    clf;
    ax(1)=subplot(411); plot(ecg); hold on; plot(line_fp_gt,'--'); grid; %xlim([1.25e5 1.31e5]);
    title('ground truth annotation'); hold on; plot(line_baseline/2,'.k'); 
    ax(2)=subplot(412); plot(ecg); hold on; plot(line_fp_der,'--r'); grid; %xlim([1.25e5 1.31e5]);
    title('derived annotation');
%     linkaxes(ax,'x');
    hold on; 
end
end

function[se] = compute_acc(res , ann, tol_win , plot_err, ecg , ax)
truth_ann = find(res(:,1));
derived_ann = find(res(:,2));

disp([ 'The no. of ground truth annotations = ' num2str( length(ann) ) ])
disp([ 'The no. of correct ground truth annotations = ' num2str( length(truth_ann) ) ])
disp([ 'The no. of identified annotations = ' num2str( length(derived_ann) ) ])

tp = 0; fn = 0; fn_arr=[]; tp_arr=[];
for i=1:length(ann)
    err = min(abs((derived_ann - ann(i) )));
    if(err < tol_win)
        tp = tp+1;
        tp_arr = [tp_arr; ann(i)];
    else
        fn = fn+1;
        fn_arr = [fn_arr; ann(i)];
    end
end
fp = max([0 length(derived_ann)-length(ann) ]);

disp(['Results : true pos - ' num2str(tp) ' -- false neg - ' num2str(fn) ' -- false pos - ' num2str(fp) ]);
disp('============================================')
se = tp/(tp+fn); 

line_tp = zeros(length(ecg),1); line_fn = zeros(length(ecg),1);
% line_fp = zeros(length(ecg),1); 
line_tp(tp_arr) = 1; line_fn(fn_arr) = 1;
% line_fp(fp_arr)=1; 
if(plot_err)
    
    ax(3)=subplot(413); plot(ecg); hold on; plot(line_tp,'--'); grid; %xlim([1.25e5 1.31e5]);
    title('true positive derived');
    ax(4)=subplot(414); plot(ecg); grid; %xlim([1.25e5 1.31e5]);
    hold on; plot(line_fn,'--g'); title('false negative derived'); 
    linkaxes(ax,'x');
end
end