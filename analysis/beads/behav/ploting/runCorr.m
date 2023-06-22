function fh = runCorr(mXsim, mXfit, conditions)


global AZred AZblue AZcactus AZsky AZriver AZsand AZmesa AZbrick

AZred = [171,5,32]/256;
AZblue = [12,35,75]/256;
AZcactus = [92, 135, 39]/256;
AZsky = [132, 210, 226]/256;
AZriver = [7, 104, 115]/256;
AZsand = [241, 158, 31]/256;
AZmesa = [183, 85, 39]/256;
AZbrick = [74, 48, 39]/256;

names       = {'Cost sample' 'softmax temperature'};
symbols     = {'Cs' '\beta'};


for c = 1:conditions

    cond_ssX = mXsim{1,c};
    cond_ffX = mXfit{1,c};

    fh{c}    = figure('Name','Paramteter correlations');
    set(fh{c}, 'Position', [811   613   600   300]);

    [~,~,~,ax] = easy_gridOfEqualFigures([0.2  0.1], [0.1 0.18 0.04]);

    % loop over parameter values
    for i = 1:size(cond_ssX,1)

        axes(ax(i)); hold on;
        plot(cond_ssX(i,:), cond_ffX(i,:), 'o', 'color', AZred, 'markersize', 8, 'linewidth', 1)
        xl = get(gca, 'xlim');
        plot(xl, xl, 'k--')

    end % end of param values loop 

%    % find 'bad' Cs values
%     thresh  = 0.4;
%     ind     = abs(cond_ssX(1,:) - cond_ffX(1,:)) > thresh;
% 
%     for i = 1:2
%         axes(ax(i));
%         plot(cond_ssX(i,ind), cond_ffX(i,ind), 'o', 'color', AZblue, 'markersize', 8, 'linewidth', 1, ...
%         'markerfacecolor', [1 1 1]*0.5)
%     end

    %set(ax(1,2),'xscale', 'log', 'yscale' ,'log')
    
    axes(ax(1)); t      = title('Cost Sample');
    axes(ax(2)); t(2)   = title('softmax temperature');
    
    axes(ax(1)); xlabel('simulated Cs'); ylabel('fit Cs'); 
    axes(ax(2)); xlabel('simulated \beta'); ylabel('fit \beta'); 

    set(ax, 'tickdir', 'out', 'fontsize', 18)
    set(t, 'fontweight', 'normal')
    addABCs(ax(1), [-0.07 0.08], 32)
    addABCs(ax(2), [-0.1 0.08], 32, 'B')
    set(ax, 'tickdir', 'out')
    
    for i = 1:size(cond_ssX,1)
        axes(ax(i));
        xl = get(gca, 'xlim');
        plot(xl, xl, 'k--')
    end

clear ind
end % end of conditions loop


end