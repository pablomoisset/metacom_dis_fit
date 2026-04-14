close all

% This is the only line that should be changed to select the input file
pattern = '../results_paper/R6_51469*.mat' ;

expression = '.*_rep0' ;


%responses = {'alpha_uM', 'beta_uM', 'gamma_uM',...
%             'alpha_wM', 'beta_wM', 'gamma_wM'} ;

%This selects the  response variables to be plotted
responses = {'alpha_uM', 'beta_uM', 'gamma_uM' } ;

% Titles for the plots
%response_titles = {'STDERR $D^u_\alpha$', 'STDERR $D^u_\beta$', 'STDERR $D^u_\gamma$' } ;
response_titles = {'$\mathcal{D}^u_\alpha$', '$\mathcal{D}^u_\beta$', '$\mathcal{D}^u_\gamma$' } ;

%responses = {'alpha_uM', 'eveness', 'gamma_uM',...
%             'alpha_wM', 'beta_JwM', 'gamma_wM'} ;
         
%responses = {'gamma_biomass'} ;

 
P1 = 'adapt.pars.kappaW' ; %y axis
P2 = 'D.pars.D0' ; %x axis
logspace_P1 = false ;
logspace_P2 = true ;

fresponse = @nanmean ;
%fresponse = @(A) std(A,"omitnan") ;
%fresponse = @(A) std(A,"omitnan")/(nanmean(A)+realmin) ;
%fresponse = @(A) std(A,"omitnan")/sqrt(sum(isfinite(A))) ;

[bmp, unique_P1s, unique_P2s,p] = load_data(pattern, expression, responses, P1, P2, fresponse);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Heatmap
%%{
if p.W.pars.nbasals == p.W.pars.nspecies
    web_type = 'Comp.' ;
else
    web_type = 'Fw.' ;
end
figure('Name', sprintf('Type=%s, beta=%g, comp=%g, Tnoise=%g, Tr=%g', web_type, p.D.beta(1), p.W.alpha(1,1), p.L.pars.T_noise, (p.L.pars.T_max-p.L.pars.T_min)/2) );

%trans_m = [1 4 2 5 3 6] ;
trans_m = [1 2 3] ;

for j = 1:numel(responses)
    i = trans_m(j) ;
    if ~isempty(responses{i})
%        subplot(3,2,j);
        subplot(1,3,j);
%        imagesc((bmp{i}-380)*63/(1090-380)+1) ;
        bmp_i = flipud(bmp{i}) ;
        
        imagesc(bmp_i) ;
        min(bmp{i}(:))
        max(bmp{i}(:))
         %zlim([300 1200])
        ax = gca;
%        ax.CLim=[1050 1090];
        %ax.CLim=[820 1050];
        %ax.CLim=[330 870];
        colorbar;
        title(response_titles{i}, 'interpreter', 'latex', 'fontsize', 17) ;
        
        y = linspace(1,size(bmp{i},1),5) ;
        yticks(y) ;
        if logspace_P1
            s = sprintf('%.1e ',round(logspace(log10(max(unique_P1s)),log10(min(unique_P1s)),numel(y))',9)) ;
        else
            s = sprintf('%.1e ',round(linspace(max(unique_P1s),min(unique_P1s),numel(y))',9)) ;
        end
        yticklabels(strsplit(s(1:end-1)));
        ax.YAxis.FontSize = 12;
%        ylabel(P1, 'interpreter', 'latex');
        ylabel('$\kappa^W$' , 'interpreter', 'latex', 'fontsize', 15);
        
        x = linspace(1,size(bmp{i},2),5) ;
        xticks(x);
        
        if logspace_P2
            s = sprintf('%.1e ',round(logspace(log10(min(unique_P2s)),log10(max(unique_P2s)),numel(x))',8)) ;
        else
            s = sprintf('%.1e ',round(linspace(min(unique_P2s),max(unique_P2s),numel(x)),8));
        end
        
        xticklabels(strsplit(s(1:end-1)));
        ax.XAxis.FontSize = 12;
        %xlabel(P2, 'interpreter', 'latex');
        xlabel('$D_0$' , 'interpreter', 'latex', 'fontsize', 15);
%        ax.XAxis.FontSize = 15;
    end
end
%set(gcf, 'Position',  [100, 100, 760, 1000])
%set(gcf, 'Position',  [100, 100, 300, 230])
%set(gcf, 'Position', [-5.2471e+02 9.4043e+02 6.6171e+02 6.4686e+02])
set(gcf, 'Position',[6 237 1882 420])
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bmp, unique_P1s, unique_P2s, first_p] = load_data(name_pattern, filter_expr, response_vars, P1, P2, fresponse)

files = dir(name_pattern) ;

P1s = [] ;
P2s = [] ;

response_mean = cell(1,numel(response_vars)) ;
%response_means = cell(1,numel(response_vars)) ;

for i = 1:numel(files)
    i
    startIndex = regexp(files(i).name, filter_expr,'ONCE');
    if ~isempty(startIndex)
        files(i).name
        load([files(i).folder '/' files(i).name] ) ;
        q = 2 ;
        
        for r = 1:numel(results.M)
            tmpM = results.M{r} ;
            tmpM(tmpM<1e-7) = 0 ;
            species_per_site = sum(tmpM>0) ;
            [Dq_alpha, Dq_beta, Dq_gamma, true_alpha_wM, true_beta_wM, true_gamma_wM, shannon] = hill_numbers(tmpM,q) ;
            results.Dq_alpha_list(r) = Dq_alpha ;
            results.Dq_beta_list(r) = Dq_beta ;
            results.Dq_gamma_list(r) = Dq_gamma ;
            results.true_alpha_wM_list(r) = true_alpha_wM ;
            results.true_beta_wM_list(r)  = true_beta_wM ;
            results.true_gamma_wM_list(r) = true_gamma_wM ;
            results.eveness_list(r) = mean(shannon./species_per_site) ;
        end
        P1s = [P1s; eval(['all_p(1).p.' P1])] ;
        P2s = [P2s; eval(['all_p(1).p.' P2])] ;
        for name_idx=1:numel(response_vars)
            response_vars{name_idx}
            if ~isempty(response_vars{name_idx})
                tmp = eval(['results.' response_vars{name_idx} '_list'])	;
                response_mean{name_idx} = [response_mean{name_idx}; fresponse(tmp)] ;
                %response_means{name_idx} = [response_means{name_idx}, tmp] ;                
            end
        end
    end
end
%p_alpha_mean = log10(D0s)
unique_P1s = unique(P1s) ;
unique_P2s = unique(P2s) ;

bmp = cell(1,numel(response_vars)) ;

for name_idx=1:numel(response_vars)
    if ~isempty(response_vars{name_idx})
        bmp{name_idx} = zeros(numel(unique_P1s),numel(unique_P2s))+min(response_mean{name_idx}(:)) ;
        
        for k = 1:numel(P1s)
            ii = find(P1s(k)==unique_P1s) ;
            jj = find(P2s(k)==unique_P2s) ;
            bmp{name_idx}(ii,jj) = response_mean{name_idx}(k);
        end
    end
end
first_p = all_p(1).p ;
end
