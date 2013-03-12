% (c) 2012-2013 M Schaub -- michael.schaub09@imperial.ac.uk
function MarkovZoomingMap(A,filename,time,mode,nr_runs)
%MARKOVZOOMINGMAP implements a version of the Markov Zooming Map equation method
% as discussed in the article:
%
% "Encoding dynamics for multiscale community detection: Markov time sweeping
%  for the Map equation", M.T.Schaub, R. Lambiotte, M. Barahona,
%  Physical Reviews E, Aug, 2012. Vol. 86(2), pp. 026112 ; 
%  see also arXiv:1109.6642
%
% Inputs:	A 			--	Adjacency Matrix of the graph to be analysed
%							(weighted, undirected)
%
%			filename	--	name of the output file
%
%			time		-- 	vector of Markov times to be analysed
%
%           mode        --  'linearised' or 'full' (default);
%                           linearising the dynamics can lead to faster
%                           clustering and will at time t=1 correspond to
%                           the original Infomap. Note that you should not
%                           use the linearised version for times t > 1.
%
%
%           nr_runs     --  number of runs you want to run infomap for each
%                           Markov time. Default: 100

convertAdjMatrixToDirectedPajekGraph(A,[filename '.net']);

if nargin < 5
    nr_runs=100;
end

% OS inputs
seed = randi(10000);
a = num2str(seed); b =[filename '.net']; c = num2str(nr_runs); dir = pwd();
% run MAP equation
system(sprintf('infomap_dir/infomap %s "%s/%s" %s', a, dir, b, c));

% entropy statistics original graph
k = sum(A,2);
D= sparse(diag(k));
pi = k'/sum(k);
M=D\A;
PI = diag(pi);
P0 = PI*M;
P0 = P0(P0 ~= 0);
M = M(M~=0);
P = -P0.*log2(M);
h=sum(P);


if nargin < 4
    mode ='full';
end

% set timeframe and new naming scheme
if nargin<3
    if strcmp(mode,'full')
        time =logspace(-1,2,100);
    else
        time =logspace(-2,0,100);
    end
end
new_name = [filename 'ZoomingMap'];
mkdir(new_name);

h_exp = zeros(1,length(time));
L_exp = h_exp;
clustering_new = zeros(length(A),length(time));
i=1;
for t = time
    
    
    % create new graph
    if strcmp(mode,'full')
        A_new = D*expm(-(eye(size(D))-D\A)*t);
    elseif strcmp(mode,'linearised')
        if(t>1)
            error(['the implemented linearised version is invalid for t>1' ...
                '(you may want to use a different linearisation).'])
        end
        A_new = D*(1-t)+A*t;
    else
        error('Please provide a valid mode')
    end
        
    convertAdjMatrixToDirectedPajekGraph(A_new,[new_name '.net']);
    
    
    % entropy statistics "dynamical" graph
    k = sum(A_new,2);
    D_new= sparse(diag(k));
    pi = k'/sum(k);
    M=D_new\A_new;
    PI = diag(pi);
    P0 = PI*M;
    P0 = P0(P0 ~= 0);
    M = M(M~=0);
    P = -P0.*log2(M);
    h_exp(i)=sum(P);
    
    
    % OS inputs
    seed = randi(10000);
    a = num2str(seed); b =[new_name '.net']; c = num2str(nr_runs);
    
    %run actual infomap program
    system(sprintf('infomap_dir/infomap %s "%s/%s" %s', a, dir, b, c));
    
    temp = load([new_name '.clu']);
    L_exp(i) = temp(1);
    clustering_new(:,i) = temp(2:end);
    % movefile([new_name '.clu'],[new_name '/' new_name num2str(t) '.clu']);
    % movefile([new_name '.tree'],[new_name '/' new_name num2str(t) '.tree']);
    
    i=i+1;
end

N_new= max(clustering_new);

% rename output files and store results
movefile([filename '.net'],[new_name '/' filename '.net']);
temp = load([filename '.clu']);
L = temp(1);
clustering = temp(2:end);
N= max(clustering);
delete([filename '.clu']);
delete([filename '.tree']);
delete([new_name '.clu']);
delete([new_name '.tree']);
delete([new_name '.net']);
% movefile([filename '.clu'],[new_name '/' filename '.clu']);
% movefile([filename '.tree'],[new_name '/' filename '.tree']);
% movefile([new_name '.net'],[new_name '/' new_name '.net']);

save([new_name '/' 'Map_clustering.mat'],'h','h_exp','time','L','L_exp','clustering','clustering_new','N','N_new')

end
