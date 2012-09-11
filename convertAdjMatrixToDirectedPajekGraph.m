% (c) 2012 M Schaub -- michael.schaub09@imperial.ac.uk
function convertAdjMatrixToDirectedPajekGraph(A,filename)
% Converts adjacency matrix of undirected network into pajek graph of an
% directed/undirected network. Assumes network to be fully connected.
% Inputs:
%           A:          Adjacency matrix of undirected graph
%       
%           filename:   string with filename of output file,
%                       e.g. 'pajekgraph.net'. If file already exists
%                       contents are overwritten
% 

nr_nodes = length(A);
% open file and write name list
fid = fopen(filename,'w+','n','ascii');
fprintf(fid,'*Vertices %i \n',nr_nodes);
for i= 1:nr_nodes
    fprintf(fid,'%i \"%i\" \n', int32(i), int32(i));
end


% print arcs
fprintf(fid,'*Arcs\n');


[i j link] = find(A);
% writeout = [j i link];
% dlmwrite(filename,writeout,'-append','delimiter', ' ');

psize = length(link);
for z= 1:psize
            % looping in column first order lower triangular hence swapping indices!
            fprintf(fid,'%i %i %e \n', i(z), j(z), link(z));          

end
fclose('all');


end
