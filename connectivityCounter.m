% Connectivity Counter Function
% This function counts the number of elements connecting to each node,
% accounting for unit cell repeatability.  IT also counts the number of
% holes (unused nodes), accounting for repeatability
function [N,numHoles] = connectivityCounter(sidenum,CA,NC,sel)
    % Initialize variables
    ND = NC./sel;
    
    % Add up counters based on nodal connectivities (sans repeatability)
    binedges = (1:1:(size(NC,1)+1)) - 0.5;
    [N,~] = histcounts(CA,binedges);
    
    % Determine number of holes, accounting for repeatability
    leftedgenodes = 2:1:(sidenum-1); 
    bottomedgenodes = (sidenum+1):sidenum:((sidenum^2)-(2*sidenum)+1);
    edgenodes = [leftedgenodes,bottomedgenodes];
    numHoles = nnz(~N);
    if N(1) == 0
        numHoles = numHoles - 3;
    end
    for j = 1:1:length(edgenodes)
        if N(edgenodes(j)) == 0
            numHoles = numHoles - 1;
        end
    end
    
    % Loop through each node to account for repeatability in connectivity
    for i = 1:1:size(NC,1)
        if (ND(i,1) == 0) && (ND(i,2) == 0) % Node in bottom left corner
            % Identify opposite node
            oppnode = i+((sidenum^2)-1);
            
            % Find elements connecting to/from oppnode
            indone = CA(:,1) == oppnode; indtwo = CA(:,2) == oppnode;
            onCAone = CA(indone,:); onCAtwo = CA(indtwo,:);
            onCA = [setdiff(onCAone,onCAtwo,'rows');...
                    setdiff(onCAtwo,onCAone,'rows')];
                
            % Identify adjacent nodes
            badjnode = i+(sidenum-1);
            aadjnode = i+(sidenum*(sidenum-1));
            
            % Find elements connecting to/from adjacent nodes
            bindone = CA(:,1) == badjnode; bindtwo = CA(:,2) == badjnode;
            bnCAone = CA(bindone,:); bnCAtwo = CA(bindtwo,:);
            bnCA = [setdiff(bnCAone,bnCAtwo,'rows');...
                    setdiff(bnCAtwo,bnCAone,'rows')];
            aindone = CA(:,1) == aadjnode; aindtwo = CA(:,2) == aadjnode;
            anCAone = CA(aindone,:); anCAtwo = CA(aindtwo,:);
            anCA = [setdiff(anCAone,anCAtwo,'rows');...
                    setdiff(anCAtwo,anCAone,'rows')];
                
            % Based on location of node i relative to adjacent nodes, 
            % eliminate duplicate elements
            for q = 1:1:size(bnCA,1)
                if (ND(bnCA(q,1),1) ~= ND(bnCA(q,2),1)) && ...
                   (ND(bnCA(q,1),2) ~= ND(bnCA(q,2),2))
                    N(i) = N(i) + 1;
                end
            end
            for q = 1:1:size(anCA,1)
                if (ND(anCA(q,1),1) ~= ND(anCA(q,2),1)) && ...
                   (ND(anCA(q,1),2) ~= ND(anCA(q,2),2))
                    N(i) = N(i) + 1;
                end
            end
            
            % Add connections from opposite node
            N(i) = N(i) + size(onCA,1); 
            
        elseif (ND(i,1) == 0) && (ND(i,2) == 1) % Node in top left corner
            % Identify opposite node
            oppnode = i+((sidenum-1)^2); 
            
            % Find elements connecting to/from oppnode
            indone = CA(:,1) == oppnode; indtwo = CA(:,2) == oppnode;
            onCAone = CA(indone,:); onCAtwo = CA(indtwo,:);
            onCA = [setdiff(onCAone,onCAtwo,'rows');...
                    setdiff(onCAtwo,onCAone,'rows')];
                
            % Identify adjacent nodes
            badjnode = i-(sidenum-1);
            aadjnode = i+(sidenum*(sidenum-1));
            
            % Find elements connecting to/from adjacent nodes
            bindone = CA(:,1) == badjnode; bindtwo = CA(:,2) == badjnode;
            bnCAone = CA(bindone,:); bnCAtwo = CA(bindtwo,:);
            bnCA = [setdiff(bnCAone,bnCAtwo,'rows');...
                    setdiff(bnCAtwo,bnCAone,'rows')];
            aindone = CA(:,1) == aadjnode; aindtwo = CA(:,2) == aadjnode;
            anCAone = CA(aindone,:); anCAtwo = CA(aindtwo,:);
            anCA = [setdiff(anCAone,anCAtwo,'rows');...
                    setdiff(anCAtwo,anCAone,'rows')];
                
            % Based on location of node i relative to adjacent nodes, 
            % eliminate duplicate elements
            for q = 1:1:size(bnCA,1)
                if (ND(bnCA(q,1),1) ~= ND(bnCA(q,2),1)) && ...
                   (ND(bnCA(q,1),2) ~= ND(bnCA(q,2),2))
                    N(i) = N(i) + 1;
                end
            end
            for q = 1:1:size(anCA,1)
                if (ND(anCA(q,1),1) ~= ND(anCA(q,2),1)) && ...
                   (ND(anCA(q,1),2) ~= ND(anCA(q,2),2))
                    N(i) = N(i) + 1;
                end
            end
            
            % Add connections from opposite node
            N(i) = N(i) + size(onCA,1); 
            
        elseif (ND(i,1) == 1) && (ND(i,2) == 0) % Node in bot. right corner
            % Identify opposite node
            oppnode = i-((sidenum-1)^2); 
            
            % Find elements connecting to/from oppnode
            indone = CA(:,1) == oppnode; indtwo = CA(:,2) == oppnode;
            onCAone = CA(indone,:); onCAtwo = CA(indtwo,:);
            onCA = [setdiff(onCAone,onCAtwo,'rows');...
                    setdiff(onCAtwo,onCAone,'rows')];
            
            % Identify adjacent nodes
            badjnode = i+(sidenum-1);
            aadjnode = i-(sidenum*(sidenum-1));
            
            % Find elements connecting to/from adjacent nodes
            bindone = CA(:,1) == badjnode; bindtwo = CA(:,2) == badjnode;
            bnCAone = CA(bindone,:); bnCAtwo = CA(bindtwo,:);
            bnCA = [setdiff(bnCAone,bnCAtwo,'rows');...
                    setdiff(bnCAtwo,bnCAone,'rows')];
            aindone = CA(:,1) == aadjnode; aindtwo = CA(:,2) == aadjnode;
            anCAone = CA(aindone,:); anCAtwo = CA(aindtwo,:);
            anCA = [setdiff(anCAone,anCAtwo,'rows');...
                    setdiff(anCAtwo,anCAone,'rows')];
                
            % Based on location of node i relative to adjacent nodes, 
            % eliminate duplicate elements
            for q = 1:1:size(bnCA,1)
                if (ND(bnCA(q,1),1) ~= ND(bnCA(q,2),1)) && ...
                   (ND(bnCA(q,1),2) ~= ND(bnCA(q,2),2))
                    N(i) = N(i) + 1;
                end
            end
            for q = 1:1:size(anCA,1)
                if (ND(anCA(q,1),1) ~= ND(anCA(q,2),1)) && ...
                   (ND(anCA(q,1),2) ~= ND(anCA(q,2),2))
                    N(i) = N(i) + 1;
                end
            end
            
            % Add connections from opposite node
            N(i) = N(i) + size(onCA,1); 
            
        elseif (ND(i,1) == 1) && (ND(i,2) == 1) % Node in top right corner
            % Identify opposite node
            oppnode = i-((sidenum^2)-1);
            
            % Find elements connecting to/from oppnode
            indone = CA(:,1) == oppnode; indtwo = CA(:,2) == oppnode;
            onCAone = CA(indone,:); onCAtwo = CA(indtwo,:);
            onCA = [setdiff(onCAone,onCAtwo,'rows');...
                    setdiff(onCAtwo,onCAone,'rows')];
            
            % Identify adjacent nodes
            badjnode = i-(sidenum-1);
            aadjnode = i-(sidenum*(sidenum-1));
            
            % Find elements connecting to/from adjacent nodes
            bindone = CA(:,1) == badjnode; bindtwo = CA(:,2) == badjnode;
            bnCAone = CA(bindone,:); bnCAtwo = CA(bindtwo,:);
            bnCA = [setdiff(bnCAone,bnCAtwo,'rows');...
                    setdiff(bnCAtwo,bnCAone,'rows')];
            aindone = CA(:,1) == aadjnode; aindtwo = CA(:,2) == aadjnode;
            anCAone = CA(aindone,:); anCAtwo = CA(aindtwo,:);
            anCA = [setdiff(anCAone,anCAtwo,'rows');...
                    setdiff(anCAtwo,anCAone,'rows')];
                
            % Based on location of node i relative to adjacent nodes, 
            % eliminate duplicate elements
            for q = 1:1:size(bnCA,1)
                if (ND(bnCA(q,1),1) ~= ND(bnCA(q,2),1)) && ...
                   (ND(bnCA(q,1),2) ~= ND(bnCA(q,2),2))
                    N(i) = N(i) + 1;
                end
            end
            for q = 1:1:size(anCA,1)
                if (ND(anCA(q,1),1) ~= ND(anCA(q,2),1)) && ...
                   (ND(anCA(q,1),2) ~= ND(anCA(q,2),2))
                    N(i) = N(i) + 1;
                end
            end
            
            % Add connections from opposite node
            N(i) = N(i) + size(onCA,1);  
            
        elseif (ND(i,1) == 0) % Edge node is on left side
            % Identify opposite node
            oppnode = i+(sidenum*(sidenum-1));
            
            % Find elements connecting to/from oppnode
            indone = CA(:,1) == oppnode; indtwo = CA(:,2) == oppnode;
            onCAone = CA(indone,:); onCAtwo = CA(indtwo,:);
            onCA = [setdiff(onCAone,onCAtwo,'rows');...
                    setdiff(onCAtwo,onCAone,'rows')];
        
            % Based on location of node i/oppnode, eliminate duplicate
            % elements
            for q = 1:1:size(onCA,1)
                if ND(onCA(q,1),1) ~= ND(onCA(q,2),1)
                    N(i) = N(i) + 1;
                end
            end
        elseif (ND(i,1) == 1) % Edge node is on right side
            % Identify opposite node
            oppnode = i-(sidenum*(sidenum-1));
            
            % Find elements connecting to/from oppnode
            indone = CA(:,1) == oppnode; indtwo = CA(:,2) == oppnode;
            onCAone = CA(indone,:); onCAtwo = CA(indtwo,:);
            onCA = [setdiff(onCAone,onCAtwo,'rows');...
                    setdiff(onCAtwo,onCAone,'rows')];
        
            % Based on location of node i/oppnode, eliminate duplicate
            % elements
            for q = 1:1:size(onCA,1)
                if ND(onCA(q,1),1) ~= ND(onCA(q,2),1)
                    N(i) = N(i) + 1;
                end
            end
        elseif (ND(i,2) == 0) % Edge node is on bottom side
            % Identify opposite node
            oppnode = i+(sidenum-1);
            
            % Find elements connecting to/from oppnode
            indone = CA(:,1) == oppnode; indtwo = CA(:,2) == oppnode;
            onCAone = CA(indone,:); onCAtwo = CA(indtwo,:);
            onCA = [setdiff(onCAone,onCAtwo,'rows');...
                    setdiff(onCAtwo,onCAone,'rows')];
        
            % Based on location of node i/oppnode, eliminate duplicate
            % elements
            for q = 1:1:size(onCA,1)
                if ND(onCA(q,1),2) ~= ND(onCA(q,2),2)
                    N(i) = N(i) + 1;
                end
            end
        elseif (ND(i,2) == 1) % Edge node is on top side
            % Identify opposite node
            oppnode = i-(sidenum-1);
            
            % Find elements connecting to/from oppnode
            indone = CA(:,1) == oppnode; indtwo = CA(:,2) == oppnode;
            onCAone = CA(indone,:); onCAtwo = CA(indtwo,:);
            onCA = [setdiff(onCAone,onCAtwo,'rows');...
                    setdiff(onCAtwo,onCAone,'rows')];
        
            % Based on location of node i/oppnode, eliminate duplicate
            % elements
            for q = 1:1:size(onCA,1)
                if ND(onCA(q,1),2) ~= ND(onCA(q,2),2)
                    N(i) = N(i) + 1;
                end
            end
        end
    end
end