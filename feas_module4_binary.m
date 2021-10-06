% FOURTH MODULE: The design must contact surrounding unit cells 
% (intercell connectivity)
% contactbool = true for a design that contacts adjacent unit cells, and
% false for a design that does not contact adjacent unit cells
function contactbool = feas_module4_binary(CA,sidenum)
    % Initialize values
    SortedCA = sortrows(CA);   
    contactboolx = 0; contactbooly = 0;

    % Initialize node category vectors
    usedpoints = unique(SortedCA);
    lenodes = 2:1:(sidenum-1); 
    renodes = ((sidenum^2)-sidenum+2):1:((sidenum^2)-1);
    tenodes = (2*sidenum):sidenum:((sidenum^2)-sidenum);
    benodes = (sidenum+1):sidenum:((sidenum^2)-(2*sidenum)+1);
    cornernodes = [1,sidenum,((sidenum^2)-sidenum+1),(sidenum^2)];
    lrnodepairs = [lenodes;renodes]';
    tbnodepairs = [benodes;tenodes]';
    
    % Determine if design has at least one repeated external contact point
    for y = 1:1:length(usedpoints)
        otherpoints = usedpoints([1:(y-1),(y+1):end]);
        if ismember(usedpoints(y),cornernodes)
            ocnode = usedpoints(y);
            if ocnode == 1
                lrcnodes = [((sidenum^2)-sidenum+1),(sidenum^2)];
                tbcnodes = [sidenum,(sidenum^2)];
                ocn1 = lrcnodes; ocn2 = tbcnodes; 
                ocn3 = [((sidenum^2)-sidenum+1),sidenum];
            elseif ocnode == sidenum
                lrcnodes = [((sidenum^2)-sidenum+1),(sidenum^2)];
                tbcnodes = [1,((sidenum^2)-sidenum+1)];
                ocn1 = lrcnodes; ocn2 = tbcnodes; 
                ocn3 = [1,(sidenum^2)];
            elseif ocnode == ((sidenum^2)-sidenum+1)
                lrcnodes = [1,sidenum];
                tbcnodes = [sidenum,(sidenum^2)];
                ocn1 = lrcnodes; ocn2 = tbcnodes; 
                ocn3 = [1,(sidenum^2)];
            elseif ocnode == (sidenum^2)
                lrcnodes = [1,sidenum];
                tbcnodes = [1,((sidenum^2)-sidenum+1)];
                ocn1 = lrcnodes; ocn2 = tbcnodes; 
                ocn3 = [((sidenum^2)-sidenum+1),sidenum];
            end
            lrcornercontact = ismember(lrcnodes,otherpoints);
            tbcornercontact = ismember(tbcnodes,otherpoints);
            ocn1contact = ismember(ocn1,otherpoints);
            ocn2contact = ismember(ocn2,otherpoints);
            ocn3contact = ismember(ocn3,otherpoints);
            if all(ocn1contact) || all(ocn2contact) || all(ocn3contact)
                contactboolx = 1; contactbooly = 1;
            elseif any(lrcornercontact)
                contactboolx = 1;
            elseif any(tbcornercontact)
                contactbooly = 1;
            end
        elseif ismember(usedpoints(y),lrnodepairs)
            [i,j] = find(lrnodepairs == usedpoints(y));
            if j == 1
                if ismember(lrnodepairs(i,2),otherpoints)
                    contactboolx = 1;
                end
            elseif j == 2
                if ismember(lrnodepairs(i,1),otherpoints)
                    contactboolx = 1;
                end
            end
        elseif ismember(usedpoints(y),tbnodepairs)
            [i,j] = find(tbnodepairs == usedpoints(y));
            if j == 1
                if ismember(tbnodepairs(i,2),otherpoints)
                    contactbooly = 1;
                end
            elseif j == 2
                if ismember(tbnodepairs(i,1),otherpoints)
                    contactbooly = 1;
                end
            end    
        end
    end
    
    % Determine if intercell connectivity exists in both x and y
    if (contactboolx == 1) && (contactbooly == 1)
        contactbool = 1;
    else
        contactbool = 0;
    end
end