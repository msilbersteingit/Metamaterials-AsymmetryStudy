% THIRD MODULE: The entire design must be connected within itself
% (intracell connecvitiy)
% flowbool = true for a unit cell that's fully connected within itself, and
% false for a unit cell that is not
function flowbool = feas_module3_binary(CA)
    % Initialize values
    SortedCA = sortrows(CA);
    flowbool = 0;

    % Find unique startpoints and endpoints
    startpoints = unique(SortedCA(:,1)); endpoints = unique(SortedCA(:,2));
    startunique = setdiff(startpoints,endpoints);
    endunique = setdiff(endpoints,startpoints);
    
    % Check if there are any unique start- or endpoints at all
    if isempty(startunique) || isempty(endunique)
        flowbool = 1;
        return
    else
        % For each unique startpoint, determine if connectivity exists to
        % all unique endpoints
        endreached = zeros(1,length(endunique));
        for v = 1:1:length(startunique)
            mCA = SortedCA(SortedCA(:,1)==startunique(v),:);
            while size(mCA,1) < size(SortedCA,1)
                nmCA = mCA;
                noneleft = 1;
                for z = 1:1:size(mCA,1)
                    tCA = SortedCA(SortedCA(:,1)== (mCA(z,2)),:);
                    repeatrows = [];
                    for q = 1:1:size(tCA,1)
                        if ismember(tCA(q,:),mCA,'rows')
                            repeatrows = [repeatrows,q];
                        end
                    end
                    tCA(repeatrows,:) = [];
                    if ~isempty(tCA)
                        noneleft = 0;
                    end
                    nmCA = [nmCA;tCA];
                end
                nmCA = unique(nmCA,'rows');
                flowdet = ismember(endunique,nmCA(:,2));
                % If the current unique startpoint leads to all unique
                % endpoints, the design is intraconnected
                if all(flowdet) == true
                    flowbool = 1;
                    return
                % If this unique startpoint only connects to one unique
                % endpoint, this portion of the design could be isolated
                % from the rest, depending on whether any other flowpaths
                % terminate in the current one
                elseif (noneleft == 1) && (sum(flowdet) == 1)
                    usednodes = unique(nmCA);
                    extflowdet = ismember(usednodes,SortedCA(:,2));
                    if ~any(extflowdet)
                        flowbool = 0;
                        disp('here');
                        return
                    else
                        endreached(flowdet==1) = endreached(flowdet==1)+1;
                        break
                    end
                % If there are no longer any members in the possible
                % flowpaths, this unique startpoint doesn't connect to all
                % unique endpoints  
                elseif noneleft == 1
                    endreached(flowdet==1) = endreached(flowdet==1) + 1;
                    break
                end
                mCA = nmCA;
            end
        end
        %{
        % If a single startpoint doesn't lead to all endpoints, determine
        % if all startpoints cumulatively lead to all endpoints
        if all(endreached) && any(endreached>1)
            flowbool = 1;
        end
        %}
    end
end