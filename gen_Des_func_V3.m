%%% GENERATIVE DESIGN FUNCTION V3
function CA_des_all = gen_Des_func_V3(sidenum,sel,n_des,...
         short_member_prob,short_members_only,numperdes,symm_type)
    
    % Initialize values
    consider_minimum_connections = true; % include the nodal connectivity condition as one of the rules
    continue_finding_more_designs = true; % stop generating designs after a satisfactory design and its children have been found

    rotational_symmetry = false; % toggle generation of rotationally symmetric designs
    vertical_symmetry = false; % toggle generation of vertically symmetric designs
    horizontal_symmetry = false; % toggle generation of horizontally symmetric designs
    full_symmetry = false; % toggle generation of fully symmetric designs (symmetry about central horizontal and vertical axes)
    % Note: Toggling both vertical and horizontal symmetry will not work for this case
    if symm_type == 1
        rotational_symmetry = true;
    elseif symm_type == 2
        vertical_symmetry = true;
    elseif symm_type == 3
        horizontal_symmetry = true;
    elseif symm_type == 4
        full_symmetry = true;
    end
    
    % Parameters for the rules
    NC = generate_NC(sel,sidenum);
    n_total_members = nchoosek(sidenum*sidenum,2);
    CA_all = generate_full_CA(sidenum);

    % Create node labels
    labels = struct;
    for i = 1:sidenum*sidenum
        labels.(strcat('label',num2str(i))) = num2str(i);
    end

    % Generating designs
    CA_des_all = {};
    des_count = 0;

    %all_child_designs_found = false; % boolean indicating that a satisfactory design and its children have been found

    while true
        CA_des_current = zeros(n_total_members,2);
        member_count = 0;
        %{
        % Visualization of design generation
        figure
        % Plot node positions
        for k = 1:size(NC,1)
            plot(NC(k,1),NC(k,2),'*r')
            hold on
            text(NC(k,1),NC(k,2),labels.(strcat('label',num2str(k))),'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15)
            hold on
            drawnow limitrate
        end
        %}
        % Add random member
        viable_member_found = false;
        while (~viable_member_found)
            %rand_member_index = randi(size(CA_all,1));
            %rand_member = CA_all(rand_member_index,:);

            % Choose random member with probability biasing towards shorter members
            first_member = true;
            prob_array = get_probabilities(CA_all,NC,sidenum,sel,short_members_only,short_member_prob,first_member);
            feas_indices = 1:1:size(CA_all,1);
            rand_index = randsample(feas_indices,1,true,prob_array);
            rand_member = CA_all(rand_index,:);

            CA_des_current(member_count+1,:) = rand_member;
            member_count = member_count + 1;

            repeat_member = get_corresponding_repeated_member(rand_member,sidenum);
            if (~isempty(repeat_member))
                CA_des_current(member_count+1,:) = repeat_member;
                member_count = member_count + 1;
            end

            % Add symmetric members (if needed)
            if rotational_symmetry
                rotated_members = get_corresponding_rotated_members(rand_member,NC,sel);
                for i = 1:3
                    CA_des_current(member_count+1,:) = rotated_members(i,:);
                    member_count = member_count + 1;
                    repeat_member = get_corresponding_repeated_member(rotated_members(i,:),sidenum);
                    if (~isempty(repeat_member))
                        CA_des_current(member_count+1,:) = repeat_member;
                        member_count = member_count + 1;
                    end
                end
            end
            if vertical_symmetry
                vert_flipped_member = get_corresponding_vertical_flip_member(rand_member,NC,sel);
                CA_des_current(member_count+1,:) = vert_flipped_member;
                member_count = member_count + 1;
                repeat_member_vert = get_corresponding_repeated_member(vert_flipped_member,sidenum);
                if (~isempty(repeat_member_vert))
                    CA_des_current(member_count+1,:) = repeat_member_vert;
                    member_count = member_count + 1;
                end
            end
            if horizontal_symmetry
                horz_flipped_member = get_corresponding_horizontal_flip_member(rand_member,NC,sel);
                CA_des_current(member_count+1,:) = horz_flipped_member;
                member_count = member_count + 1;
                repeat_member_horz = get_corresponding_repeated_member(horz_flipped_member,sidenum);
                if (~isempty(repeat_member_horz))
                    CA_des_current(member_count+1,:) = repeat_member_horz;
                    member_count = member_count + 1;
                end
            end
            if full_symmetry
                fully_symmetric_members = get_corresponding_fully_symmetric_members(rand_member,NC,sel);
                for i = 1:3
                    CA_des_current(member_count+1,:) = fully_symmetric_members(i,:);
                    member_count = member_count + 1;
                    repeat_member_full = get_corresponding_repeated_member(fully_symmetric_members(i,:),sidenum);
                    if (~isempty(repeat_member_full))
                        CA_des_current(member_count+1,:) = repeat_member_full;
                        member_count = member_count + 1;
                    end
                end
            end

            % Check for intersection and overlap satisfaction
            no_inters = feas_module1_binary(CA_des_current(1:member_count,:),NC,sel);
            no_overlaps = feas_module2_binary(CA_des_current(1:member_count,:),NC,sel);

            if (no_inters && no_overlaps)
                %{
                for i = 1:member_count
                    % Plotting current member
                    x1 = NC(CA_des_current(i,1),1);
                    y1 = NC(CA_des_current(i,1),2);
                    x2 = NC(CA_des_current(i,2),1);
                    y2 = NC(CA_des_current(i,2),2);
                    line = plot([x1,x2],[y1,y2],'-b','LineWidth',2);
                    drawnow limitrate
                    hold on
                end
                %}

                viable_member_found = true;
            else
                % Remove members since they are not feasible
                CA_des_current(1:member_count,:) = zeros(member_count,2);
                member_count = 0;
            end

        end

        % Start generative design process
        random_member_old = rand_member;
        satisfactory_design = false;
        first_member = false;

        while ~satisfactory_design

            % Add feasible connected member at random
            conn_members_rand = get_feasible_connected_members(random_member_old,CA_des_current(1:member_count,:),CA_all,NC,sidenum,sel,short_members_only,rotational_symmetry,vertical_symmetry,horizontal_symmetry,full_symmetry);

            % If no feasible connected member is found, find feasible connected
            % members at node with lowest connections
            if (isempty(conn_members_rand))

                [n_conn_nodes_current, ~] = connectivityCounter(sidenum,CA_des_current(1:member_count,:),NC,sel);

                % Add up counters based on nodal connectivities (sans repeatability)
                %binedges = (1:1:(size(NC,1)+1)) - 0.5;
                %[n_conn_nodes_current,~] = histcounts(CA_des_current(1:member_count,:),binedges);

                unused_nodes = n_conn_nodes_current == 0;
                n_conn_nodes_current(unused_nodes) = 50; 
                [~,min_index] = min(n_conn_nodes_current);
                conn_members_rand = get_feasible_connected_members(min_index,CA_des_current(1:member_count,:),CA_all,NC,sidenum,sel,short_members_only,rotational_symmetry,vertical_symmetry,horizontal_symmetry,full_symmetry);
            end

            iter_member_count = 0;

            %rand_member_index = randi(size(conn_members_rand,1));
            %random_member_new = conn_members_rand(rand_member_index,:);

            % Choose random member with probability biasing towards shorter members
            prob_array = get_probabilities(conn_members_rand,NC,sidenum,sel,short_members_only,short_member_prob,first_member);
            feas_indices = 1:1:size(conn_members_rand,1);
            rand_index = randsample(feas_indices,1,true,prob_array);
            random_member_new = conn_members_rand(rand_index,:);

            CA_des_current(member_count+iter_member_count+1,:) = random_member_new;
            iter_member_count = iter_member_count + 1;

            % Add repeat member (if present)
            repeat_member_new = get_corresponding_repeated_member(random_member_new,sidenum);
            if (~isempty(repeat_member_new))

                CA_des_current(member_count+iter_member_count+1,:) = repeat_member_new;
                iter_member_count = iter_member_count + 1;

            end

            % Add symmetric members and corresponding repeated (if needed and applicable)
            if rotational_symmetry
                rot_members = get_corresponding_rotated_members(random_member_new,NC,sel);
                CA_des_current(member_count+iter_member_count+1:member_count+iter_member_count+size(rot_members,1),:) = rot_members;
                iter_member_count = iter_member_count + size(rot_members,1);

                for i = 1:size(rot_members,1)
                    repeat_member_rot = get_corresponding_repeated_member(rot_members(i,:),sidenum);
                    if (~isempty(repeat_member_rot))
                        if (~ismember(repeat_member_rot,CA_des_current(1:member_count+iter_member_count,:),'rows'))
                            CA_des_current(member_count+iter_member_count+1,:) = repeat_member_rot;
                            iter_member_count = iter_member_count + 1;
                        end
                    end
                end 

            end

            if vertical_symmetry
                flip_member_vert = get_corresponding_vertical_flip_member(random_member_new,NC,sel);
                if (~ismember(flip_member_vert,CA_des_current(1:member_count+iter_member_count,:),'rows'))
                    CA_des_current(member_count+iter_member_count+1:member_count+iter_member_count+size(flip_member_vert,1),:) = flip_member_vert;
                    iter_member_count = iter_member_count + size(flip_member_vert,1);

                    repeat_member_vert = get_corresponding_repeated_member(flip_member_vert,sidenum);
                    if (~isempty(repeat_member_vert))
                        if (~ismember(repeat_member_vert,CA_des_current(1:member_count+iter_member_count,:),'rows'))
                            CA_des_current(member_count+iter_member_count+1,:) = repeat_member_vert;
                            iter_member_count = iter_member_count + 1;
                        end
                    end
                end

            end

            if horizontal_symmetry
                flip_member_horz = get_corresponding_horizontal_flip_member(random_member_new,NC,sel);
                if (~ismember(flip_member_horz,CA_des_current(1:member_count+iter_member_count,:),'rows'))
                    CA_des_current(member_count+iter_member_count+1:member_count+iter_member_count+size(flip_member_horz,1),:) = flip_member_horz;
                    iter_member_count = iter_member_count + size(flip_member_horz,1);

                    repeat_member_horz = get_corresponding_repeated_member(flip_member_horz,sidenum);
                    if (~isempty(repeat_member_horz))
                        if (~ismember(repeat_member_horz,CA_des_current(1:member_count+iter_member_count,:),'rows'))
                            CA_des_current(member_count+iter_member_count+1,:) = repeat_member_horz;
                            iter_member_count = iter_member_count + 1;
                        end
                    end
                end

            end
            if full_symmetry
                fullsym_members = get_corresponding_fully_symmetric_members(random_member_new,NC,sel);

                for i = 1:size(fullsym_members,1)
                    if (~ismember(fullsym_members(i,:),CA_des_current(1:member_count+iter_member_count,:),'rows'))
                        CA_des_current(member_count+iter_member_count+1,:) = fullsym_members(i,:);
                        iter_member_count = iter_member_count + 1;
                        repeat_member_fullsym = get_corresponding_repeated_member(fullsym_members(i,:),sidenum);
                        if (~isempty(repeat_member_fullsym))
                            if (~ismember(repeat_member_fullsym,CA_des_current(1:member_count+iter_member_count,:),'rows'))
                                CA_des_current(member_count+iter_member_count+1,:) = repeat_member_fullsym;
                                iter_member_count = iter_member_count + 1;
                            end
                        end
                    end
                end 

            end

            % Check for intersection and overlap satisfaction
            no_inters = feas_module1_binary(CA_des_current(1:member_count,:),NC,sel);
            no_overlaps = feas_module2_binary(CA_des_current(1:member_count,:),NC,sel);

            if (no_inters && no_overlaps)
                %{
                % Plotting members
                for i = 1:iter_member_count
                    x1 = NC(CA_des_current(member_count+i,1),1);
                    y1 = NC(CA_des_current(member_count+i,1),2);
                    x2 = NC(CA_des_current(member_count+i,2),1);
                    y2 = NC(CA_des_current(member_count+i,2),2);
                    line2_corr = plot([x1,x2],[y1,y2],'-b','LineWidth',2);
                    drawnow limitrate
                    hold on
                end
                %}
            else
                % Remove previous members since they are not feasible
                CA_des_current(member_count+1:member_count+iter_member_count,:) = zeros(iter_member_count,2);
                continue

            end

            % Check for intercell and intracell connectivity and nodal
            % connection rules satisfaction
            intracell_conn = feas_module3_binary(CA_des_current(1:member_count,:));
            intercell_conn = feas_module4_binary(CA_des_current(1:member_count,:),sidenum);

            nodes_rule_sat = true;
            if consider_minimum_connections           
                [n_conn_nodes, n_holes] = connectivityCounter(sidenum,CA_des_current(1:member_count,:),NC,sel);

    %             for j = 1:length(n_conn_nodes)
    %                 if (n_conn_nodes(j) <= 2 && n_conn_nodes(j) > 0)
    %                     nodes_rule_sat = false;
    %                     break
    %                 end
    %             end

                if(any(n_conn_nodes == 1))
                    nodes_rule_sat = false;
                end
            end

            % If the three rules are satisfied, deem the design as satisfactory
            if(intercell_conn && intracell_conn && nodes_rule_sat)
            %if(intercell_conn && intracell_conn)
                satisfactory_design = true;
                disp('Satisfactory design found')
                hold off
            end

            random_member_old = random_member_new;
            member_count = member_count + iter_member_count;
            disp(strcat('Number of members added: ',num2str(member_count)))
        end

        % Add the generated design to the design structure
        CA_des_sorted = sortrows(CA_des_current(1:member_count,:));
        des_count = des_count + 1;
        CA_des_all(des_count) = {CA_des_sorted};

        % Keep adding feasible members till list of feasible members is exhausted
        additional_feas_members = get_all_feasible_connections(CA_des_sorted,CA_all,NC,sidenum,sel,short_members_only,rotational_symmetry,vertical_symmetry,horizontal_symmetry,full_symmetry);
        all_feas_exhausted = false;
        if isempty(additional_feas_members)
            all_feas_exhausted = true;
        end

        CA_des_add = {};
        des_created = 0;
        while ~all_feas_exhausted
            % Select first member from array
            %current_feas_member = additional_feas_members(1,:);

            % Select member at random from array
            prob_array = get_probabilities(additional_feas_members,NC,sidenum,sel,short_members_only,short_member_prob,first_member);
            feas_indices = 1:1:size(additional_feas_members,1);
            rand_index = randsample(feas_indices,1,true,prob_array);
            current_feas_member = additional_feas_members(rand_index,:);

            % Add member to last design from series
            CA_des_new = vertcat(CA_des_sorted,current_feas_member);
            repeat_member = get_corresponding_repeated_member(current_feas_member,sidenum);
            if (~isempty(repeat_member))
                CA_des_new = vertcat(CA_des_new,repeat_member);
            end
            if rotational_symmetry
                rot_members = get_corresponding_rotated_members(current_feas_member,NC,sel);
                CA_des_new = vertcat(CA_des_new,rot_members);
                for i = 1:size(rot_members,1)
                    repeat_rot_member = get_corresponding_repeated_member(rot_members(i,:),sidenum);
                    if (~isempty(repeat_rot_member))
                        if (~ismember(repeat_rot_member,CA_des_new,'rows'))
                            CA_des_new = vertcat(CA_des_new,repeat_rot_member);
                        end
                    end
                end
            end
            if vertical_symmetry
                flip_member_vert = get_corresponding_vertical_flip_member(current_feas_member,NC,sel);
                if (~ismember(flip_member_vert,CA_des_new,'rows'))
                    CA_des_new = vertcat(CA_des_new,flip_member_vert);
                    repeat_vert_member = get_corresponding_repeated_member(flip_member_vert,sidenum);
                    if (~isempty(repeat_vert_member))
                        if (~ismember(repeat_vert_member,CA_des_new,'rows'))
                            CA_des_new = vertcat(CA_des_new,repeat_vert_member);
                        end
                    end
                end
            end
            if horizontal_symmetry
                flip_member_horz = get_corresponding_horizontal_flip_member(current_feas_member,NC,sel);
                if (~ismember(flip_member_horz,CA_des_new,'rows'))
                    CA_des_new = vertcat(CA_des_new,flip_member_horz);
                    repeat_horz_member = get_corresponding_repeated_member(flip_member_horz,sidenum);
                    if (~isempty(repeat_horz_member))
                        if (~ismember(repeat_horz_member,CA_des_new,'rows'))
                            CA_des_new = vertcat(CA_des_new,repeat_horz_member);
                        end
                    end
                end
            end
            if full_symmetry
                fullsym_members = get_corresponding_fully_symmetric_members(current_feas_member,NC,sel);
                for i = 1:size(fullsym_members,1)
                    if (~ismember(fullsym_members(i,:),CA_des_new,'rows'))
                        CA_des_new = vertcat(CA_des_new,fullsym_members(i,:));
                        repeat_fullsym_member = get_corresponding_repeated_member(fullsym_members(i,:),sidenum);
                        if (~isempty(repeat_fullsym_member))
                            if (~ismember(repeat_fullsym_member,CA_des_new,'rows'))
                                CA_des_new = vertcat(CA_des_new,repeat_fullsym_member);
                            end
                        end
                    end
                end
            end
            %visualize_CA(CA_des_new,NC,labels)

            % Save design
            des_created = des_created + 1;
            CA_des_sorted = sortrows(CA_des_new);
            CA_des_add(des_created) = {CA_des_sorted};

            disp('Additional candidate design created')
            additional_feas_members = get_all_feasible_connections(CA_des_sorted,CA_all,NC,sidenum,sel,short_members_only,rotational_symmetry,vertical_symmetry,horizontal_symmetry,full_symmetry);

            % Complete loop if all feasible members are used up, and select
            % some of the new designs to save permanently
            if isempty(additional_feas_members)
                all_feas_exhausted = true;
                if des_created < numperdes
                    for j = 1:1:des_created
                        des_count = des_count + 1;
                        CA_des_all(des_count) = CA_des_add(des_created);
                    end
                else
                    idxbin = randperm(des_created,numperdes);
                    for k = 1:1:length(idxbin)
                        des_count = des_count + 1;
                        CA_des_all(des_count) = CA_des_add(idxbin(k));
                    end
                end
            end
        end
        %all_child_designs_found = true;
        %if (all_child_designs_found && ~continue_finding_more_designs)
        if ~continue_finding_more_designs
            break
        end
        if (des_count >= n_des)
            break
        end
    end
end

%%% REQUISITE SUBFUNCTIONS

function NC = generate_NC(sel,sidenum)
    notchvec = linspace(0,1,sidenum);
    NC = [];
    for i = 1:1:sidenum
        for j = 1:1:sidenum
            NC = [NC;notchvec(i),notchvec(j)];
        end
    end
    NC = sel.*NC;
end

function CA_full = generate_full_CA(sidenum)
    member_count = 0;
    n_total_members = nchoosek(sidenum*sidenum,2);
    CA_full = zeros(n_total_members,2);
    for i = 1:((sidenum*sidenum)-1)
        for j = i+1:(sidenum*sidenum)
            CA_full(member_count+1,:) = [i,j];
            member_count = member_count + 1;
        end
    end
end

function repeat_member = get_corresponding_repeated_member(member,sidenum)
    repeat_member = [];
    left_edge_nodes = 1:1:sidenum;
    bottom_edge_nodes = 1:sidenum:(sidenum*sidenum - sidenum)+1;
    top_edge_nodes = sidenum:sidenum:sidenum*sidenum+1;
    right_edge_nodes = (sidenum*sidenum - sidenum)+1:1:sidenum*sidenum;
    
    if (all(ismember(member,left_edge_nodes))) % member is a left edge member
        repeat_member = member + (sidenum*sidenum - sidenum);
    elseif (all(ismember(member,bottom_edge_nodes))) % member is a bottom edge member
        repeat_member = member + sidenum - 1;
    elseif (all(ismember(member,right_edge_nodes))) % member is a right edge member
        repeat_member = member - (sidenum*sidenum - sidenum);
    elseif (all(ismember(member,top_edge_nodes))) % member is a top edge member
        repeat_member = member - (sidenum - 1);
    end
end

function prob_array = get_probabilities(feas_members,NC,sidenum,sel,only_short_members,prob_short_member,first_member_addition)
    prob_array = zeros(1,size(feas_members,1));
    
    if (only_short_members && ~first_member_addition)
        prob_array = 1/size(feas_members,1).*ones(1,size(feas_members,1));
    elseif (~only_short_members || first_member_addition)
        short_member_indices = zeros(1,size(feas_members,1));
        long_member_indices = zeros(1,size(feas_members,1));
        n_short_members = 0;
        n_long_members = 0;
        for i = 1:size(feas_members,1)
            current_member = feas_members(i,:);
            node1_pos = [NC(current_member(1),1), NC(current_member(1),2)];
            node2_pos = [NC(current_member(2),1), NC(current_member(2),2)];
            member_length = sqrt((node1_pos(1) - node2_pos(1))^2 + (node1_pos(2) - node2_pos(2))^2);
            if member_length <= (sqrt(2)+1e-5)*(sel/(sidenum - 1)) % additional 1e-5 to resolve precision issue
                n_short_members = n_short_members + 1;
                short_member_indices(n_short_members) = i;
            else
                n_long_members = n_long_members + 1;
                long_member_indices(n_long_members) = i;
            end
        end
        short_member_indices_red = short_member_indices(1:n_short_members);
        long_member_indices_red = long_member_indices(1:n_long_members);
        
        if only_short_members
            prob_array(short_member_indices_red) = 1/n_short_members;
            prob_array(long_member_indices_red) = 0;
        else
            prob_array(short_member_indices_red) = prob_short_member/n_short_members;
            prob_array(long_member_indices_red) = (1 - prob_short_member)/n_long_members;
        end
    end
    
end

function rot_members_red = get_corresponding_rotated_members(member,nodal_conn,sel)
    rot_members_red = zeros(3,2);
    x1 = nodal_conn(member(1,1),1);
    y1 = nodal_conn(member(1,1),2);
    x2 = nodal_conn(member(1,2),1);
    y2 = nodal_conn(member(1,2),2);
    % 3 rotated members for 5x5 node grid corresponding to "member" to be
    % found
    for i = 1:1:3
        ry1 = round(4*x1,2)/4; ry2 = round(4*x2,2)/4;
        rx1 = round(4*(-y1 + (sel)),2)/4;
        rx2 = round(4*(-y2 + (sel)),2)/4;
        Qx = find(nodal_conn(:,1)==rx1); 
        Qy = find(nodal_conn(:,2)==ry1);
        Q = intersect(Qx,Qy);
        Dx = find(nodal_conn(:,1)==rx2); 
        Dy = find(nodal_conn(:,2)==ry2);
        D = intersect(Dx,Dy);
        rot_members_red(i,:) = [min(Q,D),max(Q,D)];
        x1 = rx1; 
        x2 = rx2; 
        y1 = ry1; 
        y2 = ry2;
    end
end

function flip_member = get_corresponding_horizontal_flip_member(member,nodal_conn,sel)
    x1 = nodal_conn(member(1,1),1);
    y1 = nodal_conn(member(1,1),2);
    x2 = nodal_conn(member(1,2),1);
    y2 = nodal_conn(member(1,2),2);
    
    % Determine horizontal flipped member 
    ry1 = y1;
    ry2 = y2;
    %rx1 = round(4*(x1 + (sel)),2)/4;
    %rx2 = round(4*(x2 + (sel)),2)/4;
    if x1 == sel/2
        rx1 = x1;
    else
        rx1 = sel - x1;
    end
    if x2 == sel/2
        rx2 = x2;
    else
        rx2 = sel - x2;
    end
    %Qx = find(nodal_conn(:,1)==rx1);
    %Qy = find(nodal_conn(:,2)==ry1);
    Qx = find_min_distance_indices(nodal_conn(:,1),rx1);
    Qy = find_min_distance_indices(nodal_conn(:,2),ry1);
    Q = intersect(Qx,Qy);
    %Dx = find(nodal_conn(:,1)==rx2);
    %Dy = find(nodal_conn(:,2)==ry2);
    Dx = find_min_distance_indices(nodal_conn(:,1),rx2);
    Dy = find_min_distance_indices(nodal_conn(:,2),ry2);
    D = intersect(Dx,Dy);
    flip_member = [min(Q,D),max(Q,D)];

end

function flip_member = get_corresponding_vertical_flip_member(member,nodal_conn,sel)
    x1 = nodal_conn(member(1,1),1);
    y1 = nodal_conn(member(1,1),2);
    x2 = nodal_conn(member(1,2),1);
    y2 = nodal_conn(member(1,2),2);
    
    % Determine vertical flipped member 
    %ry1 = round(4*(y1 + (sel)),2)/4;
    %ry2 = round(4*(y2 + (sel)),2)/4;
    if y1 == sel/2
        ry1 = y1;
    else
        ry1 = sel - y1;
    end
    if y2 == sel/2
        ry2 = y2;
    else
        ry2 = sel - y2;
    end
    rx1 = x1;
    rx2 = x2;
    %Qx = find(nodal_conn(:,1)==rx1);
    %Qy = find(nodal_conn(:,2)==ry1);
    Qx = find_min_distance_indices(nodal_conn(:,1),rx1);
    Qy = find_min_distance_indices(nodal_conn(:,2),ry1);
    Q = intersect(Qx,Qy);
    %Dx = find(nodal_conn(:,1)==rx2);
    %Dy = find(nodal_conn(:,2)==ry2);
    Dx = find_min_distance_indices(nodal_conn(:,1),rx2);
    Dy = find_min_distance_indices(nodal_conn(:,2),ry2);
    D = intersect(Dx,Dy);
    flip_member = [min(Q,D),max(Q,D)];
end

function flipped_members = get_corresponding_fully_symmetric_members(member,nodal_conn,sel)
    flipped_members = zeros(3,2);
    
    % First find horizontal flipped member
    flip_member_horz = get_corresponding_horizontal_flip_member(member,nodal_conn,sel);
    flipped_members(1,:) = flip_member_horz;
    
    % Then, find vertically flipped members for both original and horizontally flipped member
    flip_member_vert_orig = get_corresponding_vertical_flip_member(member,nodal_conn,sel);
    flip_member_vert_horz = get_corresponding_vertical_flip_member(flip_member_horz,nodal_conn,sel);
    
    flipped_members(2,:) = flip_member_vert_orig;
    flipped_members(3,:) = flip_member_vert_horz;
   
end

function conn_members_red = get_feasible_connected_members(member,CA_des,CA_full,nodal_conn,sidenum,sel,only_short_members,rot_symm,vert_symm,horz_symm,full_symm)
    % CA_des is the current CA without the connections
    feas_conn_members = zeros(size(CA_full,1),2);
    n_feas_conn_members = 0;
    for i = 1:size(CA_full,1)
        if any(ismember(member,CA_full(i,:)))
            current_member = CA_full(i,:);
            
            % Determine symmetric members (if needed)
            symmetric_members = zeros(5,2);
            num_symmetric_members = 0;
            if rot_symm
                rot_members = get_corresponding_rotated_members(current_member,nodal_conn,sel);
                symmetric_members(num_symmetric_members+1:num_symmetric_members+size(rot_members,1),:) = rot_members;
                num_symmetric_members = num_symmetric_members + size(rot_members,1);
            end
            if vert_symm
                vert_flip_members = get_corresponding_vertical_flip_member(current_member,nodal_conn,sel);
                symmetric_members(num_symmetric_members+1:num_symmetric_members+size(vert_flip_members,1),:) = vert_flip_members;
                num_symmetric_members = num_symmetric_members + 1;
            end
            if horz_symm
                horz_flip_members = get_corresponding_horizontal_flip_member(current_member,nodal_conn,sel);
                symmetric_members(num_symmetric_members+1:num_symmetric_members+size(horz_flip_members,1),:) = horz_flip_members;
                num_symmetric_members = num_symmetric_members + 1;
            end    
            if full_symm
                full_flip_members = get_corresponding_fully_symmetric_members(current_member,nodal_conn,sel);
                symmetric_members(num_symmetric_members+1:num_symmetric_members+size(full_flip_members,1),:) = full_flip_members;
                num_symmetric_members = num_symmetric_members + size(full_flip_members,1);
            end
            symmetric_members_red = symmetric_members(1:num_symmetric_members,:);
            symmetric_members_comb = vertcat(current_member,symmetric_members_red);
            symmetric_members_comb = unique(symmetric_members_comb,'rows');
            
            CA_des_current = vertcat(CA_des,symmetric_members_comb);
            
            % Check for intersection and overlap satisfaction
            no_inters = feas_module1_binary(CA_des_current,nodal_conn,sel);
            no_overlaps = feas_module2_binary(CA_des_current,nodal_conn,sel);
        
            if (no_inters && no_overlaps)
                if only_short_members
                    node1_pos = [nodal_conn(current_member(1),1), nodal_conn(current_member(1),2)];
                    node2_pos = [nodal_conn(current_member(2),1), nodal_conn(current_member(2),2)];
                    member_length = sqrt((node1_pos(1) - node2_pos(1))^2 + (node1_pos(2) - node2_pos(2))^2);
                    if member_length <= (sqrt(2)+1e-5)*(sel/(sidenum - 1)) % additional 1e-5 to resolve precision issue
                        feas_conn_members(n_feas_conn_members+1,:) = current_member;
                        n_feas_conn_members = n_feas_conn_members + 1;
                    end
                else
                    feas_conn_members(n_feas_conn_members+1,:) = current_member;
                    n_feas_conn_members = n_feas_conn_members + 1;
                end
            end
        end
    end
    conn_members_red = feas_conn_members(1:n_feas_conn_members,:);
end

function all_feasible_members_un = get_all_feasible_connections(CA_sat,CA_complete,nodal_conn,sidenum,sel,only_short_members,rotational_symmetry,vertical_symmetry,horizontal_symmetry,full_symmetry)
    all_feasible_members = zeros(size(CA_complete,1),2);
    n_feas_conns = 0;
    for i = 1:size(CA_sat,1)
        current_member = CA_sat(i,:);
        current_conn_members = get_feasible_connected_members(current_member,CA_sat,CA_complete,nodal_conn,sidenum,sel,only_short_members,rotational_symmetry,vertical_symmetry,horizontal_symmetry,full_symmetry);
        all_feasible_members(n_feas_conns+1:n_feas_conns+size(current_conn_members,1),:) = current_conn_members; 
        n_feas_conns = n_feas_conns + size(current_conn_members,1);
    end
    all_feasible_members_red = all_feasible_members(1:n_feas_conns,:);
    all_feasible_members_un = unique(all_feasible_members_red,'rows');
end

% Function creates a figure visualizing the corresponding CA
function [] = visualize_CA(CA_current,NC,labels_struct)
    figure
    % Plot node positions
    for k = 1:size(NC,1)
        plot(NC(k,1),NC(k,2),'*r')
        hold on
        text(NC(k,1),NC(k,2),labels_struct.(strcat('label',num2str(k))),'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15)
        hold on
        drawnow limitrate
    end
    % Plot the members of CA_current
    for m = 1:size(CA_current,1)
        x1 = NC(CA_current(m,1),1);
        y1 = NC(CA_current(m,1),2);
        x2 = NC(CA_current(m,2),1);
        y2 = NC(CA_current(m,2),2);
        plot([x1,x2],[y1,y2],'-b','LineWidth',2);
        drawnow limitrate
        hold on
    end
    hold off
end

% Function finds the indices of the array that have closest value to given value
% (substitute to directly using find() function)
function min_indices = find_min_distance_indices(array,val)
    diff_array = abs(array - val);
    min_diff = min(diff_array);
    min_indices = find(diff_array == min_diff);
end
    