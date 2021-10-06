% FUNCTION FOR CHARACTERISTIC #6 (DENSITY BIAS)
% ASSUMES ODD-VALUED SIDENUM
function [lrdiff,tbdiff] = densityBiasChecker(NC,CA,sidenum,sel,r)
    % Initialize values
    SortedCA = sortrows((sort(CA'))');
    
    % Determine elements belonging to the top, bottom, left, and right
    % halves of the unit cell
    tbdivstart = (sidenum/2)+0.5;
    tbdivend = ((sidenum^2)-tbdivstart+1);
    lrdivstart = (((sidenum^2)-sidenum+1)/2)+0.5;
    lrdivend = ((sidenum^2)-lrdivstart+1);
    
    % Find all interior and exterior edge elements
    lhedgenodes = [(1:1:sidenum),(sidenum:sidenum:lrdivend),...
                   (lrdivstart:1:lrdivend),(1:sidenum:lrdivstart)];
    lhedgenodes = unique(lhedgenodes);
    rhedgenodes = [(lrdivstart:1:lrdivend),(lrdivend:sidenum:(sidenum^2)),...
                   (lrdivstart:sidenum:((sidenum^2)-sidenum+1)),...
                   (((sidenum^2)-sidenum+1):1:(sidenum^2))];
    rhedgenodes = unique(rhedgenodes);
    bhedgenodes = [(1:sidenum:((sidenum^2)-sidenum+1)),(1:1:tbdivstart),...
                   (((sidenum^2)-sidenum+1):1:tbdivend),...
                   (tbdivstart:sidenum:tbdivend)];
    bhedgenodes = unique(bhedgenodes);
    thedgenodes = [(tbdivstart:sidenum:tbdivend),(tbdivstart:1:sidenum),...
                   (tbdivend:1:(sidenum^2)),(sidenum:sidenum:(sidenum^2))];
    thedgenodes = unique(thedgenodes);
    
    plhedgemembers = (ismember(SortedCA(:,1),lhedgenodes) & ismember(SortedCA(:,2),lhedgenodes));
    prhedgemembers = (ismember(SortedCA(:,1),rhedgenodes) & ismember(SortedCA(:,2),rhedgenodes));
    pthedgemembers = (ismember(SortedCA(:,1),thedgenodes) & ismember(SortedCA(:,2),thedgenodes));
    pbhedgemembers = (ismember(SortedCA(:,1),bhedgenodes) & ismember(SortedCA(:,2),bhedgenodes));
    
    lCAedgenodes = SortedCA.*[plhedgemembers,plhedgemembers];
    rCAedgenodes = SortedCA.*[prhedgemembers,prhedgemembers];
    tCAedgenodes = SortedCA.*[pthedgemembers,pthedgemembers];
    bCAedgenodes = SortedCA.*[pbhedgemembers,pbhedgemembers];
    lCAedgenodes = lCAedgenodes(any(lCAedgenodes,2),:);
    rCAedgenodes = rCAedgenodes(any(rCAedgenodes,2),:);
    tCAedgenodes = tCAedgenodes(any(tCAedgenodes,2),:);
    bCAedgenodes = bCAedgenodes(any(bCAedgenodes,2),:);
    
    lx1 = NC(lCAedgenodes(:,1),1); lx2 = NC(lCAedgenodes(:,2),1);
    ly1 = NC(lCAedgenodes(:,1),2); ly2 = NC(lCAedgenodes(:,2),2);
    lL = sqrt(((lx2-lx1).^2)+((ly2-ly1).^2));
    langles = rad2deg(abs(acos((lx2-lx1)./lL)));
    rx1 = NC(rCAedgenodes(:,1),1); rx2 = NC(rCAedgenodes(:,2),1);
    ry1 = NC(rCAedgenodes(:,1),2); ry2 = NC(rCAedgenodes(:,2),2);
    rL = sqrt(((rx2-rx1).^2)+((ry2-ry1).^2));
    rangles = rad2deg(abs(acos((rx2-rx1)./rL)));
    tx1 = NC(tCAedgenodes(:,1),1); tx2 = NC(tCAedgenodes(:,2),1);
    ty1 = NC(tCAedgenodes(:,1),2); ty2 = NC(tCAedgenodes(:,2),2);
    tL = sqrt(((tx2-tx1).^2)+((ty2-ty1).^2));
    tangles = rad2deg(abs(acos((tx2-tx1)./tL)));
    bx1 = NC(bCAedgenodes(:,1),1); bx2 = NC(bCAedgenodes(:,2),1);
    by1 = NC(bCAedgenodes(:,1),2); by2 = NC(bCAedgenodes(:,2),2);
    bL = sqrt(((bx2-bx1).^2)+((by2-by1).^2));
    bangles = rad2deg(abs(acos((bx2-bx1)./bL)));
    
    lCAedgy = [];
    for i = 1:1:size(lCAedgenodes,1)
        x1 = lx1(i); x2 = lx2(i); y1 = ly1(i); y2 = ly2(i); 
        if (langles(i) == 0) && (x1 >= 0) && (x2 <= (0.5*sel))
            if ((y1 == 0) && (y2 == 0)) || ((y1 == sel) && (y2 == sel))
                lCAedgy = [lCAedgy;lCAedgenodes(i,:)];
            end
        elseif (langles(i) == 90) && (y1 >= 0) && (y2 <= sel)
            if ((x1 == 0) && (x2 == 0)) || ((x1 == (0.5*sel)) && (x2 == (0.5*sel)))
                lCAedgy = [lCAedgy;lCAedgenodes(i,:)];
            end
        end
    end
    rCAedgy = [];
    for i = 1:1:size(rCAedgenodes,1)
        x1 = rx1(i); x2 = rx2(i); y1 = ry1(i); y2 = ry2(i);  
        if (rangles(i) == 0)  && (x1 >= (0.5*sel)) && (x2 <= sel)
            if ((y1 == 0) && (y2 == 0)) || ((y1 == sel) && (y2 == sel))
                rCAedgy = [rCAedgy;rCAedgenodes(i,:)];
            end
        elseif (rangles(i) == 90) && (y1 >= 0) && (y2 <= sel)
            if ((x1 == sel) && (x2 == sel)) || ((x1 == (0.5*sel)) && (x2 == (0.5*sel)))
                rCAedgy = [rCAedgy;rCAedgenodes(i,:)];
            end
        end
    end
    tCAedgy = [];
    for i = 1:1:size(tCAedgenodes,1)
        x1 = tx1(i); x2 = tx2(i); y1 = ty1(i); y2 = ty2(i);  
        if (tangles(i) == 0) && (x1 >= 0) && (x2 <= 1)
            if ((y1 == (0.5*sel)) && (y2 == (0.5*sel))) || ((y1 == sel) && (y2 == sel))
                tCAedgy = [tCAedgy;tCAedgenodes(i,:)];
            end
        elseif (tangles(i) == 90) && (y1 >= (0.5*sel)) && (y2 <= sel)
            if ((x1 == sel) && (x2 == sel)) || ((x1 == 0) && (x2 == 0))
                tCAedgy = [tCAedgy;tCAedgenodes(i,:)];
            end
        end
    end
    bCAedgy = [];
    for i = 1:1:size(bCAedgenodes,1)
        x1 = bx1(i); x2 = bx2(i); y1 = by1(i); y2 = by2(i);  
        if (bangles(i) == 0)  && (x1 >= 0) && (x2 <= sel)
            if ((y1 == (0.5*sel)) && (y2 == (0.5*sel))) || ((y1 == 0) && (y2 == 0))
                bCAedgy = [bCAedgy;bCAedgenodes(i,:)];
            end
        elseif (bangles(i) == 90) && (y1 >= 0) && (y2 <= (0.5*sel))
            if ((x1 == sel) && (x2 == sel)) || ((x1 == 0) && (x2 == 0))
                bCAedgy = [bCAedgy;bCAedgenodes(i,:)];
            end
        end
    end
                          
    % Find elements that lie fully within the halves
    lhnodes = 1:1:lrdivend;
    lhCAone = SortedCA(ismember(SortedCA(:,1),lhnodes),:);
    lhCAtwo = SortedCA(ismember(SortedCA(:,2),lhnodes),:);
    lhCA = intersect(lhCAone,lhCAtwo,'rows');
    if ~isempty(lCAedgy)
        lhintmembers = lhCA(~ismember(lhCA,lCAedgy,'rows'),:);
    else
        lhintmembers = lhCA;
    end
    rhnodes = lrdivstart:1:(sidenum^2);
    rhCAone = SortedCA(ismember(SortedCA(:,1),rhnodes),:);
    rhCAtwo = SortedCA(ismember(SortedCA(:,2),rhnodes),:);
    rhCA = intersect(rhCAone,rhCAtwo,'rows');
    if ~isempty(rCAedgy)
        rhintmembers = rhCA(~ismember(rhCA,rCAedgy,'rows'),:);
    else
        rhintmembers = rhCA;
    end
    
    thnodes = []; bhnodes = [];
    for i = 1:1:tbdivstart
        bhnodes = [bhnodes,(i:sidenum:((sidenum^2)-sidenum+i))];
        thnodes = [thnodes,((i+tbdivstart-1):sidenum:...
                   ((sidenum^2)-sidenum+(i+tbdivstart-1)))];
    end
    thCAone = SortedCA(ismember(SortedCA(:,1),thnodes),:);
    thCAtwo = SortedCA(ismember(SortedCA(:,2),thnodes),:);
    thCA = intersect(thCAone,thCAtwo,'rows');
    if ~isempty(tCAedgy)
        thintmembers = thCA(~ismember(thCA,tCAedgy,'rows'),:);
    else
        thintmembers = thCA;
    end
    bhCAone = SortedCA(ismember(SortedCA(:,1),bhnodes),:);
    bhCAtwo = SortedCA(ismember(SortedCA(:,2),bhnodes),:);
    bhCA = intersect(bhCAone,bhCAtwo,'rows');
    if ~isempty(bCAedgy)
        bhintmembers = bhCA(~ismember(bhCA,bCAedgy,'rows'),:);
    else
        bhintmembers = bhCA;
    end
    
    ilx1 = NC(lhintmembers(:,1),1); ilx2 = NC(lhintmembers(:,2),1);
    ily1 = NC(lhintmembers(:,1),2); ily2 = NC(lhintmembers(:,2),2);
    ilL = sqrt(((ilx2-ilx1).^2)+((ily2-ily1).^2));
    irx1 = NC(rhintmembers(:,1),1); irx2 = NC(rhintmembers(:,2),1);
    iry1 = NC(rhintmembers(:,1),2); iry2 = NC(rhintmembers(:,2),2);
    irL = sqrt(((irx2-irx1).^2)+((iry2-iry1).^2));
    itx1 = NC(thintmembers(:,1),1); itx2 = NC(thintmembers(:,2),1);
    ity1 = NC(thintmembers(:,1),2); ity2 = NC(thintmembers(:,2),2);
    itL = sqrt(((itx2-itx1).^2)+((ity2-ity1).^2));
    ibx1 = NC(bhintmembers(:,1),1); ibx2 = NC(bhintmembers(:,2),1);
    iby1 = NC(bhintmembers(:,1),2); iby2 = NC(bhintmembers(:,2),2);
    ibL = sqrt(((ibx2-ibx1).^2)+((iby2-iby1).^2));
    intCA = [lhintmembers;rhintmembers;thintmembers;bhintmembers];
    
    % Find edge members that are in multiple halves, tag their split
    %find all edge nodes
    alledgenodes = [lhedgenodes,rhedgenodes,thedgenodes,bhedgenodes];
    %find all members connecting to/from these nodes
    alledgemembers = (ismember(SortedCA(:,1),alledgenodes) & ismember(SortedCA(:,2),alledgenodes));
    allCAedgenodes = SortedCA.*[alledgemembers,alledgemembers];
    allCAedgenodes = allCAedgenodes(any(allCAedgenodes,2),:);
    %find all edge members by filtering by angle/position
    allx1 = NC(allCAedgenodes(:,1),1); allx2 = NC(allCAedgenodes(:,2),1);
    ally1 = NC(allCAedgenodes(:,1),2); ally2 = NC(allCAedgenodes(:,2),2);
    allL = sqrt(((allx2-allx1).^2)+((ally2-ally1).^2));
    allangles = rad2deg(abs(acos((allx2-allx1)./allL)));
    allCAedgy = [];
    for i = 1:1:size(allCAedgenodes,1)
        x1 = allx1(i); x2 = allx2(i); y1 = ally1(i); y2 = ally2(i);  
        if (allangles(i) == 0) 
            if ((y1 == sel) && (y2 == sel)) || ((y1 == (0.5*sel)) && (y2 == (0.5*sel)))...
                    || ((y1 == 0) && (y2 == 0))
                allCAedgy = [allCAedgy;allCAedgenodes(i,:)];
            end
        elseif (allangles(i) == 90)
            if ((x1 == sel) && (x2 == sel)) || ((x1 == (0.5*sel)) && (x2 == (0.5*sel)))...
                    || ((x1 == 0) && (x2 == 0))
                allCAedgy = [allCAedgy;allCAedgenodes(i,:)];
            end
        end
    end
    %members in this new edge member vector not present in prior edge
    %member vectors are in multiple halves
    properCAedgy = [lCAedgy;rCAedgy;tCAedgy;bCAedgy];
    outstretchedgylog = ~ismember(allCAedgy,properCAedgy,'rows');
    if ~isempty(allCAedgy)
        outstretchededgeCA = allCAedgy.*[outstretchedgylog,outstretchedgylog];
        outstretchededgeCA = outstretchededgeCA(any(outstretchededgeCA,2),:);
        %find these outstretched edge members' locations and tag their area
        %splits accordingly
        edgelrsplit = []; edgetbsplit = []; vertosCAedgy = []; horosCAedgy = [];
        osx1 = NC(outstretchededgeCA(:,1),1); 
        osx2 = NC(outstretchededgeCA(:,2),1);
        osy1 = NC(outstretchededgeCA(:,1),2); 
        osy2 = NC(outstretchededgeCA(:,2),2);
        osL = sqrt(((osx2-osx1).^2)+((osy2-osy1).^2));
        osangles = rad2deg(abs(acos((osx2-osx1)./osL)));
        for i = 1:1:size(outstretchededgeCA,1)
            x1 = osx1(i); x2 = osx2(i); y1 = osy1(i); y2 = osy2(i);
            if (osangles(i) == 0) 
                L1 = (0.5*sel) - x1; L2 = x2 - (0.5*sel);
                lsplitL = ((osL(i)-L1)/osL(i))*osL(i);
                rsplitL = ((osL(i)-L2)/osL(i))*osL(i);
                edgelrsplit = [edgelrsplit;lsplitL,rsplitL];
                horosCAedgy = [horosCAedgy;outstretchededgeCA(i,:)];
            elseif (osangles(i) == 90)
                L1 = (0.5*sel) - y1; L2 = y2 - (0.5*sel);
                bsplitL = ((osL(i)-L1)/osL(i))*osL(i);
                tsplitL = ((osL(i)-L2)/osL(i))*osL(i);
                edgetbsplit = [edgetbsplit;bsplitL,tsplitL];
                vertosCAedgy = [vertosCAedgy;outstretchededgeCA(i,:)];
            end
        end
    else
        edgelrsplit = []; edgetbsplit = [];
    end
    
    % Find interior members that are in multiple halves, tag their split
    %members not part of any prior group are outstretched interior members
    allpriormembers = [allCAedgy;intCA];
    osintmems = ~ismember(SortedCA,allpriormembers,'rows');
    osintCA = SortedCA.*[osintmems,osintmems];
    osintCA = osintCA(any(osintCA,2),:);
    %find these outstretched interior members' locations and tag their area
    %splits accordingly
    lrsplit = []; tbsplit = []; vertosCA = []; horosCA = [];
    iosx1 = NC(osintCA(:,1),1); 
    iosx2 = NC(osintCA(:,2),1);
    iosy1 = NC(osintCA(:,1),2); 
    iosy2 = NC(osintCA(:,2),2);
    iosL = sqrt(((iosx2-iosx1).^2)+((iosy2-iosy1).^2));
    iosangles = rad2deg(abs(acos((iosx2-iosx1)./iosL)));
    for i = 1:1:size(osintCA,1)
        x1 = iosx1(i); x2 = iosx2(i); y1 = iosy1(i); y2 = iosy2(i);
        xL1 = (0.5*sel) - x1; xL2 = x2 - (0.5*sel);
        lsplitL = ((iosL(i)-xL1)/iosL(i))*iosL(i);
        rsplitL = ((iosL(i)-xL2)/iosL(i))*iosL(i);
        lrsplit = [lrsplit;lsplitL,rsplitL];
        horosCA = [horosCA;osintCA(i,:)];
        yL1 = (0.5*sel) - y1; yL2 = y2 - (0.5*sel);
        bsplitL = ((iosL(i)-yL1)/iosL(i))*iosL(i);
        tsplitL = ((iosL(i)-yL2)/iosL(i))*iosL(i);
        tbsplit = [tbsplit;bsplitL,tsplitL];
        vertosCA = [vertosCA;osintCA(i,:)];
    end

    % Determine volume fraction of members in all 4 halves
    %for each half: volume = (pi*r^2)*(length of all interior members) +
    %(0.5*pi*r^2)*(length of all edge members) + (pi*r^2)*(split length of
    %outstretched int members) + (0.5*pi*r^2)*(split length of outstretched
    %edge members)
    A1 = (pi*(r^2)); A2 = (0.5*pi*(r^2));
    if isempty(lrsplit)
        lrval1 = 0; lrval2 = 0;
    else
        lrval1 = (A1*sum(lrsplit(:,1))); lrval2 = (A1*sum(lrsplit(:,2)));
    end
    if isempty(tbsplit)
        tbval1 = 0; tbval2 = 0;
    else
        tbval1 = (A1*sum(tbsplit(:,1))); tbval2 = (A1*sum(tbsplit(:,2)));
    end
    if isempty(edgelrsplit)
        edgelrval1 = 0; edgelrval2 = 0;
    else
        edgelrval1 = (A2*sum(edgelrsplit(:,1))); 
        edgelrval2 = (A2*sum(edgelrsplit(:,2)));
    end
    if isempty(edgetbsplit)
        edgetbval1 = 0; edgetbval2 = 0;
    else
        edgetbval1 = (A2*sum(edgetbsplit(:,1))); 
        edgetbval2 = (A2*sum(edgetbsplit(:,2)));
    end
    lhVol = (A1*sum(ilL))+(A2*sum(lL))+lrval1+edgelrval1;
    rhVol = (A1*sum(irL))+(A2*sum(rL))+lrval2+edgelrval2;
    thVol = (A1*sum(itL))+(A2*sum(tL))+tbval2+edgetbval2;
    bhVol = (A1*sum(ibL))+(A2*sum(bL))+tbval1+edgetbval1;
    solidVol = r*(sel^2);
    lhvf = lhVol/solidVol; rhvf = rhVol/solidVol; 
    thvf = thVol/solidVol; bhvf = bhVol/solidVol;
    
    % Determine the differences in vf
    lrdiff = rhvf - lhvf; tbdiff = thvf - bhvf;
end