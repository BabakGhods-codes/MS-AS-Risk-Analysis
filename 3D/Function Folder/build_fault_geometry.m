function FG = build_fault_geometry( ...
    Fault_Type, W, L, Hf, Dip, dL, ...
    Xs, Ys, Site1_loc)

%%% Distance Calculation
%%% Discretization of W and L and Calculation 
W_vec = 0: dL : W; 
L_vec = 0: dL : L;

%%% Calcualtion of the Grid points Distance to the Site
X_poits  = Xs + L_vec';
Y_points = Ys + W_vec'*cos(Dip);  
Y_points_repmat = repmat(Y_points,1,length(L_vec));

Z_points = -Hf - W_vec'*sin(Dip); 
Z_points_repmat = repmat(Z_points,1,length(L_vec));

points_loc = [ ...
    repmat(X_poits,length(W_vec),1), ...
    reshape(Y_points_repmat',[],1), ...
    reshape(Z_points_repmat',[],1)];

%%% Calculation of Various "R"s and "Distance Parameters"  for GMPMs
Rrup_points = sum((repmat(Site1_loc,length(W_vec)*length(L_vec),1) - points_loc).^2,2).^0.5;

points_loc_Projec_JB = points_loc; 
points_loc_Projec_JB(:,3) = 0;

Rjb_points = sum((repmat(Site1_loc,length(W_vec)*length(L_vec),1) - points_loc_Projec_JB).^2,2).^0.5;

Ztor_points = points_loc(:,3);
 
%%% Calculation of Azimuth Alpha and Rx based on Kaklamanos et. al. 2011
for i = 1 : length(L_vec) * length(W_vec)

    % Step 1: Determination of the "9 Cases" for each point
    D_dL_point_X =  [points_loc_Projec_JB(i,1)-dL/2, points_loc_Projec_JB(i,1)+dL/2];
    D_dL_point_Y =  [points_loc_Projec_JB(i,2)-dL/2, points_loc_Projec_JB(i,2)+dL/2];

    if Site1_loc(2) < D_dL_point_Y(1)*0.9999 % Cases 1, 4, and 7

        if Site1_loc(1) <= D_dL_point_X(2)*1.0001 && Site1_loc(1) >= D_dL_point_X(1)*0.9999
            Case(i) = 4;
        elseif Site1_loc(1) > D_dL_point_X(2)*1.0001
            Case(i) = 1;
        elseif Site1_loc(1) < D_dL_point_X(1)*0.9999
            Case(i) = 7;    
        else
            error('Something is wrong Babak!!!')
        end

    elseif  Site1_loc(2) > D_dL_point_Y(2)*1.0001 % Cases 3, 6, and 9

        if Site1_loc(1) <= D_dL_point_X(2)*1.0001 && Site1_loc(1) >= D_dL_point_X(1)*0.9999
            Case(i) = 6;
        elseif Site1_loc(1) > D_dL_point_X(2)*1.0001
            Case(i) = 3;
        elseif Site1_loc(1) < D_dL_point_X(1)*0.9999
            Case(i) = 9;    
        else
            error('Something is wrong Babak!!!')
        end

    elseif  Site1_loc(2) <= D_dL_point_Y(2)*1.0001 && Site1_loc(2) >= D_dL_point_Y(1)*0.9999 % Cases 2, 5, and 8

        if Site1_loc(1) <= D_dL_point_X(2)*1.0001 && Site1_loc(1) >= D_dL_point_X(1)*0.9999
            Case(i) = 5;
        elseif Site1_loc(1) > D_dL_point_X(2)*1.0001
            Case(i) = 2;
        elseif Site1_loc(1) < D_dL_point_X(1)*0.9999
            Case(i) = 8;    
        else
            error('Something is wrong Babak!!!')
        end    
    end

    % Step 2: Determination of alpha and Rx for each case
    angle = atan2(points_loc_Projec_JB(i,2) - Site1_loc(2), points_loc_Projec_JB(i,1) - Site1_loc(1));
    angle = abs(angle);  
    absolute_alpha = min(angle, pi - angle);  

    if Case(i) == 1
        alpha = -absolute_alpha;
        theta = alpha;        
    elseif Case(i) == 2 || Case(i) == 3
        alpha = absolute_alpha;
        theta = alpha;
    elseif Case(i) == 4
        alpha = -pi/2;
        theta = alpha;
    elseif Case(i) == 5 || Case(i) == 6
        alpha = pi/2; 
        theta = alpha;
    elseif Case(i) == 7 
        alpha = -pi+absolute_alpha;
        theta = alpha;
    elseif Case(i) == 8 || Case(i) == 9
        alpha = pi-absolute_alpha; 
        theta = absolute_alpha;
    else
        error('Some thing is wrong Babak!')
    end  

    Rx_points(i) = Rjb_points(i)*sin(theta);
    Cases_points(i) = Case(i);
    alpha_points(i) = alpha;
    
end

Ry_points = points_loc(:,1) - Site1_loc(1);

%%% Calcualtion of the Grid points location Reletive to Fault
W_vec_repmat = repmat(W_vec',1,length(L_vec));

points_fault_loc = [ ...
    repmat(L_vec',length(W_vec),1), ...
    reshape(W_vec_repmat',[],1)];

%% Store outputs
FG = struct();

FG.Fault_Type = Fault_Type;

FG.W = W;
FG.L = L;
FG.Hf = Hf;
FG.Dip = Dip;
FG.dL = dL;

FG.Xs = Xs;
FG.Ys = Ys;
FG.Site1_loc = Site1_loc;

FG.W_vec = W_vec;
FG.L_vec = L_vec;

FG.points_loc = points_loc;
FG.points_loc_Projec_JB = points_loc_Projec_JB;
FG.points_fault_loc = points_fault_loc;

FG.Rrup_points = Rrup_points;
FG.Rjb_points = Rjb_points;
FG.Ztor_points = Ztor_points;
FG.Rx_points = Rx_points;
FG.Cases_points = Cases_points;
FG.alpha_points = alpha_points;
FG.Ry_points = Ry_points;

end