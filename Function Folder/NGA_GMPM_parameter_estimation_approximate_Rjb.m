function [rup ] = NGA_GMPM_parameter_estimation_approximate_Rjb(M, Rjb, Dip, Vs30, Region, Fault_Type, FVS30, fas, Fhw)
% The approximation is done based on the Kaklamanos et. al. 2011: Estimating Unknown Input Parameters
% when Implementing the NGA Ground-Motion Prediction Equations in Engineering Practice
% It is assumed that we only know the M, Rjb, Fhw, Fault_Type, and Vs30. Also
% FVS30 is assumed to be 1 (Vs30 in infered) and fas is considered as 0.
% For Dip, it can be assinged directly, otherwise it is estimated
% The region is known too
    rup.M = M;
    
    if strcmp(Fault_Type,'Unspecified')
        rup.Fault_Type = 0;
        Lambda = pi/4;
        a_RW = -1.01; b_RW = 0.32; sigm_RW = 0.15;
        Zhyp = 7.08 + 0.61*M;
    elseif strcmp(Fault_Type,'Strike-Slip')
        rup.Fault_Type = 1;
        Lambda = 0;
        a_RW = -0.76; b_RW = 0.27; sigm_RW = 0.14;
        Zhyp = 5.63 + 0.68*M;
    elseif strcmp(Fault_Type,'Normal')
        rup.Fault_Type = 2;
        Lambda = -pi/2;
        a_RW = -1.14; b_RW = 0.35; sigm_RW = 0.12;
        Zhyp = 11.24 - 0.2*M;
    elseif strcmp(Fault_Type,'Reverse')
        rup.Fault_Type = 3;
        Lambda = pi/2;
        a_RW = -1.61; b_RW = 0.41; sigm_RW = 0.15;
        Zhyp = 11.24 - 0.2*M;
    else
       error('The fault type is not defined or is misspelled') 
    end

    if strcmp(Region,'global')
        rup.Region = 0;
    elseif strcmp(Region,'California')
        rup.Region = 1;
    elseif strcmp(Region,'Japan')
        rup.Region = 2;
    elseif strcmp(Region,'China or Turkey')
        rup.Region = 3;
    elseif strcmp(Region,'Italy')
        rup.Region = 4;    
    else
       error('The Region type is not defined or is misspelled') 
    end
    %%% Estimation of dip %%%
    if isempty(Dip)
        if strcmp(Fault_Type,'Unspecified')
            Dip = pi * 60/180;
        elseif strcmp(Fault_Type,'Strike-Slip')
            Dip = pi/2;
        elseif strcmp(Fault_Type,'Normal')
            Dip = pi * 50/180;
        elseif strcmp(Fault_Type,'Reverse')
            Dip = pi * 40/180;
        else
           error('The fault type is not defined or is misspelled') 
        end 
    end
    %%% Estimation of alpha %%%
    if Fhw == 1
        rup.alpha = pi * 50/180; % Based on Kaklamanos et. al. 2011 recommendation 
    else
        rup.alpha = pi * -50/180; % Based on Kaklamanos et. al. 2011 recommendation 
    end
    
    
    rup.Dip_degree     = Dip/pi*180;
    rup.lambda_degree  = Lambda/pi*180;
    rup.Fhw            = Fhw;

    rup.W       = 10 ^(a_RW + b_RW * M); % Rupture width is calculated based on the Wells and Coppersmith 1994
    rup.Zbot    = []; % We do not need this because we have the rupture width, "W"
    rup.Zhyp    = Zhyp;
    rup.Ztor    = max([Zhyp-0.6*rup.W*sin(Dip),0]);
    rup.FVS30   = FVS30;  % FVS30 = 0 for Vs30 is inferred from geology; FVS30 = 1 for measured  Vs30
    rup.fas     = fas;  % Flag for aftershocks
    
    rup.Vs30    = Vs30;
    if Vs30 < 180
        rup.Z1o0_AS = exp(6.745);
    elseif Vs30 >= 180 && Vs30 <= 500
        rup.Z1o0_AS = exp(6.745 - 1.35 * log(Vs30/180));
    elseif Vs30 > 500
        rup.Z1o0_AS = exp(5.394 - 4.48 * log(Vs30/500));
    end
    rup.Z1o0_CY = exp(28.5 - 3.82/8 * log(Vs30^8 + 378.7^8));
    
    rup.Z25 = 519 + 3.595 * rup.Z1o0_AS;
    
    
    rup.Z1o0_AS = rup.Z1o0_AS /1000; % (km)
    rup.Z1o0_CY = rup.Z1o0_CY /1000; % (km)
    rup.Z25     = rup.Z25 /1000; % (km)
    
    % Calculation of Tjb, Rx, and Ry with back-calculation of the equations    
    rup.Rjb     = Rjb;
    [Rrup, Rx_estimate, Ry_estimate] = distance_Distance_Kaklamanos2011(Rjb, Dip, Fault_Type, rup);
    rup.Rx      = Rx_estimate;
    rup.Ry      = Ry_estimate;     
    rup.Rrup    = Rrup;
    % Just for testing the accuracy 
    % Rjb_Test = back_calculate_Rjb(Rrup, Dip, Fault_Type, rup);
    % error_Rjb = abs(Rjb_Test - Rjb)/Rjb;
    % error_Rjb
    % if error_Rjb > 0.001
    %     error('Babak: Input and back-calculated Rjb do not match')
    % end
    
function [Rrup] = distance_Rjb_Kaklamanos2011(Rjb, Dip, Fault_Type, rup)
    % Calculation of Rx (only equation 7, 8 and 12 are considered becasue the alpha is only 50 or -50)
    if ~strcmp(Fault_Type,'Strike-Slip')
        if rup.alpha > 0
            if rup.W*cos(Dip) >= Rjb * abs(tan(rup.alpha)) % Case 2
                 Rx = Rjb * abs(tan(rup.alpha)); % Eq. 7
            else % Case 3
                 Rx = Rjb*tan(rup.alpha)*cos(rup.alpha-asin(rup.W*cos(Dip)*cos(rup.alpha)/Rjb)); % Eq. 8
            end
        elseif rup.alpha < 0 % Cases 1, 4, 7
            Rx =  Rjb * sin(rup.alpha);
        else
            error('It is not specified !!!')
        end
    else
        Rx =  Rjb * sin(rup.alpha); % Eq. 13
    end 
    % Calculation of Rrup_p (Eqs 14 to 17)
    if Rx < rup.Ztor*tan(Dip)
        Rrup_p = (Rx^2 + rup.Ztor^2)^0.5; % Eq. 15
    elseif Rx >= rup.Ztor*tan(Dip) && Rx <= rup.Ztor*tan(Dip) + rup.W*sec(Dip)
        Rrup_p = Rx*sin(Dip) + rup.Ztor*cos(Dip); % Eq. 16
    elseif Rx > rup.Ztor*tan(Dip) + rup.W*sec(Dip)
        Rrup_p = ((Rx-rup.W*cos(Dip))^2  + (rup.Ztor+rup.W*sin(Dip))^2)^0.5; % Eq. 17
    else
        error('Something is wrong Babak!!!')
    end
    % Calculation of Ry (only Eq. 20 is considered because other alpha values are not considered)
    Ry = abs(Rx * cot(rup.alpha));
    % Rrup
    Rrup = (Rrup_p^2 + Ry^2)^0.5; % Eq. 14
    
end

function Rjb = back_calculate_Rjb(Rrup_target, Dip, Fault_Type, rup)
    % Objective function: minimize the difference between calculated and target Rrup
    fun = @(Rjb) abs(distance_Rjb_Kaklamanos2011(Rjb, Dip, Fault_Type, rup) - Rrup_target);

    % Initial guess for Rjb
    Rjb0 = Rrup_target*0.8; % Replace with appropriate initial value

    % Optimization options
    options = optimoptions('fmincon','Display','off');

    % Lower and upper bounds for Rjb (if applicable)
    lb = [0]; % Lower bound
    ub = []; % Upper bound

    % Perform optimization
    Rjb = fmincon(fun, Rjb0, [], [], [], [], lb, ub, [], options);
end

function [Rrup, Rx, Ry] = distance_Distance_Kaklamanos2011(Rjb, Dip, Fault_Type, rup)
    % Calculation of Rx (only equation 7, 8 and 12 are considered becasue the alpha is only 50 or -50)
    if ~strcmp(Fault_Type,'Strike-Slip')
        if rup.alpha > 0
            if rup.W*cos(Dip) >= Rjb * abs(tan(rup.alpha)) % Case 2
                 Rx = Rjb * abs(tan(rup.alpha));
            else % Case 3
                 Rx = Rjb*tan(rup.alpha)*cos(rup.alpha-asin(rup.W*cos(Dip)*cos(rup.alpha)/Rjb));
            end
        elseif rup.alpha < 0 % Cases 1, 4, 7
            Rx =  Rjb * sin(rup.alpha);
        else
            error('It is not specified !!!')
        end
    else
        Rx =  Rjb * sin(rup.alpha);
    end 
    % Calculation of Rrup_p (Eqs 14 to 17)
    if Rx < rup.Ztor*tan(Dip)
        Rrup_p = (Rx^2 + rup.Ztor^2)^0.5;
    elseif Rx >= rup.Ztor*tan(Dip) && Rx <= rup.Ztor*tan(Dip) + rup.W*sec(Dip)
        Rrup_p = Rx*sin(Dip) + rup.Ztor*cos(Dip);
    elseif Rx > rup.Ztor*tan(Dip) + rup.W*sec(Dip)
        Rrup_p = ((Rx-rup.W*cos(Dip))^2 + (rup.Ztor+rup.W*sin(Dip))^2)^0.5;
    else
        error('Something is wrong Babak!!!')
    end
    % Calculation of Ry (only Eq. 20 is considered because other alpha values are not considered)
    Ry = abs(Rx * cot(rup.alpha));
    % Rrup
    Rrup = (Rrup_p^2 + Ry^2)^0.5;
    
end
%     if ~strcmp(Fault_Type,'Strike-Slip')
%         if alpha > 0
%             Rx_over_Rjb = abs(tan(alpha));
%         elseif alpha < 0
%             Rx_over_Rjb = sin(alpha);
%         else
%             error('It is not specified !!!')
%         end
%     else
%         Rx_over_Rjb = sin(alpha);
%     end 
%     
%     Ry_over_Rjb = cot(alpha) * Rx_over_Rjb;
%     
%     % Estimating some Rrup to see which one is correct
%     Rrup_p_1 = Rrup / (Rx_over_Rjb^2 + rup.Ztor^2 + Ry_over_Rjb^2);% Eq. 15 Kaklamanos et. al. 2011
%     Rrup_p_2 = 
%     
    
end