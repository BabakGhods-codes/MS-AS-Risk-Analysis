function [rup] = NGA_GMPM_parameter_estimation_accurate(M, Rrup, Rjb, Rx, Ry, Ztor, Zhyp, Dip, W, Lambda, Vs30, Z1o0, Z2o5, Region, Fault_Type, FVS30, fas, Fhw)
    rup.M = M;
    
    if strcmp(Fault_Type,'Unspecified')
        rup.Fault_Type = 0;
    elseif strcmp(Fault_Type,'Strike Slip')
        rup.Fault_Type = 1;
    elseif strcmp(Fault_Type,'Normal')
        rup.Fault_Type = 2;
    elseif strcmp(Fault_Type,'Reverse')
        rup.Fault_Type = 3;
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


    rup.Dip_degree     = Dip/pi*180;
    rup.lambda_degree  = Lambda/pi*180;
    rup.Fhw            = Fhw;

    rup.W       = W;
    rup.Rrup    = Rrup;
    rup.Rjb     = Rjb;
    rup.Rx      = Rx;
    rup.Ry      = Ry; 
    rup.Ztor    = Ztor;
    rup.Zhyp    = Zhyp;
    rup.Zbot    = []; % We do not need this because we have the rupture width, "W"

    rup.Vs30    = Vs30;
    rup.Z1o0    = Z1o0; 
    rup.Z25     = Z2o5; 

    rup.FVS30   = FVS30;  % FVS30 = 0 for Vs30 is inferred from geology; FVS30 = 1 for measured  Vs30
    rup.fas     = fas;  % Flag for aftershocks
end