function [lMT, MA] = CalculateMuscleTendonLengthAndMomentArms(x, data_exp, coeff_LMT_ma)
% Calculate MuscleTendonLength and Moment arms 

% Input
offset   = data_exp.offset;  
m_offset = mean(offset);

lMT_ext = coeff_LMT_ma(1,1) + coeff_LMT_ma(2,1)*(x+m_offset) + coeff_LMT_ma(3,1)*(x+m_offset).^2 + coeff_LMT_ma(4,1)*(x+m_offset).^3;
MA_ext  = -coeff_LMT_ma(2,1) + -coeff_LMT_ma(3,1)*(x+m_offset) + -coeff_LMT_ma(4,1)*(x+m_offset).^2;
lMT_flex= coeff_LMT_ma(1,2) + coeff_LMT_ma(2,2)*(x+m_offset) + coeff_LMT_ma(3,2)*(x+m_offset).^2 + coeff_LMT_ma(4,2)*(x+m_offset).^3;
MA_flex = -coeff_LMT_ma(2,2) + -coeff_LMT_ma(3,2)*(x+m_offset) + -coeff_LMT_ma(4,2)*(x+m_offset).^2;

lMT = [lMT_ext; lMT_flex];
MA  = [MA_ext; MA_flex];
end

