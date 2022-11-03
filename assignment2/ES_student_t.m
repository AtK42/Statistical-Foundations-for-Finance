function ES = ES_student_t(alpha, df, location, scale)
    c01=tinv(alpha , df); 
    % location zero, scale 1 ES for student t, calculating it using first the analytic exact expression: 
    ES_01_analytic = -tpdf(c01,df)/tcdf(c01,df) * (df+c01^2)/(df-1);
    ES = location + scale * ES_01_analytic;
end
