function et_score=timeE_Tsallis_entropy(time_pow, q)

if sum(time_pow)>1
    time_pow=time_pow/sum(time_pow);
end

if q~=1
 et_score=0;
    for j=1:length(time_pow)
        if time_pow(j)>0
            et_score=et_score+(time_pow(j))^q;
        end
    end
    et_score=(1/(q-1))*(1-et_score);
    t_max=(length(time_pow(time_pow>0))^(1-q)-1)/(1-q);
    et_score=et_score/t_max;

else
        et_score=0;
        for j=1:length(time_pow)
            if time_pow(j)>0
                et_score=et_score-(time_pow(j))*log(time_pow(j));
            end
        end
        t_max=(length(time_pow(time_pow>0)));
        et_score=et_score/log(t_max);
end


end