% https://en.wikipedia.org/wiki/Framingham_Risk_Score

MDD = readtable('TNM_Subject_IDs.csv');
C = readtable('Predictors.tsv','FileType','text');
C = standardizeMissing(C,-3); C = standardizeMissing(C,-1); C = standardizeMissing(C,-818); C = standardizeMissing(C,-121); C = standardizeMissing(C,-7);
[~,ID,IC] = intersect(MDD.eid,C.eid);
D = readtable('Subgroups.tsv','FileType','text');
D = standardizeMissing(D,-3); D = standardizeMissing(D,-1); D = standardizeMissing(D,-818); D = standardizeMissing(D,-121); D = standardizeMissing(D,-7);
sex = D.x31_0_0;
age = D.x21003_2_0;
sex = sex(ID);
age = age(ID);
clear D
cholesterol = C.x30690_0_0(IC);
smoking = C.x20116_2_0(IC);
hdl = C.x30760_0_0(IC);
systolic = C.x4080_2_0(IC);
bpmedfem = readtable('bpmedfem.csv');
bpmeds_female = bpmedfem.x6153_2_0; % female version (okay this is here
% cause you need to know if they're on medication or not
bpmeds_male = C.x6177_2_0(IC); % male version
% WMH = WMH(MDD); WMH = WMH(ID);

Framingham = zeros(length(IC),1);
for n = 1:length(IC)
    if sex(n) == 0 % Female
        % Points for age
        if age>=20 & age<=34
            Framingham(n) = Framingham(n)-7;
        elseif age>=35 & age<=39
            Framingham(n) = Framingham(n)-3;
        elseif age>=40 & age<=44
            Framingham(n) = Framingham(n);
        elseif age>=45 & age<=49
            Framingham(n) = Framingham(n)+3;
        elseif age>=50 & age<=54
            Framingham(n) = Framingham(n)+6;
        elseif age>=55 & age<=59
            Framingham(n) = Framingham(n)+8;
        elseif age>=60 & age<=64
            Framingham(n) = Framingham(n)+10;
        elseif age>=65 & age<=69
            Framingham(n) = Framingham(n)+12;
        elseif age>=70 & age<=74
            Framingham(n) = Framingham(n)+14;
        elseif age>=75
            Framingham(n) = Framingham(n)+16;
        elseif isnan(age(n))
            Framingham(n) = NaN;
        end
        % Points for total cholesterol(n)
        if age>=20 & age<=39
            if cholesterol(n)>=160 & cholesterol(n)<200
                Framingham(n) = Framingham(n)+4;
            elseif cholesterol(n)>=200 & cholesterol(n)<240
                Framingham(n) = Framingham(n)+8;
            elseif cholesterol(n)>=240 & cholesterol(n)<280
                Framingham(n) = Framingham(n)+11;
            elseif cholesterol(n)>280
                Framingham(n) = Framingham(n)+13;
            elseif isnan(cholesterol(n))
            Framingham(n) = NaN;
            end
            if smoking(n) == 2
                Framingham(n) = Framingham(n)+9;
            elseif isnan(smoking(n))
                Framingham(n) = NaN;
            end
        elseif age>=40 & age<=49
            if cholesterol(n)>=160 & cholesterol(n)<200
                Framingham(n) = Framingham(n)+3;
            elseif cholesterol(n)>=200 & cholesterol(n)<240
                Framingham(n) = Framingham(n)+6;
            elseif cholesterol(n)>=240 & cholesterol(n)<280
                Framingham(n) = Framingham(n)+8;
            elseif cholesterol(n)>280
                Framingham(n) = Framingham(n)+10;
            elseif isnan(cholesterol(n))
                Framingham(n) = NaN;
            end
            if smoking(n) == 2
                Framingham(n) = Framingham(n)+7;
            end
        elseif age>=50 & age<=59
            if cholesterol(n)>=160 & cholesterol(n)<200
                Framingham(n) = Framingham(n)+2;
            elseif cholesterol(n)>=200 & cholesterol(n)<240
                Framingham(n) = Framingham(n)+4;
            elseif cholesterol(n)>=240 & cholesterol(n)<280
                Framingham(n) = Framingham(n)+5;
            elseif cholesterol(n)>280
                Framingham(n) = Framingham(n)+7;
            elseif isnan(cholesterol(n))
                Framingham(n) = NaN;
            end
            if smoking(n) == 2
                Framingham(n) = Framingham(n)+4;
            end
        elseif age>=60 & age<=69
            if cholesterol(n)>=160 & cholesterol(n)<200
                Framingham(n) = Framingham(n)+1;
            elseif cholesterol(n)>=200 & cholesterol(n)<240
                Framingham(n) = Framingham(n)+2;
            elseif cholesterol(n)>=240 & cholesterol(n)<280
                Framingham(n) = Framingham(n)+3;
            elseif cholesterol(n)>280
                Framingham(n) = Framingham(n)+4;
            elseif isnan(cholesterol(n))
                Framingham(n) = NaN;
            end
            if smoking(n) == 2
                Framingham(n) = Framingham(n)+2;
            elseif isnan(smoking(n))
                Framingham(n) = NaN;
            end
        elseif age>=70 
            if cholesterol(n)>=160 & cholesterol(n)<200
                Framingham(n) = Framingham(n)+1;
            elseif cholesterol(n)>=200 & cholesterol(n)<240
                Framingham(n) = Framingham(n)+1;
            elseif cholesterol(n)>=240 & cholesterol(n)<280
                Framingham(n) = Framingham(n)+2;
            elseif cholesterol(n)>280
                Framingham(n) = Framingham(n)+2;
            elseif isnan(cholesterol(n))
                Framingham(n) = NaN;
            end
            if smoking(n) == 2
                Framingham(n) = Framingham(n)+1;
            elseif isnan(smoking(n))
                Framingham(n) = NaN;
            end
        end
        if hdl(n)>60
            Framingham(n) = Framingham(n)-1;
        elseif hdl(n)>=40 & hdl(n)<50
            Framingham(n) = Framingham(n)+1;
        elseif hdl(n)<40
            Framingham(n) = Framingham(n)+2;
        elseif isnan(hdl(n))
            Framingham(n) = NaN;
        end
        if bpmeds_female == 2
            if systolic(n)>=120 & systolic(n)<130
               Framingham(n) = Framingham(n)+3;
            elseif systolic(n)>=130 & systolic(n)<140
               Framingham(n) = Framingham(n)+4;
            elseif systolic(n)>=140 & systolic(n)<160
               Framingham(n) = Framingham(n)+5;
            elseif systolic(n)>=120
                Framingham(n) = Framingham(n)+6;
            elseif isnan(systolic(n))
                Framingham(n) = NaN;
            end
        else
            if systolic(n)>=120 & systolic(n)<130
               Framingham(n) = Framingham(n)+1;
            elseif systolic(n)>=130 & systolic(n)<140
               Framingham(n) = Framingham(n)+2;
            elseif systolic(n)>=140 & systolic(n)<160
               Framingham(n) = Framingham(n)+3;
            elseif systolic(n)>=120
                Framingham(n) = Framingham(n)+4;
            elseif isnan(systolic(n))
                Framingham(n) = NaN;
            end
        end      
    
    elseif sex(n) == 1 % Male
        if age>=20 & age<=34
            Framingham(n) = Framingham(n)-9;
        elseif age>=35 & age<=39
            Framingham(n) = Framingham(n)-4;
        elseif age>=40 & age<=44
            Framingham(n) = Framingham(n);
        elseif age>=45 & age<=49
            Framingham(n) = Framingham(n)+3;
        elseif age>=50 & age<=54
            Framingham(n) = Framingham(n)+6;
        elseif age>=55 & age<=59
            Framingham(n) = Framingham(n)+8;
        elseif age>=60 & age<=64
            Framingham(n) = Framingham(n)+10;
        elseif age>=65 & age<=69
            Framingham(n) = Framingham(n)+11;
        elseif age>=70 & age<=74
            Framingham(n) = Framingham(n)+12;
        elseif age>=75
            Framingham(n) = Framingham(n)+13;
        elseif age(n) == NaN
            Framingham(n) = NaN;
        end
        % Points for total cholesterol(n)
        if age>=20 & age<=39
            if cholesterol(n)>=160 & cholesterol(n)<200
                Framingham(n) = Framingham(n)+4;
            elseif cholesterol(n)>=200 & cholesterol(n)<240
                Framingham(n) = Framingham(n)+7;
            elseif cholesterol(n)>=240 & cholesterol(n)<280
                Framingham(n) = Framingham(n)+9;
            elseif cholesterol(n)>280
                Framingham(n) = Framingham(n)+11;
            elseif isnan(cholesterol(n))
                Framingham(n) = NaN;
            end
            if smoking(n) == 2
                Framingham(n) = Framingham(n)+8;
            elseif isnan(smoking(n))
                Framingham(n) = NaN;
            end
        elseif age>=40 & age<=49
            if cholesterol(n)>=160 & cholesterol(n)<200
                Framingham(n) = Framingham(n)+3;
            elseif cholesterol(n)>=200 & cholesterol(n)<240
                Framingham(n) = Framingham(n)+5;
            elseif cholesterol(n)>=240 & cholesterol(n)<280
                Framingham(n) = Framingham(n)+6;
            elseif cholesterol(n)>280
                Framingham(n) = Framingham(n)+8;
            elseif isnan(cholesterol(n))
                Framingham(n) = NaN;
            end
            if smoking(n) == 2
                Framingham(n) = Framingham(n)+5;
            elseif isnan(smoking(n))
                Framingham(n) = NaN;
            end
        elseif age>=50 & age<=59
            if cholesterol(n)>=160 & cholesterol(n)<200
                Framingham(n) = Framingham(n)+2;
            elseif cholesterol(n)>=200 & cholesterol(n)<240
                Framingham(n) = Framingham(n)+3;
            elseif cholesterol(n)>=240 & cholesterol(n)<280
                Framingham(n) = Framingham(n)+4;
            elseif cholesterol(n)>280
                Framingham(n) = Framingham(n)+5;
            elseif isnan(cholesterol(n))
                Framingham(n) = NaN;
            end
            if smoking(n) == 2
                Framingham(n) = Framingham(n)+3;
            elseif isnan(smoking(n))
                Framingham(n) = NaN;
            end
        elseif age>=60 & age<=69
            if cholesterol(n)>=160 & cholesterol(n)<200
                Framingham(n) = Framingham(n)+1;
            elseif cholesterol(n)>=200 & cholesterol(n)<240
                Framingham(n) = Framingham(n)+1;
            elseif cholesterol(n)>=240 & cholesterol(n)<280
                Framingham(n) = Framingham(n)+2;
            elseif cholesterol(n)>280
                Framingham(n) = Framingham(n)+3;
            elseif isnan(cholesterol(n))
                Framingham(n) = NaN;
            end
            if smoking(n) == 2
                Framingham(n) = Framingham(n)+1;
            elseif isnan(smoking(n))
                Framingham(n) = NaN;
            end
        elseif age>=70 
            if cholesterol(n)>=160 & cholesterol(n)<200
                Framingham(n) = Framingham(n)+0;
            elseif cholesterol(n)>=200 & cholesterol(n)<240
                Framingham(n) = Framingham(n)+0;
            elseif cholesterol(n)>=240 & cholesterol(n)<280
                Framingham(n) = Framingham(n)+1;
            elseif cholesterol(n)>280
                Framingham(n) = Framingham(n)+1;
            elseif isnan(cholesterol(n))
                Framingham(n) = NaN;
            end
            if smoking(n) == 2
                Framingham(n) = Framingham(n)+1;
            elseif isnan(smoking(n))
                Framingham(n) = NaN;
            end
        end
        if hdl(n)>60
            Framingham(n) = Framingham(n)-1;
        elseif hdl(n)>=40 & hdl(n)<50
            Framingham(n) = Framingham(n)+1;
        elseif hdl(n)<40
            Framingham(n) = Framingham(n)+2;
         elseif isnan(hdl(n))
            Framingham(n) = NaN;
        end
         if bpmeds_male == 2
            if systolic(n)>=120 & systolic(n)<130
               Framingham(n) = Framingham(n)+1;
            elseif systolic(n)>=130 & systolic(n)<140
               Framingham(n) = Framingham(n)+2;
            elseif systolic(n)>=140 & systolic(n)<160
               Framingham(n) = Framingham(n)+2;
            elseif systolic(n)>=120
                Framingham(n) = Framingham(n)+3;
            elseif isnan(systolic(n))
                Framingham(n) = NaN;
            end
        else
            if systolic(n)>=120 & systolic(n)<130
               Framingham(n) = Framingham(n)+0;
            elseif systolic(n)>=130 & systolic(n)<140
               Framingham(n) = Framingham(n)+1;
            elseif systolic(n)>=140 & systolic(n)<160
               Framingham(n) = Framingham(n)+1;
            elseif systolic(n)>=120
                Framingham(n) = Framingham(n)+2;
            elseif isnan(systolic(n))
                Framingham(n) = NaN;
            end
        end
    end
end

%%
mtrx = [6153];
writematrix(mtrx, 'bpmedfem.csv')

%%
%fix the index
[~, ia] = intersect(MDD.eid, anhedonia_eid(:,2)); Fram{1} = Framingham(ia, :);
[~, im] = intersect(MDD.eid, mood_eid(:,2)); Fram{2} = Framingham(im, :);
[~, is] = intersect(MDD.eid, somatic_eid(:,2)); Fram{3} = Framingham(is, :);
[~, ic] = intersect(MDD.eid, chronic_eid(:,2)); Fram{4} = Framingham(ic, :);
[~, il] = intersect(MDD.eid, lateonset_eid(:,2)); Fram{5} = Framingham(il, :);
[~, ise] = intersect(MDD.eid, severe_eid(:,2)); Fram{6} = Framingham(ise, :);
%remove missing data
for n = 1:6
    nans{n} = isnan(Fram{n});
end
for n = 1:6
    fram_nanless{n} = Fram{n}(~any(nans{n},2),:);
    kms_fram{n} = kms{n}(~any(nans{n},2),:);
end
%run anova
for i=1:6
        p_fram_6(i) = anova1(fram_nanless{i}, kms_fram{i},  'display', 'off');
end
%none of them were significantly different
