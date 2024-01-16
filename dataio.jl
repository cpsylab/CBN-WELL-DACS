using StatFiles, XLSX, CSV, DataFrames, Dates
using Pipe, Glob
using Statistics

""" 
    sas_to_csv(fdir::String="data") 

    Function to convert SAS files to CSVs 
    Will enable us to look at files in spreadsheet program visually 
"""
function sas_to_csv(fdir::String="data") 
    for f ∈ glob(fdir*"/*.sas7bdat")
        df = @pipe f |> load |> DataFrame
        fout = fdir * "/" * basename(f) * ".csv"
        CSV.write(fout, df)
    end
end

""" 
    extract_questionnaires_to_csv()

    Load the questionnaires and extract them to csv files for analysis.
"""
function extract_questionnaires_to_csv()
    qs_df = @pipe DataFrame(load("data/qs.sas7bdat")) |> 
        select(_, Not([
            :STUDYID, :DOMAIN, :QSSEQ, :QSSPID, :QSSCAT, 
            :QSSTRESU, :QSSTAT, :QSREASND, :QSNAM, :QSBLFL,	
            :QSEVAL, :VISITNUM, :EPOCH, :VISITDY, :QSDY, 
            :QSEVLINT, :QSEVINTX, :QSSTRESN, :QSORRESU,
            :QSORRES, :QSTEST, :QSGRPID
        ])) 
    qs_df.VISIT = map(x -> replace(x, r"\.[0-9].*$" => "") , qs_df.VISIT)
    qs_df = qs_df[.!nonunique(qs_df, collect(1:size(qs_df, 2))), :]
    qs_df.QSDTC = map(x -> replace(x, r"T.*$" => "") , qs_df.QSDTC)

    # QIDS
    extractscale(scalename; combine=only) = @pipe qs_df |> 
        _[.!nonunique(_, collect(1:size(_, 2))), :] |> 
        filter(:QSCAT => x -> x == scalename, _) |> 
        unstack(_, [:USUBJID, :QSDTC, :VISIT], :QSTESTCD, :QSSTRESC; combine=combine)

    bdi = extractscale("BDI")
    brian = extractscale("BRIAN")
    bsi53 = extractscale("BSI-53")
    cgi = extractscale("CGI")
    chrt = extractscale("CHRT")
    eq5d5l = extractscale("EQ-5D-5L")
    etisrsf = extractscale("ETISR-SF")
    gad7 = extractscale("GAD-7")
    ipaq = extractscale("IPAQ-SF PHONE VERSION")
    leaps = extractscale("LEAPS")
    madrs = extractscale("MADRS")
    mos = extractscale("MOS SLEEP REVISED")
    pfibs = extractscale("P-FIBS")
    paq = extractscale("PAQ")
    psqi = extractscale("PSQI")
    pss = extractscale("PSS")
    qlesq = extractscale("Q-LES-Q-SF")
    qids = extractscale("QIDS-SR")
    rlcst = extractscale("RLCST")
    rsat = extractscale("RSAT")
    shaps = extractscale("SHAPS")
    whodas = extractscale("WHODAS 12-ITEM SELF")
    whoqol = extractscale("WHOQOL-BREF")
    fsr = extractscale("FSR-SA")
    mini = extractscale("MINI 700"; combine=first)
    sds = extractscale("SDS"; combine=first)

    # Write scales to CSVs 
    CSV.write("data/bdi.csv", bdi)
    CSV.write("data/brian.csv", brian)
    CSV.write("data/bsi53.csv", bsi53)
    CSV.write("data/cgi.csv", cgi)
    CSV.write("data/chrt.csv", chrt)
    CSV.write("data/eq5d5l.csv", eq5d5l)
    CSV.write("data/etisrsf.csv", etisrsf)
    CSV.write("data/gad7.csv", gad7)
    CSV.write("data/ipaq.csv", ipaq)
    CSV.write("data/leaps.csv", leaps)
    CSV.write("data/madrs.csv", madrs)
    CSV.write("data/mos.csv", mos)
    CSV.write("data/pfibs.csv", pfibs)
    CSV.write("data/paq.csv", paq)
    CSV.write("data/psqi.csv", psqi)
    CSV.write("data/pss.csv", pss)
    CSV.write("data/qlesq.csv", qlesq)
    CSV.write("data/qids.csv", qids)
    CSV.write("data/rlcst.csv", rlcst)
    CSV.write("data/rsat.csv", rsat)
    CSV.write("data/shaps.csv", shaps)
    CSV.write("data/whodas.csv", whodas)
    CSV.write("data/whoqol.csv", whoqol)
    CSV.write("data/fsr.csv", fsr)
    CSV.write("data/mini.csv", mini)
    CSV.write("data/sds.csv", sds)

    # Reload qids and convert missings to zeros then re-save 
    qids = CSV.read("data/qids.csv", missingstring=["NA", ""], DataFrame)
    for i ∈ 4:size(qids, 2)
        qids[ismissing.(qids[:, i]), i] .= 0
    end
    CSV.write("data/qids.csv", qids)
end

"""
    load_subject_ids()

    Load the subject IDs from the data directory and find those who are excluded. 
    Will also identify the CBN1 participants. 
    Saves a file with excluded subjects, and a file with included subjects.
    Returns a dataframe with the included subjects.
"""
function load_subject_ids()
    # Load and convert the NA values to empty, then save as csv. 
    subject_list = DataFrame(XLSX.readtable("data/Wellness.DEID.Subject.List.xlsx", "Wellness.DEID.Subject.List"))[:, 
    [:subject_id, :SUBJLABEL, :Relapse_Status, :Inclusion, :Status]
    ]
    CSV.write("data/subject_list.csv", subject_list)
    
    # Find the excluded subjects
    subject_list = CSV.read("data/subject_list.csv", missingstring="NA", DataFrame)
    subject_list.Inclusion = map(x -> replace(x, r" - .*" => ""), subject_list.Inclusion)
    included_subjects = subject_list[subject_list.Inclusion .== "Include", :]
    excluded_subjects = subject_list[subject_list.Inclusion .== "Exclude", :]

    # Save the excluded subjects
    CSV.write("results/excluded_subjects.csv", excluded_subjects)
    CSV.write("results/included_subjects.csv", included_subjects)

    return included_subjects
end

"""
    load_relapse_data()
"""
function load_relapse_data() 
    return @pipe DataFrame(XLSX.readtable("data/20210706_CBN_Wellness_Relapser.xlsx", "Sheet1")) |> 
        select(_, Not([
            :STUDYID, :TRTP, :PARAM, :PARAMCD, :AGE, :AGEU, :AGEGR1, :AGEGR1N,
            :RACE, :COUNTRY, :Comment, :ENRLFL, :FASFL, :SEX, :ETHNIC])) |> 
        rename(_, :SUBJID => :subject_id)
end 

""" 
    load_demographics()

    Load the demographics/baseline data. 
    
    For the socioeconomics data, we have the following codes: 
        EDULEVEL	Level of Education Attained
        EMPSTAT	    Employment Status
        HANDDOM	    Dominant Hand
        INCMLVL	    Income Level
        JOBCLAS	    Employee Job Class
        MARISTAT	Marital Status

    The `fa` (psychiatric history) file has the following codes: 
        OCCUR	Occurrence Indicator
        FAODTC	Onset Date
        FAONAGE	Age of Onset
        OCCUR	Occurrence Indicator
        DUR	    Duration
        FAODTC	Onset Date
        TRTREC	Treatment Received
        PRVMDEN	Number of Previous MDD Episodes
        VRFMDEN	Verified Number of Previous MDE

"""
function load_demographics()

    # Load age, sex and ethnicity as well as subject ids
    demog = @pipe CSV.read("data/dm.sas7bdat.csv", missingstring=["NA", ""], DataFrame)[:, 
            [:SUBJID, :USUBJID, :AGE, :SEX, :ETHNIC]
        ] |>
        rename(_, :SUBJID => :subject_id)  

    # Load the socioeconomics data and handedness
    sc = @pipe CSV.read("data/sc.sas7bdat.csv", missingstring=["NA", ""], DataFrame) |> 
        select(_, [:USUBJID, :SCTESTCD, :SCORRES]) |> 
        unstack(_, [:USUBJID], :SCTESTCD, :SCORRES)
    demog = leftjoin(demog, sc, on=:USUBJID) 

    # Psychiatric historical objects
    fa = @pipe CSV.read("data/fa.sas7bdat.csv", missingstring=["NA", ""], DataFrame)

    # Extract age of onset
    age_of_onset = @pipe filter(:FATEST => x -> x == "Age of Onset", fa)[:, [
            :USUBJID, :FAORRES
        ]] |> 
        rename(_, :FAORRES => :age_of_onset)
    demog = leftjoin(demog, age_of_onset, on=:USUBJID)

    # Extract family history of psychological illness (binary)
    famhx = @pipe fa[(fa.FATESTCD .== "OCCUR") .& (fa.FAOBJ .== "FAMILY HISTORY OF PSYCHOLOGICAL ILLNESS"), [:USUBJID, :FAORRES]] |> 
        rename(_, :FAORRES => :famhx)
    demog = leftjoin(demog, famhx, on=:USUBJID)

    # Extract number of previous major depressive episodes 
    # mde_num = @pipe fa[(fa.FAOBJ .== "MAJOR DEPRESSIVE DISORDER EPISODE") .& (fa.FATESTCD .== "PRVMDEN"), [:USUBJID, :FAORRES]] ## ORIGINAL, PRE COMPLETION OF MDE NUMBER
    # mde_num.mde_num = parse.(Int, mde_num.FAORRES) ## ORIGINAL, PRE COMPLETION OF MDE NUMBER
    mde_num = DataFrame(XLSX.readtable("data/Wellness MDE FINAL 20231101.xlsx", "Sheet1", infer_eltypes=true))[:, [:USUBJID, :MDE_NUM_FINAL]]
    mde_num_temp = Vector{Union{Missing, Int64}}(undef, length(mde_num.MDE_NUM_FINAL))
    mde_num_temp[mde_num.MDE_NUM_FINAL .!= 9999] .= mde_num.MDE_NUM_FINAL[mde_num.MDE_NUM_FINAL .!= 9999]
    mde_num.mde_num = mde_num_temp
    mde_num = select(mde_num, Not(:MDE_NUM_FINAL))
    demog = leftjoin(demog, mde_num[:,[:USUBJID, :mde_num]], on=:USUBJID)

    # Get the current MDE start date
    currmdestart = @pipe fa[(fa.FATEST .== "Onset Date") .& (fa.FAOBJ .== "MAJOR DEPRESSIVE DISORDER"), [:USUBJID, :FAORRES]] |> 
        rename(_, :FAORRES => :currmdestart)   
    demog = leftjoin(demog, currmdestart, on=:USUBJID)

    # Load comorbidities
    repl_checked(x) = replace(x, "CHECKED" => "Yes", "NOT CHECKED" => "No")
    mini = CSV.read("data/mini.csv", missingstring=["NA", ""], DataFrame)
    mini_rev = DataFrame(
        USUBJID = mini.USUBJID,
        panic_disorder_curr = repl_checked(mini.MINI0123),
        agoraphobia = repl_checked(mini.MINI0125),
        social_phobia = repl_checked(mini.MINI0126),
        ocd = repl_checked(mini.MINI0127),
        ptsd = repl_checked(mini.MINI0128),
        gad = repl_checked(mini.MINI0138),
        etoh = repl_checked(mini.MINI0129),
        drug = repl_checked(mini.MINI0130),
        lifetime_psychotic = repl_checked(mini.MINI0131),
        personality_disorder = repl_checked(mini.MINI0140),
        anorexia = repl_checked(mini.MINI0135),
        bulimia = repl_checked(mini.MINI0136),
        binge_eating = repl_checked(mini.MINI0137),
    )
    demog = leftjoin(demog, mini_rev, on=:USUBJID)

    # Load baseline scale scores
    #     madrs 
    madrs = @pipe CSV.read("data/madrs.csv", missingstring=["NA", ""], DataFrame) |> 
        select(_, Not([:QSALL, :MADR102A])) |> 
        filter(:VISIT => x -> x == "BASELINE", _) |> 
        rename(_, :MADRS111 => :madrsbl) |> 
        select(_, Not([:QSDTC, :VISIT]))
    demog = leftjoin(demog, madrs, on=:USUBJID)

    #    qids
    qids = @pipe CSV.read("data/qids.csv", missingstring=["NA", ""], DataFrame) |>
    filter(:VISIT => x -> x == "BASELINE", _)

    qids_mtx = qids[:, 4:end] |> Matrix 
    qids.qidsbl .= maximum(qids_mtx[:, 1:4], dims=2) .+ 
        qids_mtx[:, 5] .+ 
        maximum(qids_mtx[:, 6:9], dims=2) .+ 
        sum(qids_mtx[:, 10:14], dims=2) .+ 
        maximum(qids_mtx[:, 15:16], dims=2)

    qids = combine(groupby(sort(qids, :QSDTC), :USUBJID), 
        :QIDS0201 => first => :QIDS0201,
        :QIDS0202 => first => :QIDS0202,
        :QIDS0203 => first => :QIDS0203,
        :QIDS0204 => first => :QIDS0204,
        :QIDS0205 => first => :QIDS0205,
        :QIDS0206 => first => :QIDS0206,
        :QIDS0207 => first => :QIDS0207,
        :QIDS0208 => first => :QIDS0208,
        :QIDS0209 => first => :QIDS0209,
        :QIDS0210 => first => :QIDS0210,
        :QIDS0211 => first => :QIDS0211,
        :QIDS0212 => first => :QIDS0212,
        :QIDS0213 => first => :QIDS0213,
        :QIDS0214 => first => :QIDS0214,
        :QIDS0215 => first => :QIDS0215,
        :QIDS0216 => first => :QIDS0216,
        :qidsbl => first => :qidsbl)
    demog = leftjoin(demog, qids, on=:USUBJID)

    #    gad7
    gad7 = @pipe CSV.read("data/gad7.csv", missingstring=["NA", ""], DataFrame) |> 
        filter(:VISIT => x -> x == "BASELINE", _) |> 
        combine(groupby(_, :USUBJID), 
            :GAD0101 => first => :GAD0101,
            :GAD0102 => first => :GAD0102,
            :GAD0103 => first => :GAD0103,
            :GAD0104 => first => :GAD0104,
            :GAD0105 => first => :GAD0105,
            :GAD0106 => first => :GAD0106,
            :GAD0107 => first => :GAD0107
        )
    gad7.gad7bl = sum(gad7[:, 4:end] |> Matrix, dims=2) |> vec
    demog = leftjoin(demog, gad7, on=:USUBJID)

    #   cgi
    cgi = @pipe CSV.read("data/cgi.csv", missingstring=["NA", ""], DataFrame) |> 
        select(_, Not([:QSDTC, :QSALL])) |>
        filter(:VISIT => x -> !ismissing(x), _) |>
        filter(:VISIT => x -> x ∈ ["SCREENING", "BASELINE"], _) |> 
        combine(groupby(_, :USUBJID), :CGI0101 => mean => :cgibl)
    demog = leftjoin(demog, cgi, on=:USUBJID)

    #   sds
    sds = @pipe CSV.read("data/sds.csv", missingstring=["NA", ""], DataFrame) |> 
    filter(:VISIT => x -> x == "BASELINE", _) |> 
    select(_, Not([:SDS0101A, :SDS0104, :SDS0105]))
    sds_mtx = sds[:, 4:end] |> Matrix
    miss_idx = @pipe sds_mtx |> findall(ismissing, _)
    sds_mtx[miss_idx] .= 0
    sds.sdsbl = sum(sds_mtx, dims=2) |> vec
    demog = leftjoin(demog, sds, on=:USUBJID)

    #   leaps
    leaps = @pipe CSV.read("data/leaps.csv", missingstring=["NA", ""], DataFrame) |> 
    filter(:VISIT => x -> x == "BASELINE", _) |> 
    select(_, Not([:LEAP901, :LEAP901A, :LEAP902, :LEAP903])) |> 
    _[.!all(ismissing.(_[:, 4:end] |> Matrix), dims=2) |> vec, :]
    leaps.leapsbl = sum(leaps[:, 4:end] |> Matrix, dims=2) |> vec
    demog = leftjoin(demog, select(leaps, Not([:QSDTC, :VISIT])), on=:USUBJID)

    #  qlesq
    qlesq = @pipe CSV.read("data/qlesq.csv", missingstring=["NA", ""], DataFrame) |> 
    filter(:VISIT => x -> x == "BASELINE", _) |>
    _[:, 1:17] 
    qlestot = @pipe qlesq[:, 4:end] |> Matrix |> sum(_, dims=2) |> vec 
    qlesq.qlesqbl = (qlestot .- 14)/56
    demog = leftjoin(demog, select(qlesq, Not([:VISIT, :QSDTC])), on=:USUBJID)

    #    whoqol  [TODO: FIGURE OUT SCORING ]
    # whoqol = @pipe CSV.read("data/whoqol.csv", missingstring=["NA", ""], DataFrame) |> 
    # filter(:VISIT => x -> x == "BASELINE", _) #|> 
    # select(_, Not([:QSDTC, :VISIT])) |> 
    # _[.!all(ismissing.(_[:, 4:end] |> Matrix), dims=2) |> vec, :]

    #    BRIAN 
    brian = @pipe CSV.read("data/brian.csv", missingstring=["NA", ""], DataFrame) |> 
    filter(:VISIT => x -> x == "BASELINE", _) |> 
    select(_, Not([:QSDTC, :VISIT])) 
    brian.brianbl = @pipe brian[:, 4:end] |> Matrix |> sum(_, dims=2) |> vec
    demog = leftjoin(demog, brian, on=:USUBJID)

    return demog
end

"""
    load_medications()

    Load the medications data. 
    We are only using those marked for depression
"""
function load_medications()

    cm = @pipe DataFrame(load("data/cm.sas7bdat")) |> 
        select(_, Not([
            :CMSEQ, :DOMAIN, :CMSPID, :CMDECOD, :CMSCAT, 
            :CMPRESP, :CMOCCUR, :CMCLAS, :CMINDC, :CMCLASCD, :CMDOSE, :CMDOSU,
            :CMDOSTXT, :CMDOSFRQ, :CMROUTE, :CMADJ, :EPOCH, :CMSTDTC, :CMENDTC,
            :CMSTDY, :CMENDY
            ])) |> 
        filter(:CMCAT => x -> x ∈ ["ANTIDEPRESSANT TREATMENT", "ANTIDEPRESSANT TREATMENT HISTORY"]) |> 
        filter(:CMEVINTX => x -> x != "LIFETIME", _)


    cm.CURRENT .= "No"
    cm.CURRENT[cm.CMENRF .== "ONGOING"] .= "Yes"
    cm.CURRENT[cm.CMEVINTX .== "CURRENT"] .= "Yes"
    cm = filter(:CURRENT => x -> x == "Yes", cm)
    cm = select(cm, Not([:CMENRF, :CMEVINTX, :CMCAT, :CMGRPID]))

    cm.CMTRT = replace(cm.CMTRT, 
        "ABILIFY" => "ARIPIPRAZOLE",  
        "AMITRIPTYLINE (ELAVIL, ENDEP)" => "AMITRIPTYLINE", 
        "APO-MOCLOBEMIDE" => "MOCLOBEMIDE", 
        "ARIPIPRAZOLE" => "ARIPIPRAZOLE", 
        "ARIPIPRAZOLE (ABILIFY)" => "ARIPIPRAZOLE",
        "BEHAVIORAL ACTIVATION THERAPY" => "PSYCHOTHERAPY",
        "BUPROPION (SLOW RELEASE)" => "BUPROPION",
        "BUPROPION (SLOW RELEASING)" => "BUPROPION",
        "BUPROPION (SR)" => "BUPROPION",
        "BUPROPION (WELLBUTRIN, WELLBUTRIN SR, WELLBUTRIN XL)" => "BUPROPION",
        "BUPROPION (XL)" => "BUPROPION",
        "BUPROPION XL" => "BUPROPION",
        "CIPRALEX" => "ESCITALOPRAM",
        "CITALOPRAM (CELEXA)" => "CITALOPRAM",
        "COGNITIVE BEHAVIORAL THERAPY" => "PSYCHOTHERAPY",
        "CYMBALTA" => "DULOXETINE",
        "DESRENLAFAXINE (PRISTIQ)" => "DESVENLAFAXINE",
        "DESVENLAFAXINE (PRISTIQ)" => "DESVENLAFAXINE",
        "DEXEDRINE" => "PSYCHOSTIMULANT",
        "DIVALPROEX" => "VPA",
        "DULOXETINE (CYMBALTA)" => "DULOXETINE",
        "ECT TYPE UNKNOWN" => "ECT",
        "EFFEXOR" => "VENLAFAXINE",
        "ELEFEXOR" => "VENLAFAXINE",
        "ESCITALOPRAM (LEXAPRO, CIPRALEX)" => "ESCITALOPRAM",
        "FETZIMA" => "LEVOMILNACIPRAN",
        "FLUOXETINE (PROZAC)" => "FLUOXETINE",
        "FLUVOXAMINE (LUVOX)" => "FLUVOXAMINE",
        "IMIPRAMINE (TOFRANIL)" => "IMIPRAMINE",
        "INTERPERSONAL THERAPY" => "PSYCHOTHERAPY",
        "LAMOTREIGENE" => "LAMOTRIGINE",
        "LAMOTRIGENE" => "LAMOTRIGINE",
        "LAMOTRIGINE" => "LAMOTRIGINE",
        "LAMOTRIGINE (LAMICTAL)" => "LAMOTRIGINE",
        "LATUDA" => "LURASIDONE",
        "LITHIUM (AS AN AUGMENTING AGENT FOR MDD)" => "LITHIUM",
        "METHYLPHENIDATE" => "PSYCHOSTIMULANT",
        "MIRTAZAPINE (REMERON)" => "MIRTAZAPINE",
        "MODAFANIL" => "MODAFINIL",
        "NORTRIPTYLINE (PAMELOR, AVENTYL)" => "NORTRIPTYLINE",
        "OLANZAPINE (ZYPREXA)" => "OLANZAPINE",
        "PARNATE" => "TRANYLCYPROMINE",
        "PAROXETINE (PAXIL)" => "PAROXETINE",
        "PAROXETINE CR (PAXIL CR)" => "PAROXETINE",
        "PAXIL" => "PAROXETINE",
        "PRAMIPEXOLE (MIRAPEX)" => "PRAMIPEXOLE",
        "PRISTIG" => "DESVENLAFAXINE",
        "PRISTIQ" => "DESVENLAFAXINE",
        "PRISTRIG" => "DESVENLAFAXINE",
        "PROZAC" => "FLUOXETINE",
        "QUETIAPINE (SEROQUEL)" => "QUETIAPINE",
        "REMERON (MIRTAZAPINE)" => "MIRTAZAPINE",
        "REXULTI" => "BREXPIPRAZOLE",
        "RISPERIDONE (RISPERDAL)" => "RISPERIDONE",
        "SECTRALINE" => "SERTRALINE",
        "SEROQUEL" => "QUETIAPINE",
        "SERTRALINE (ZOLOFT)" => "SERTRALINE",
        "TRAZADONE" => "TRAZODONE",
        "TRAZODONE (DESYREL)" => "TRAZODONE",
        "TRINTELLIX" => "VORTIOXETINE",
        "VALPROIC ACID" => "VPA",
        "VARTIOXETINE" => "VORTIOXETINE",
        "VELAFAXINE (XR)" => "VENLAFAXINE",
        "VENLAFAXINE" => "VENLAFAXINE",
        "VENLAFAXINE (EFFEXOR AND EFFEXOR XR)" => "VENLAFAXINE",
        "VENLAFAXINE (EFFEXOR)" => "VENLAFAXINE",
        "VENLAFAXINE (XR)" => "VENLAFAXINE",
        "VENLAFAXINE XR" => "VENLAFAXINE",
        "VORTIOXEFINE" => "VORTIOXETINE",
        "VORTIOXETINE" => "VORTIOXETINE",
        "WELBUTRIN" => "BUPROPION",
        "WELBUTRIN/BUPROPION" => "BUPROPION",
        "WELLBUTRIN" => "BUPROPION",
        "WELLBUTRIN SR" => "BUPROPION",
        "WELLBUTRIN XL" => "BUPROPION",
        "ZOLOFT" => "SERTRALINE",
        "ZOLPIDEM" => "Z-DRUG",
        "ZOPICLONE" => "Z-DRUG",
        "ZOPIDOL" => "Z-DRUG"
    )

    meds = DataFrame()
    cm_subjects = unique(cm.USUBJID)
    for s ∈ cm_subjects
        meds_s = cm[cm.USUBJID .== s, :CMTRT] |> unique 
        n_meds = length(meds_s)
        meds = vcat(meds, 
            DataFrame(
                USUBJID=s, 
                meds=join(meds_s, ";"), 
                n_meds=n_meds))
    end

    return meds
end

"""
    load_and_preprocess_clinical_data()

    Loads and preprocesses clinical data and preprocess it. It also gets other data ready 
    for various analyses 
"""
function load_and_preprocess_clinical_data()
    
    # Convert all the SAS files to CSVs
    sas_to_csv()

    # Extract the questionnaires to CSVs
    extract_questionnaires_to_csv()

    # Load the ID's and the relapse data 
    subject_list = load_subject_ids() # 98 rows
    relapse_data = load_relapse_data() # 97 rows 
    clin_data = innerjoin(subject_list, relapse_data, on=:subject_id) # 97 rows

    # Load demographics and baseline information 
    demog = load_demographics() # 103 rows 
    clin_data = leftjoin(clin_data, demog, on=[:subject_id, :USUBJID]) # 97 rows
    
    # Load the medications data
    meds = load_medications() # 97 rows
    clin_data = leftjoin(clin_data, meds, on=:USUBJID) # 97 rows

    # Write the clinical data to a CSV
    CSV.write("data/clinical_data.csv", clin_data)

    # Reload for better formatting of columns 
    clin_data = CSV.read("data/clinical_data.csv", missingstring=["NA", ""], DataFrame)

    # Create relapse column 
    clin_data.relapse .= 0
    #clin_data.relapse[clin_data.Relapse_Status .== "Relapse"] .= 1 # Relapse_Status is the pre-adjudication status?
    clin_data.relapse[clin_data.CNSR .== 0] .= 1

    # Compute current MDE duration 
    currmdedur_days = clin_data.STARTDT .- clin_data.currmdestart
    currmdedur = Vector{Union{Missing, Int64}}(missing, length(clin_data.STARTDT))
    currmdedur[.!ismissing.(currmdedur_days)] .= currmdedur_days[.!ismissing.(currmdedur_days)] .|> Dates.value
    clin_data.currmdedur = currmdedur

    # Re-write the clinical data to a CSV now that dates have been handled
    CSV.write("data/clinical_data.csv", clin_data)

    # Save subset for creation of table1 
    tableonedata = @pipe clin_data[:, [
        "relapse", 
        "SITEID", 
        "AGE", "SEX", "ETHNIC", "EDULEVEL", "EMPSTAT", "HANDDOM", "INCMLVL", "JOBCLAS", "MARISTAT", 
        "age_of_onset", "mde_num", "currmdedur", "n_meds", 
        "panic_disorder_curr", "agoraphobia", "social_phobia", "ocd", "ptsd", "gad", "etoh", "drug", 
        "lifetime_psychotic", "personality_disorder", "anorexia", "bulimia", "binge_eating", "famhx", 
        "madrsbl", "qidsbl", "gad7bl", "cgibl", "sdsbl", "leapsbl", "qlesqbl", "brianbl"]] |> 
        CSV.write("data/wellness-tableonedata.csv", _)

    return clin_data
end

data = load_and_preprocess_clinical_data()


