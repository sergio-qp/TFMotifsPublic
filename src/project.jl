using DataFrames, DataFramesMeta, CSV, Glob, GenomicFeatures, Plots, StatsPlots, Measures
using BioSequences, FASTX, StatsBase, Distributions
using MotifScanner, GenomicIntersections

#global setpalette = :tol_sunset
theme(:wong)

# this is the begenning of the convenience function section----------------------------------------------------
function showwide(table)
c = ENV["COLUMNS"]
ENV["COLUMNS"] = "10000"
display(table)
ENV["COLUMNS"] = c;
nothing;
end

function regexmaker(vector)
    outp=Regex(join(vector,"|"))
    outp
end

macro Name(arg)
    string(arg)
end

# convert every instance of NKX2-8 into NKX2-2 just for the innerjoin later, so its important to remember u have done this
# this is very unsustainable, bad
function nkxchange(tf)
    if tf == "NKX2-8"
    tf = "NKX2-2"
    end
    tf
end

function annotatecoll!(tableA, tableB, colb_in, cola_out ; default=0)
    ivA = GenomicIntersections.intervals(tableA)::IntervalCollection{Int64}
    ivB = GenomicIntersections.intervals(tableB)::IntervalCollection{Int64}
    tableA[!, cola_out] = fill(default, size(tableA, 1))
    for (ia, ib) in eachoverlap(ivA, ivB)
        tableA[metadata(ia), cola_out] = tableB[metadata(ib), colb_in]
    end
    tableA
end

function pairlistmaker(tfvec)
        list = Vector{String}[]
    for i in 1:length(tfvec)
        for n in 1:(length(tfvec)-1)
            push!(list, [tfvec[i],  tfvec[n+1]])
        end
    end
    list=filter(x -> length(x) > 1,unique.(unique(sort(sort.(list)))))
    list
end
# i can make these two into one function by just making the vector whatever type the dict it storing, id have to figure out the function for that tho
function dicttovec(dict::Dict{String,DataFrame},keys)
    vector=NamedTuple[]
    for key in keys
        push!(vector, dict[key])
    end
    vector
end

function dicttovec(dict::Dict{String,Plots.Plot{Plots.GRBackend}},keys)
    vector=Plots.Plot{Plots.GRBackend}[]
    for key in keys
        push!(vector, dict[key])
    end
    vector
end

function stackdataframe(df, emptyvector=Int64[])
    collect(eachcol(df))
    for list in vcat(collect(eachcol(df)))
        append!(emptyvector,list) 
    end
    emptyvector
end

function bettersubset(df, column, conditional)
    output = filter(x -> conditional.(getindex(x,column)),df)
    output
end

function dictplotter(plotdict)
    for key in sort(collect(keys(plotdict)))
        display(plot(plotdict[key]))
    end
end

# it fucking worked, i dont believe it
function rangesubset(x,y)
    function rangesubset(inpt)
        if inpt >= x && inpt <= y
            true
        else
            false
        end
    end
end

function ismorer(x)
    function ismorer(inpt)
        if inpt >= x
            true
        else
            false
        end
    end
end


# this is the end of the convenience function section----------------------------------------------------

peakintersectFC(labels::Vector{T}, peakdict::Dict{V, DataFrame}) where {T, V} = peakintersectFC(labels, [peakdict[l] for l in labels])
function peakintersectFC(labels::Vector{T}, peaks::Vector{DataFrame}) where {T}
    allpeaks = reduce(vcat, peaks)

    sort!(allpeaks, [:chrom, :start, :stop, :Origin]);
    chroms    = allpeaks.chrom #::Vector{String}
    locations = [(a+1):b for (a, b) in zip(allpeaks.start, allpeaks.stop)] #::Vector{UnitRange{Int}}
    allpeaks.Group = overlappinglocations(chroms, locations);


    combpeaks = combine(groupby(allpeaks, :Group),
                    :chrom => first => :chrom,
                    :start => minimum => :start,
                    :stop => maximum => :stop,
                    nrow => :TotalPeaks,
                    :score => mean => :score,
                    :Origin => Ref => :Origin,
                    :FC => mean => :FC)

    membership =  DataFrame([1(l .∈ combpeaks.Origin) for l in labels], Symbol.(labels))

    [combpeaks[!, Not([:Group, :Origin])] membership]
end

# make 1 based inclusive, add some names and unique snp ID column
function loadgwas!(gwas)
    gwas.start .+= 1;
    gwas.ID = 1:size(gwas, 1)
    loc = split.(String.(gwas.name), r"[_]")
    gwas.ref = getindex.(loc,2)
    gwas.alt = getindex.(loc,3)
    rename!(gwas, :name => :snp)
    gwas
end

# filters the gvat set down to just the tfs ur interested in, adds unique snp ID col and performs nkxchange (i dont know if we still need that, but i think we do, there's 
# probably a way to generalize it by giving it specific tfs u want to substitute and what u want to substitute them with in a vector)
function loadgvat(rawgvat,tfvector)
    tfreg = regexmaker(tfvector)
    gvat = @subset(rawgvat, occursin.(tfreg,:TF))
    loc = split.(String.(gvat.snp), r"[_]")
    gvat.chrom = first.(loc)
    gvat.start = parse.(Int, getindex.(loc,2))
    gvat.stop = parse.(Int, getindex.(loc,2))
    gvat.ref = getindex.(loc,3)
    gvat.alt = getindex.(loc,4)
    insertcols!(gvat, 4, :alt_auc => gvat.oligo_auc .- gvat.pbs)
    sort!(gvat, [:chrom, :start, :stop])
    uniquefilt = unique(gvat, :snp)
    uniquefilt.ID = 1:size(uniquefilt, 1)
    inefficiency= innerjoin(gvat, uniquefilt, on=:snp, makeunique=true)
    gvat.ID = inefficiency.ID

    gvat.TF = nkxchange.(gvat.TF);
    gvat
end

# this gets the seqs of the peaks of one tf peak file and appends them to the same table (its peakdict table)
function identifyseqs!(peakdictfile, genomefile)
    tfset = peakdictfile
    file  = genomefile


    reader = open(FASTA.Reader, file, index=string(file, ".fai"))
    seqs = LongSequence{DNAAlphabet{4}}[]
    fields = split.(string.(tfset.chrom, "+", tfset.start, "+", tfset.stop), r"[+]")
    chroms = first.(fields)
    iv(a,b) = a:b
    w = 0
    locs = iv.(parse.(Int, getindex.(fields, 2)) .+ w,parse.(Int, getindex.(fields, 3)) .- w)
    for (c, l) in zip(chroms, locs)
        seq = FASTA.extract(reader, DNAAlphabet{4}(), c, l)
        push!(seqs, seq)
    end
    close(reader)
    peakdictfile.seqs = seqs
    peakdictfile
end


function loadreplicates(filelocation, genomefile, tfvector)

    #does not include any shuffling functions, can prolly add em on with an if statement asking if u want a shuffled table
    #original loadreplicates
    #-----------------------------------------------------------------------------------------------------------
    files = glob("macs_*.narrowPeak", filelocation)
    filter!(f -> !occursin(r"H[23]", f), files)
    samples = replace.(basename.(files), r"macs_|_peaks.narrowPeak" => "")
    factor = replace.(replace.(samples, r"_HI_[0-9]*" => ""), r"_" => "-")
    meta = DataFrame(Study="Islets", Sample=samples, Factor=factor, PeakFile=files)
    meta.NumPeaks = countlines.(meta.PeakFile)
    peaks = CSV.read.(meta.PeakFile, DataFrame, skipto=1, header=[:chrom, :start, :stop, :name, :score, :strand, :FC, :nlog10p, :nlog10q, :summit]);
    for (p, l) in zip(peaks, meta.Sample)
    p[!, :Origin] .= l;
        p.chrom = string.(p.chrom)
    end

    meta.Index = 1:size(meta, 1)
    peakdict = Dict{String, DataFrame}()
    for sdf in groupby(meta, :Factor)
        sdf.Factor[1]=replace(sdf.Factor[1], r"_" => "-")
        peaktable = peakintersectFC(string.(sdf.Sample), peaks[sdf.Index])
        peaktable[!, :Origin] .= sdf.Factor[1]
        # display("============================================")
        peakdict[sdf.Factor[1]] = @subset(peaktable, :TotalPeaks .== size(sdf, 1))[!, Not(r"_HI_")]
        peakdict[sdf.Factor[1]][:,:ind] = collect(1:size(peakdict[sdf.Factor[1]], 1)) 
        #this next line filters out all the weird chroms
        peakdict[sdf.Factor[1]] = @subset(peakdict[sdf.Factor[1]], .!occursin.(r"_|chrY",:chrom))

    end

    tfs = sort(unique(meta.Factor))
    tfcombo = peakintersectFC(tfs, peakdict)
    tfcombo[!,:peakID] = 1:length(tfcombo.chrom)
    # countintersection!(tfcombo, gwas, :gwasSNPCount)
    # countintersection!(tfcombo, unique(gvat, [:snp]), :gvatSNPCount)
    #-------------------------------------------------------------------------------------------------------------
    # this does an identifyseq for every tf in a vector, just to get them all onto the peakdict file
    for tf in tfvector
        identifyseqs!(peakdict[tf],genomefile)
    end

    peakdict, tfcombo
end

function mymotifmaker(motiffilelocation, tfvector)
    tfreg = regexmaker(tfvector)
    mymotifs = Dict{String, NamedTuple}()
     # make motif index and clean them up for use with pkg
     motifmeta, motifs = loadtransfac(motiffilelocation);
     motifmeta.Index = 1:size(motifmeta, 1)
     motifmeta.MotifName = uppercase.(motifmeta.MotifName)
     motiffilt = @subset(motifmeta, occursin.(tfreg,:MotifName))
     #  index = hcat(String.(last.(split.((motiffilt.MotifName), r"[.]"))),motiffilt.Index) i have no idea what this does or why its here
     for motif in motifs[motiffilt.Index]
        i=1
        if haskey(mymotifs, motif.name)
            while haskey(mymotifs,string(uppercase(motif.name), "_", string(i)))
                i+=1
            end
            mymotifs[string(uppercase(motif.name), "_", string(i))] = motif
            
        else
            mymotifs[uppercase(motif.name)] = motif
        end
     end
     mymotifs
end

#this is such a scuffed function, but it basically just scans the insides of the seqs for possible instances of the motif
#and then keeps only the ones that have a likelihood score above a set number 
#THERE HAS BEEN A MODIFICATION, :start and :stop are now the peak coords and :motifstart/stop are the abs_coords, this is to do an annotatecoll! on their peaks
#the indexing has changed, it used to be an incrementing variable, but now its just what gets given by the eachrow/seq thing, which has fixed the issue
#this version works by taking all the points above a certain cutoff, i want to try and get a version that works off of maximal points
function singlemotifscan(set, motif, cutoff, tf)
    df=@subset(scanmotstats(motif,set.seqs[1]), :Motif .== 0)
    insertcols!(df, 4, :peakind => 0)
    insertcols!(df, 4, :motifstop => 0)
    insertcols!(df, 4, :motifstart => 0)
    insertcols!(df, 2, :chrom => "")
    for seq in eachrow(set)
        seqdf=@subset(scanmotstats(motif,seq.seqs), :prmax .> cutoff)
        # ive replaced the manual indexing with just copying the in from the seq line
        seqdf[!,:peakind] .= seq.ind
        seqdf[!,:chrom] .= seq.chrom
        seqdf[!,:motifstart] .= seqdf[!,:start] .+ seq.start .- 1
        seqdf[!,:motifstop] .= seqdf[!,:stop] .+ seq.start .- 1
        seqdf[!,:start] .= seq.start
        seqdf[!,:stop] .= seq.stop
        append!(df, seqdf)
    end
    rename!(df, string.(tf, names(df)))
    df
end


# one function to automate the singlescan motifs with tf and motif vectors
function multiplemotifscan(peakdict,tfvector,motifvector,cutoff)

    motifscandict = Dict{String, DataFrame}()
    for (tf, motif) in zip(tfvector, motifvector)
        motifscandict[tf] = singlemotifscan(peakdict[tf], motif, cutoff, tf)
    end
    motifscandict
end

#this is a version that works off of only the maximal points
function singlemaxmotifscan(set, motif, tf)
    df=@subset(scanmotstats(motif,set.seqs[1]), :Motif .== 0)
    insertcols!(df, 4, :peakind => 0)
    insertcols!(df, 4, :motifstop => 0)
    insertcols!(df, 4, :motifstart => 0)
    insertcols!(df, 4, :cdist => 0)
    insertcols!(df, 2, :chrom => "")
    for seq in eachrow(set)
        seqdf=scanmotstats(motif,seq.seqs)
        seqdf[!,:peakind] .= seq.ind
        seqdf[!,:chrom] .= seq.chrom
        seqdf[!,:motifstart] .= seqdf[!,:start] .+ seq.start .- 1
        seqdf[!,:motifstop] .= seqdf[!,:stop] .+ seq.start .- 1
        seqdf[!,:start] .= seq.start
        seqdf[!,:stop] .= seq.stop
        if !isempty(seqdf)
            maxval=maximum(seqdf.prmax)
            maxrows=@subset(seqdf, :prmax .== maxval)
            maxrows[!,:cdist] = abs.(div.(maxrows[:,:start] .+ maxrows[:,:stop],2) .- div.(maxrows[:,:start] .+ maxrows[:,:stop],2))
            mind=minimum(maxrows.cdist)
            maxrows=@subset(maxrows, :cdist .== mind)
            push!(df, maxrows[1,:])
        end
    end
    rename!(df, string.(tf, names(df)))
    df
end

# the relative multiple scan function
function multiplemaxmotifscan(peakdict,tfvector,motifvector)

    maxmotifscandict = Dict{String, DataFrame}()
    for (tf, motif) in zip(tfvector, motifvector)
        maxmotifscandict[tf] = singlemaxmotifscan(peakdict[tf], motif, tf)
    end
    maxmotifscandict
end

# one function to automate annotatecoll! executions to match respective peaks for a pair
function peakmatch(peakdict, tfvector)
    #this does a flat peakintersectFC between the peakdict files you give it, then subsets it to get just the cases where both have peaks, like a bedtools merge(?)
    peakcombo = peakintersectFC(tfvector, peakdict)
    total = vec(sum(Matrix(peakcombo[!, Symbol.(tfvector)]), dims=2))
    ind = total .> (length(tfvector) - 1)
    df = peakcombo[ind, :]

    #then this performs an annotatecoll! for loop, watch out tho, cuz its a ! function
    for tf in tfvector
        annotatecoll!(df, peakdict[tf], :ind, Symbol(tf,"peakind"))
    end
    df.ID = collect(1:size(df, 1))
    df
end

function scanmachine(snpset,mymotifs,genomelocation)
   
    uniquefilt = unique(snpset, [:snp, :chrom, :start, :stop])
    seqs = loadrefseqs(uniquefilt, genomelocation, 25);
    seqs.ID = uniquefilt.ID
    scanned = motifscanall(seqs, mymotifs)
    scanned.MotifName = uppercase.(scanned.MotifName);
    scanned = innerjoin(scanned,seqs, on =:ID);

    if sum(occursin.("TF",names(snpset))) > 0
        outputdf = innerjoin(snpset,scanned, on =[:ID, :TF => :MotifName])
    else
        outputdf = innerjoin(snpset,scanned, on =:ID)
    end
    outputdf
end

function repeatscanmachine(snpsetvec,snpnames,mymotifs,genomelocation)
    snpscandict = Dict{String, DataFrame}()
    for (snpset,snpname) in zip(snpsetvec, snpnames)
        snpscandict[snpname] = scanmachine(snpset,mymotifs,genomelocation)
    end
    snpscandict
end

function snpidentifier!(maindf, dfvector, dfnames, tfvector)
    tfreg=regexmaker(tfvector)
    for (df, dfname) in zip(dfvector,dfnames)
        if sum(occursin.("TF",names(df))) > 0
            df=@subset(df, occursin.(tfreg,:TF))
            annotatecoll!(maindf, df, :ID, Symbol(dfname,"ID"))
        else
            annotatecoll!(maindf, df, :ID, Symbol(dfname,"ID"))
        end
    end
    maindf
end


#one function to perform the innerjoins that will match corresponding tf peaks
function innerjoiner(maindf, motifscandict, tfvector)
    for tf in tfvector
        maindf = innerjoin(maindf,motifscandict[tf], on=Symbol(tf,"peakind"))
        # return maindf
    end
    maindf
end



# well, the sanity checker worked very well, surprisingly, since it showed that the problem persists
# there is still tons of instances of peaks being matched to the place they dont belong
# its interesting also that there are cases where the peak matches neiter of the two tfs, you'd think itd at least match one of them
# the problem has now been fixed, it was with indexing assigning in singlemotifscan
function sanitychecker(df, tfvector)
    startbitvec=BitVector[]
    stopbitvec=BitVector[]

    for tf in tfvector
        tfdf = df[!, [Symbol.(tf,"start"),Symbol.(tf,"stop")]]
        startbitvec = df.start .> tfdf[!,1]
        stopbitvec = df.stop .< tfdf[!,2]
    end
    @show sum(startbitvec), sum(stopbitvec)
    # @show stopbitvec
end

# gets the motif distances based off of the middle of where the motif is
# also added a line to make the "cross" column 
function motifdistancer!(tf1,tf2,df)
    df[:,Symbol("motifdistance",tf1,"_",tf2)] .= div.(df[:,Symbol(tf1, "motifstart")] .+ df[:,Symbol(tf1, "motifstop")], 2) .- div.(df[:,Symbol(tf2, "motifstart")] .+ df[:,Symbol(tf2, "motifstop")], 2)
    df[:,:strandrelation] .= ifelse.(df[:,Symbol(tf1,"strand")] .== df[:,Symbol(tf2,"strand")],"same","opposite")
end

function motifdistancerloop!(matchingpeaks, tfvector; centertfs=nothing)
    if centertfs == nothing
        list = Vector{String}[]
        for i in 1:length(tfvector)
            for n in 1:(length(tfvector)-1)
                push!(list, [tfvector[i],  tfvector[n+1]])
            end
        end
        list=filter(x -> length(x) > 1,unique.(unique(sort(sort.(list)))))
        list

        for pair in list
            motifdistancer!(pair[1],pair[2],matchingpeaks)
        end
    else
        list = Vector{String}[]
        for ctf in centertfs
            for dtf in tfvector
                push!(list, [ctf,  dtf])
            end
        end
        list = filter(x -> x[1] != x[2],list)
        list

        for pair in list
            motifdistancer!(pair[1],pair[2],matchingpeaks)
        end
    end
end

function snpchecker!(df, tfvector)
    for tf in tfvector
        df[!,Symbol(tf,"snpcheck")] .= (df[!,:snploc] .< df[:,Symbol(tf, "motifstop")]) .& (df[!,:snploc] .> df[:,Symbol(tf, "motifstart")])
    end
    df
end
# makes sure the resulting distances are smaller than the peaks they're in
function distancechecker(df)
    cols = occursin.("motifdistance",names(df))
    lengthbitvec=BitVector[]
    errorsum=0
    tdf = (Matrix(df[!, cols]))
    # tdf = Matrix(df[!, cols])

    # return size(tdf, 2)
    for col in 1:size(tdf, 2)
        lengthbitvec = tdf[:,col] .> (df.stop .- df.start)
        errorsum = sum(lengthbitvec)
    end
    @show errorsum    
end

#you should be inputting here the matchingpeaks table and snp file info
# this is a bit of a bottleneck(?) for us, since its innerjoining everything, which may cut out cases where a pair of tfs show up even if every other tf doesnt
function snpdictmaker(maindf, dfvector, dfnames, tfvector, scorecolvec, cutoffvec)
    snpdict = Dict{String, DataFrame}()
    tfreg=regexmaker(tfvector)
    for (df, dfname, scorecol, cutoff) in zip(dfvector,dfnames, scorecolvec, cutoffvec)
        if sum(occursin.("TF",names(df))) > 0
            df=@subset(df, occursin.(tfreg,:TF))
        end
        # return df
        df = rename(df, Dict(:start => "snploc"))
        snpdict[dfname] = innerjoin(maindf,df[!,Not([:chrom,:stop])], on= Symbol(dfname,"ID") => :ID)
        # gonna need to make some changes so that the disrupted score isnt hard coded and for gvat
        snpdict[dfname][!,:disrupted] = snpdict[dfname][:,scorecol] .> cutoff
        snpdict[dfname][!,:disrupted] = string.(snpdict[dfname][!,:disrupted])
        snpdict[dfname][!,:disrupted] = replace.(snpdict[dfname][!,:disrupted], "true" => "disrupted")
        snpdict[dfname][!,:disrupted] = replace.(snpdict[dfname][!,:disrupted], "false" => "normal")
        # snpchecker!(snpdict[dfname],tfvector)
    end
    snpdict
end

# we want it to have a big if statement at the start that makes it so that it either checks for the presence or absence of certain scoredmotifs
# this is so that for something like 1, where u want to get peaks with no snps, you check for the absence of ones rather than the presence of a zero in both
# thinking about it, when you check for individual snps, you're gonna want a presence of one and an absence of 1s from the other
# to check for absence, you just check to make sure they are greater than 0
# 3 states, absence 0 is presence for both, absence 1 is 1st is presence, second is absence, absence 2 is absence of both
# presence is if x > 1, absence is if x == 0
# i didnt figure out how to the the sorting thing he mentioned, its not optimal but it works
function simplesort(dataset, groupcol, checkcols, checkvalues; absence=0, hardpresence=false)
    gdf = groupby(dataset, groupcol)
    outputdf = @subset(dataset, :chrom .== "dog")

    # this works but is dogshit, lets see if we can improve on it
    # nm, not worth the time
    if absence == 0
        if hardpresence == false
            for df in gdf
                if (sum(df[:,checkcols[1]] .== checkvalues[1]) > 1) .& (sum(df[:,checkcols[2]].== checkvalues[2]) > 1)
                    append!(outputdf,df)
                end
            end
        else
            for df in gdf
                if (sum(df[:,checkcols[1]] .== checkvalues[1]) == size(df,1)) .& (sum(df[:,checkcols[2]].== checkvalues[2]) == size(df,1))
                    append!(outputdf,df)
                end
            end
        end        
    elseif absence == 1
        for df in gdf
            if (sum(df[:,checkcols[1]] .== checkvalues[1]) > 1) .& (sum(df[:,checkcols[2]].== checkvalues[2]) == 0)
                append!(outputdf,df)
            end
        end
    elseif absence == 2
        for df in gdf
            if (sum(df[:,checkcols[1]] .== checkvalues[1]) == 0) .& (sum(df[:,checkcols[2]].== checkvalues[2]) == 0)
                append!(outputdf,df)
            end
        end
    end


    outputdf
end


# we now want a function that will perform all the steps for an analysis like we just did given a set of tfs, peakfiles, a main tf/tfvec from which distances will be calculated,
# and a snp set or set of snpsets, inputs should be the peakfiles, genomefile, motiffile, tfvector, maximal or above certain cutoff, cutoff=x, [snpset list], [snp zscore cols],
# [snpset names], 

#do this thing to check if the mymotifs is a dict or array
# sumdiv(x::Vector{T}) where {T} = sumdiv(sum(x))
# sumdiv(total) = total/2

# sumdiv([1, 2, 34])
# sumdiv(10)

# function multiplemaxmotifscan(...)
#    ### works with vectors
# end

# function multiplemaxmotifscan(peakdict, tfvector, motifdict::Dict{T, V}) where {T, V}

#    ### convert motifdict to vector 
#    motifscandict=multiplemaxmotifscan(peakdict,tfvector,motifvec)
    
# end

function identifymotifs(motiffile, tfvector, snpsetvec, peakdict, maxonly, cutoff, snpnamevec, centertfs)
    # makes the mymotifs dict
    mymotifs = mymotifmaker(motiffile, tfvector)
    # makes the mymotifs vector
    motifvec=NamedTuple[]
    for tf in tfvector
        push!(motifvec, mymotifs[tf])
    end
    # annotatecol magic identifies which peaks of a given set overlap and where
    pairedpeaks=peakmatch(peakdict,tfvector)
    # takes the now matched peaks and performs similar magic to match those peaks with snps from a set of snp sets
    snpidentifier!(pairedpeaks, snpsetvec, snpnamevec, tfvector)

    # performs independent scanmotstats on tfs from a tfvector. These are important because...? scanmachine returns motif distances too, better ones, in fact, but we have
    # to mess with them a lil
    if maxonly==true
        motifscandict=multiplemaxmotifscan(peakdict,tfvector,motifvec)
    else
        motifscandict=multiplemotifscan(peakdict,tfvector,motifvec,cutoff)
    end


    # this innerjoins every motif scanned onto the paired peak table from each tf given
    matchingpeaks = innerjoiner(pairedpeaks,motifscandict,tfvector)
    #this plainly compares the midpoint distances of the motifs given for each same peak
    motifdistancerloop!(matchingpeaks, tfvector, centertfs=centertfs)
    # adds on the full sequence for the peak
    identifyseqs!(matchingpeaks,genomefile)

    mymotifs, motifvec, motifscandict, pairedpeaks, matchingpeaks

end

function associatesnps(pairedpeaks, snpnamevec, tfvector, zscorevec, disruptcutoffvec, snpsetvec, motifvec, genomefile)
    # makes a numbered list if no snpset names are given
    if snpnamevec==nothing
        snpnamevec = string.(1:size(snpsetvec,1))
    end

    # 
    snpscandict = repeatscanmachine(snpsetvec,snpnamevec,motifvec,genomefile)

    snpscanvec=DataFrame[]
    for snpname in snpnamevec
        push!(snpscanvec, snpscandict[snpname])
    end

    finalsnpdict= snpdictmaker(pairedpeaks, snpscanvec, snpnamevec, tfvector, zscorevec, disruptcutoffvec)

    finalsnpdict, snpscandict
end



#split this function up a bit more into motif stuff, scanning stuff and crossing stuff
function throughput(genomefile, peakfiles, motiffile, tfvector, snpsetvec, zscorevec, disruptcutoffvec; maxonly=true, snpnamevec=nothing, centertfs=nothing, cutoff=0)
    # flattens the replicate trials into one set of peaks for each tf, adds ID col for peaks in file, puts in peakdict. also makes tfcombo as a side thing
    (peakdict, tfcombo) = loadreplicates(peakfiles, genomefile, tfvector)

    # 
    (mymotifs, motifvec, motifscandict, pairedpeaks) = identifymotifs(motiffile, tfvector, snpsetvec, peakdict, maxonly, cutoff, snpnamevec, centertfs)
    
    finalsnpdict, snpscandict = associatesnps(pairedpeaks, snpnamevec, tfvector, zscorevec, disruptcutoffvec, snpsetvec, motifvec, genomefile)

    # for snpname in snpnamevec
    #     sanitychecker(finalsnpdict[snpname], tfvector)
    #     distancechecker(finalsnpdict[snpname])
    # end
    
    finalsnpdict
end


# this is for making that datatable of the most frequent distances
function distancefreqs(df::DataFrame,key::String)
    output=sort(combine(groupby(df, Symbol("motifdistance",key)), nrow => Symbol("count",key)),Symbol("count",key), rev=true)
    output[!,:ID] = 1:size(output,1)
    output
end

function distancefreqs(dict::Dict{String, DataFrame})
    countdict = Dict{String, DataFrame}()
    outputdf = DataFrame(ID = 1:2000)
    for key in keys(dict)
        countdict[key] = distancefreqs(dict[key], key)
        outputdf = outerjoin(outputdf,countdict[key], on=:ID)
    end
    countdict,outputdf
end

# you want this function to output a dictionary of matchingpeaks tables on for each unique pair, so there should be about 15?
# this is the biggest bottleneck, i think, cuz it just runs the whole idmot thing for eachone of the pairs to make sure the interactions marked are limited to each pair
function matchingpeakpairs(motiffile, tfpairvec, snpsetvec, peakdict, maxonly, cutoff, snpnamevec)
    matchingdict = Dict{String, DataFrame}()
    for tfpair in tfpairvec
        (x,y,z,xx,matchingdict[join(tfpair,"_")]) = identifymotifs(motiffile, tfpair, snpsetvec, peakdict, maxonly, cutoff, snpnamevec, [tfpair[1]])
    end

    matchingdict
end

# here for posterity if we need to look at these again
function originalappendsnpscores!(gvat,dict)
    for key in keys(dict)
        dict[key] = leftjoin(dict[key], @subset(gvat, occursin.(:TF,key))[:,[:ID,:TF,:pbs,:pval, :oligo_auc,:oligo_pval]], on=:gvatID=>:ID)
    end
end

function snpscorer(snpset,dict,pairs; snpt="none",tf="either",disrupt=false,mymotifs=nothing,genomelocation=nothing)
    outputdict = Dict{String, DataFrame}()
    if snpt ==  "none"
        for pair in pairs
            key = join(pair,"_")
            outputdict[key] = @subset(@subset(dict[key], :gvatID .== 0), :gwasID .== 0)
        end
    elseif snpt == "gvat"
        if tf == "either"
            for pair in pairs
                key = join(pair,"_")
                outputdict[key] = innerjoin(dict[key], @subset(snpset, occursin.(:TF,key))[:,[:ID,:TF,:pbs,:pval, :oligo_auc,:oligo_pval]], on=:gvatID=>:ID)
                if disrupt == true
                    outputdict[key] = @subset(outputdict[key],:pval .< 0.05)
                end
            end
        elseif tf == "tf1"
            # i worry that this is including snps on sites that have low binding scores, but we can probably add a filter at some point after if we need it
            # this is only getting the distance distributions for only the lines subset by tf
            # wait, does this mean that every line in the default gvat happens twice, once for the strongest snp for each tf? does that even matter?
            for pair in pairs
                key = join(pair,"_")
                outputdict[key] = @subset(innerjoin(dict[key], @subset(snpset, occursin.(:TF,key))[:,[:ID,:TF,:pbs,:pval, :oligo_auc,:oligo_pval]], on=:gvatID=>:ID), :TF .== pair[1])
                if disrupt == true
                    outputdict[key] = @subset(outputdict[key],:pval .< 0.05)
                end
            end
        elseif tf == "tf2"
            for pair in pairs
                key = join(pair,"_")
                outputdict[key] = @subset(innerjoin(dict[key], @subset(snpset, occursin.(:TF,key))[:,[:ID,:TF,:pbs,:pval,:oligo_auc,:oligo_pval]], on=:gvatID=>:ID), :TF .== pair[2])
                if disrupt == true
                    outputdict[key] = @subset(outputdict[key],:pval .< 0.05)
                end
            end
        elseif tf == "both"
            for pair in pairs
                key = join(pair,"_")
                outputdict[key] = leftjoin(dict[key], @subset(snpset, occursin.(:TF,pair[1]))[:,[:ID,:TF,:pbs,:pval, :oligo_auc,:oligo_pval]], on=:gvatID=>:ID)
                outputdict[key] = innerjoin(outputdict[key], @subset(snpset, occursin.(:TF,pair[2]))[:,[:ID,:TF,:pbs,:pval, :oligo_auc,:oligo_pval]], on=:gvatID=>:ID, makeunique=true)
            end
        end      
    elseif snpt == "gwas"
        gwasrep = scanmachine(snpset,mymotifs,genomelocation)
        # adding an extra subset to only include prrefalts that aren't 0
        if tf == "either"
            for pair in pairs
                key = join(pair,"_")
                outputdict[key] = innerjoin(dict[key], @subset(gwasrep, occursin.(:MotifName,key))[:,[:ID,:snp,:MotifName,:LR_RefAlt,:PR_RefAlt]], on=:gwasID=>:ID)
                if disrupt == true
                    outputdict[key] = @subset(outputdict[key],:PR_RefAlt .!= 0.0)
                end
            end
        elseif tf == "tf1"
            for pair in pairs
                key = join(pair,"_")
                outputdict[key] = innerjoin(dict[key], @subset(gwasrep, occursin.(:MotifName,pair[1]))[:,[:ID,:snp,:MotifName,:LR_RefAlt,:PR_RefAlt]], on=:gwasID=>:ID)
                if disrupt == true
                    outputdict[key] = @subset(outputdict[key],:PR_RefAlt .!= 0.0)
                end            
            end
        elseif tf == "tf2"
            for pair in pairs
                key = join(pair,"_")
                outputdict[key] = innerjoin(dict[key], @subset(gwasrep, occursin.(:MotifName,pair[2]))[:,[:ID,:snp,:MotifName,:LR_RefAlt,:PR_RefAlt]], on=:gwasID=>:ID)
                if disrupt == true
                    outputdict[key] = @subset(outputdict[key],:PR_RefAlt .!= 0.0)
                end
            end
        elseif tf == "both"
            for pair in pairs
                key = join(pair,"_")
                outputdict[key] = leftjoin(dict[key], @subset(gwasrep, occursin.(:MotifName,pair[1]))[:,[:ID,:snp,:MotifName,:LR_RefAlt,:PR_RefAlt]], on=:gwasID=>:ID)
                outputdict[key] = innerjoin(outputdict[key], @subset(gwasrep, occursin.(:MotifName,pair[2]))[:,[:ID,:snp,:MotifName,:LR_RefAlt,:PR_RefAlt]], on=:gwasID=>:ID, makeunique=true)
                if disrupt == true
                    outputdict[key] = @subset(outputdict[key],:PR_RefAlt .!= 0.0)
                end
            end
        end  
    end
    outputdict
end

# this corrects all the motif distance stuff to flip what negative and positive numbers represent and how they interpret strands
function distancecorrection(outputdf::DataFrame, pair; rand=nothing)
  if rand === nothing
    # this is just to make +10 mean tf2 is ahead on the strand; this is not the cause for the divide, they look different in graphs cuz of binnig i suppose
    outputdf[!,Symbol("motifdistance",pair[1],"_",pair[2])] = outputdf[!,Symbol("motifdistance",pair[1],"_",pair[2])].*-1
    # this is the method he suggested for the strand flipping
    outputdf[!,Symbol("motifdistance",pair[1],"_",pair[2])] = ifelse.(outputdf[!,Symbol(pair[1],"strand")] .== "-", -1, 1) .* outputdf[!,Symbol("motifdistance",pair[1],"_",pair[2])]
  else
    outputdf[!,Symbol("motifdistance",pair[1],"_",pair[2],"_rand",rand)] = outputdf[!,Symbol("motifdistance",pair[1],"_",pair[2],"_rand",rand)].*-1
    outputdf[!,Symbol("motifdistance",pair[1],"_",pair[2],"_rand",rand)] = ifelse.(outputdf[!,Symbol(pair[1],"strand")] .== "-", -1, 1) .* outputdf[!,Symbol("motifdistance",pair[1],"_",pair[2],"_rand",rand)]
  end
  outputdf
end

function distancecorrection(dict::Dict{String, DataFrame}, pairs::Vector{Vector{String}}; randrange=nothing)
    outputdict = deepcopy(dict)
    for pair in pairs
        name = join(pair,"_")
        outputdict[name] = distancecorrection(outputdict[name], pair)
        if randrange !== nothing
            for i in 1:randrange
                distancecorrection(outputdict[name], pair, rand=i)
            end
        end
    end    
outputdict
end

#just to not clutter randomizerer
function shortr!(df,pair,i,loops)
    offlength = first(unique(df[!,Symbol(pair[i],"motifstop")] .- df[!,Symbol(pair[i],"motifstart")] .+ 1))
    df[!,Symbol(pair[i],string("motifstartrand",loops))] = [rand(row[Symbol(pair[i],"start")]:(row[Symbol(pair[i],"stop")] - offlength + 1)) for row in eachrow(df)]
    df[!,Symbol(pair[i],string("motifstoprand",loops))] = df[!,Symbol(pair[i],string("motifstartrand",loops))] .+ offlength .- 1;
end

# this returns the same tables with the added randomized start, stop and distance columns
function randomizerer!(matchingpeaks::DataFrame, pair::Vector{String}; loops="", tf=nothing)
    if tf === nothing
        range = 1:2
        tf1s = string("rand",loops)
        tf2s = string("rand",loops)
    elseif tf == 1
        range = 1:1
        tf1s = string("rand",loops)
        tf2s = ""
    elseif tf == 2
        range = 2:2
        tf1s = ""
        tf2s = string("rand",loops)
    end

    for i in range
        shortr!(matchingpeaks,pair,i,loops)
    end
    matchingpeaks[!,Symbol("motifdistance",pair[1],"_",pair[2],"_rand",loops)] = div.(matchingpeaks[!,Symbol(pair[1],"motifstart",tf1s)] .+ matchingpeaks[!,Symbol(pair[1],"motifstop",tf1s)], 2) .- div.(matchingpeaks[!,Symbol(pair[2],"motifstart",tf2s)] .+ matchingpeaks[!,Symbol(pair[2],"motifstop",tf2s)], 2);
end

# branch for using the function with lists of dataframes
function randomizerer!(matchingdict::Dict{String, DataFrame}, pairs::Vector{Vector{String}}; loops=50, tf=nothing)
    for pair in pairs
        name = join(pair,"_")
        for loop in 1:loops
            randomizerer!(matchingdict[name], pair, loops=string(loop), tf=tf)
        end
    end
end

# just a lil something to increase the instances of snps to see them better graphed
function instanceincreaser(df,n)
    loaddf = copy(df)
        for n in 1:n
            append!(loaddf,copy(df))
        end
        loaddf
end

function xticker!(output)
    old_xticks = xticks(output[1]) # grab xticks of the 1st subplot
    new_xticks = ([147, -147, 150, -150],["147", "-147"])
    vline!(new_xticks[1][1:2])
    keep_indices = findall(x -> all(x .≠ new_xticks[1]), old_xticks[1])
    merged_xticks = (old_xticks[1][keep_indices] ∪ new_xticks[1], old_xticks[2][keep_indices] ∪ new_xticks[2])
    xticks!(merged_xticks)
end


# outputs the wiry plots of normal vs rand
function referencegraph(matchingpeaks::DataFrame, pair::Vector{String}; bins=((-1000:10:1000) .+ 0.5), xlims=(-1000, 1000),snpinstances=5,absolute=false)
    # bins = (-1000:10:1000) .+ 0.5
    plot()
    name = join(pair,"_")
    if absolute==false
        stephist!((matchingpeaks[!,Symbol("motifdistance",pair[1],"_",pair[2])]), title=string(name," distances ", "(n=",size(matchingpeaks,1),")"), bins=bins, labels="standard", normalize=true)
        plot!(fmt=:png)
        stephist!(stackdataframe(select(matchingpeaks[!,:],r"_rand")), title=string(name," distances ", "(n=",size(matchingpeaks,1),")"), bins=bins, labels="randomized", normalize=true)
        plot!(fmt=:png)
    else
        stephist!(abs.(matchingpeaks[!,Symbol("motifdistance",pair[1],"_",pair[2])]), title=string(name," distances ", "(n=",size(matchingpeaks,1),")"), bins=bins, labels="standard", normalize=true)
        plot!(fmt=:png)
        stephist!(abs.(stackdataframe(select(matchingpeaks[!,:],r"_rand"))), title=string(name," distances ", "(n=",size(matchingpeaks,1),")"), bins=bins, labels="randomized", normalize=true)
        plot!(fmt=:png)
    end
    # this down here getting the distance frequency of sites with ascertained binding snps and sites with significantly preferentially binding snps
    # undo this later, just focusing on the distances rn
    #histogram!((instanceincreaser(@subset(@subset(matchingpeaks,:gvatID .> 0),:oligo_pval .< 0.05),snpinstances)[!,Symbol("motifdistance",pair[1],"_",pair[2])]), title=string(name," distances ", "(n=",size(matchingpeaks,1),")"), bins=bins, labels="with snps x5 incidence")
    #plot!(fmt=:png)
    #histogram!((instanceincreaser(@subset(@subset(@subset(matchingpeaks,:gvatID .> 0),:pval .< 0.05),:oligo_pval .< 0.05),snpinstances)[!,Symbol("motifdistance",pair[1],"_",pair[2])]), title=string(name," distances ", "(n=",size(matchingpeaks,1),")"), bins=bins, labels=string("with pbsnps x",snpinstances," incidence"), normalize=true)
    #plot!(fmt=:png)

    # this should produces marks at +/-147 just to see if the nucleosome length is relevant
    output = plot!(size=(1200, 300), xlims=xlims, xticks=25)
    #xticker!(output)
    output
end

# sticks those outputs in a dictionary so we have each pair's output on hand
function referencegraph(matchingdict::Dict{String, DataFrame}, pairs::Vector{Vector{String}}; bins=((-1000:10:1000) .+ 0.5), xlims=(-1000, 1000),snpinstances=5,absolute=false)
    plotdict = Dict{String, Plots.Plot{Plots.GRBackend}}()
    for pair in pairs
        name = join(pair,"_")
        plotdict[name] = referencegraph(matchingdict[name], pair, bins=bins, xlims=xlims,snpinstances=snpinstances,absolute=absolute)
    end
    plotdict
end

function significancegraph(matchingpeaks::DataFrame, pair::Vector{String}; bins=(-1500:10:1500), xlims=(-1500, 1500),snpinstances=1,absolute=false,statsdict=statsdict)
    name = join(pair,"_")
    if absolute==false
        truecount = fit(Histogram, (matchingpeaks[!,Symbol("motifdistance",pair[1],"_",pair[2])]), bins).weights
        randcount = fit(Histogram, (stackdataframe(select(matchingpeaks[!,:],r"_rand"))), bins).weights
    else
        truecount = fit(Histogram, abs.(matchingpeaks[!,Symbol("motifdistance",pair[1],"_",pair[2])]), bins).weights
        randcount = fit(Histogram, abs.(stackdataframe(select(matchingpeaks[!,:],r"_rand"))), bins).weights
    end
    randcount = sum(truecount) .* randcount ./ sum(randcount)
    #truecount = fit(Histogram, (matchingpeaks[!,Symbol("motifdistance",pair[1],"_",pair[2])]),weights((matchingpeaks[!,Symbol("motifdistance",pair[1],"_",pair[2])])), bins)
    #randcount = fit(Histogram, (stackdataframe(select(matchingpeaks[!,:],r"_rand"))),weights((stackdataframe(select(matchingpeaks[!,:],r"_rand")))), bins)
    stats = DataFrame(Dist=bins[1:end-1], Count=truecount, RandCount=randcount)
    # return stats
    stats.FC = stats.Count./stats.RandCount
    stats.PoissonPvalue = @with stats ccdf.(Poisson.(:RandCount), :Count)
    statsdict[name] = stats
    bar(stats.Dist, -log10.(stats.PoissonPvalue),title=string(name," distances ", "(n=",size(matchingpeaks,1),")"),bins=bins,size=(1200, 300), xlims=xlims, ylims=(0,20),labels="-log10 p-value")
    plot!(fmt=:png)

    # copying all the stuff from the previous graphing function
    # this down here getting the distance frequency of sites with ascertained binding snps and sites with significantly preferentially binding snps
    #histogram!((instanceincreaser(@subset(@subset(matchingpeaks,:gvatID .> 0),:oligo_pval .< 0.05),snpinstances)[!,Symbol("motifdistance",pair[1],"_",pair[2])]), title=string(name," distances ", "(n=",size(matchingpeaks,1),")"),fillalpha=0.4, bins=bins, labels="with snps x5 incidence")
    #plot!(fmt=:png)
    # remember that the oligo_pval filter means its likely the tf will bind to this oligo and the pval filter means its likely there is a significant difference in the binding affinities of the ref and alt oligos, meaning its a disrupting snp?
    # that is to say, if one of these lands on a specific distance frequency peak, it means cobound peaks where the motifs are this distance apart are more likely to encounter/be affected by a pbsnp(?)
    # if absolute==false
    #     histogram!((instanceincreaser(@subset(@subset(@subset(matchingpeaks,:gvatID .> 0),:pval .< 0.05),:oligo_pval .< 0.05),snpinstances)[!,Symbol("motifdistance",pair[1],"_",pair[2])]), title=string(name," distances ", "(n=",size(matchingpeaks,1),")"),fillalpha=0.9, bins=bins, labels=string("with pbsnps x",snpinstances," incidence"))
    #     plot!(fmt=:png)
    # else
    #     histogram!(abs.(instanceincreaser(@subset(@subset(@subset(matchingpeaks,:gvatID .> 0),:pval .< 0.05),:oligo_pval .< 0.05),snpinstances)[!,Symbol("motifdistance",pair[1],"_",pair[2])]), title=string(name," distances ", "(n=",size(matchingpeaks,1),")"),fillalpha=0.9, bins=bins, labels=string("with pbsnps x",snpinstances," incidence"))
    #     plot!(fmt=:png)
    # end
    # this should produces marks at +/-147 and +/-200, just to see if the nucleosome length is relevant
    output = plot!(size=(1200, 300), xlims=xlims, xticks=25)

    #xticker!(output)
    output
end

function significancegraph(matchingdict::Dict{String, DataFrame}, pairs::Vector{Vector{String}}; bins=(-1500:10:1500), xlims=(-1500, 1500),snpinstances=1,absolute=false)
    plotdict = Dict{String, Plots.Plot{Plots.GRBackend}}()
    statsdict = Dict{String, DataFrame}()
    for pair in pairs
        name = join(pair,"_")
        plotdict[name] = significancegraph(matchingdict[name], pair, bins=bins, xlims=xlims,snpinstances=snpinstances,absolute=absolute,statsdict=statsdict)
    end
    plotdict,statsdict
end

# this subsets for different strands
# oppositestrand is just there if u want to pull out the places where they're on opposite ones
function strandsplit(df::DataFrame,pair,direction; oppositestrand=false)
    if oppositestrand == false
        outputdf = bettersubset(bettersubset(df,Symbol(pair[1],"strand"), isequal(string(direction))),Symbol(pair[2],"strand"),isequal(string(direction)))
    else
        outputdf = bettersubset(bettersubset(df,Symbol(pair[1],"strand"), isequal(string(direction))),Symbol(pair[2],"strand"),!isequal(string(direction)))
    end
    outputdf
end

function strandsplit(dict::Dict{String, DataFrame},pairs,direction; oppositestrand=false)
    outputdict = Dict{String, DataFrame}()
    for pair in pairs
        name = join(pair,"_")
        outputdict[name] = strandsplit(dict[name],pair,direction,oppositestrand=oppositestrand)
    end
    outputdict
end

# these are copies of the previous functions altered to accomodate the strand direction pairs-----------------------------------------------

# outputs the wiry plots of normal vs rand
function pairgraph(matchingpeaks1::DataFrame, matchingpeaks2::DataFrame, pair::Vector{String}; bins=((-1000:10:1000) .+ 0.5), xlims=(-1000, 1000),snpinstances=5,absolute=false)
    # bins = (-1000:10:1000) .+ 0.5
    plot()
    name = join(pair,"_")
    if absolute==false
        stephist!((matchingpeaks1[!,Symbol("motifdistance",pair[1],"_",pair[2])]), title=string(name," distances ", "(n=",size(matchingpeaks1,1),")"), bins=bins, labels="plus", normalize=true)
        plot!(fmt=:png)
        stephist!(stackdataframe(select(matchingpeaks1[!,:],r"_rand")), title=string(name," distances ", "(n=",size(matchingpeaks1,1),")"), bins=bins, labels="plusrand", normalize=true)
        plot!(fmt=:png)

        stephist!((matchingpeaks2[!,Symbol("motifdistance",pair[1],"_",pair[2])]), title=string(name," distances ", "(n=",size(matchingpeaks2,1),")"), bins=bins, labels="minus", normalize=true)
        plot!(fmt=:png)
        stephist!(stackdataframe(select(matchingpeaks2[!,:],r"_rand")), title=string(name," distances ", "(n=",size(matchingpeaks2,1),")"), bins=bins, labels="minusrand", normalize=true)
        plot!(fmt=:png)
    else
        stephist!(abs.(matchingpeaks1[!,Symbol("motifdistance",pair[1],"_",pair[2])]), title=string(name," distances ", "(n=",size(matchingpeaks1,1),")"), bins=bins, labels="plus", normalize=true)
        plot!(fmt=:png)
        stephist!(abs.(stackdataframe(select(matchingpeaks1[!,:],r"_rand"))), title=string(name," distances ", "(n=",size(matchingpeaks1,1),")"), bins=bins, labels="plusrand", normalize=true)
        plot!(fmt=:png)
        stephist!(abs.(matchingpeaks2[!,Symbol("motifdistance",pair[1],"_",pair[2])]), title=string(name," distances ", "(n=",size(matchingpeaks2,1),")"), bins=bins, labels="minus", normalize=true)
        plot!(fmt=:png)
        stephist!(abs.(stackdataframe(select(matchingpeaks2[!,:],r"_rand"))), title=string(name," distances ", "(n=",size(matchingpeaks2,1),")"), bins=bins, labels="minusrand", normalize=true)
        plot!(fmt=:png)
    end
    
    
    
    # this down here getting the distance frequency of sites with ascertained binding snps and sites with significantly preferentially binding snps
    # undo this later, just focusing on the distances rn
    #histogram!((instanceincreaser(@subset(@subset(matchingpeaks,:gvatID .> 0),:oligo_pval .< 0.05),snpinstances)[!,Symbol("motifdistance",pair[1],"_",pair[2])]), title=string(name," distances ", "(n=",size(matchingpeaks,1),")"), bins=bins, labels="with snps x5 incidence")
    #plot!(fmt=:png)
    #histogram!((instanceincreaser(@subset(@subset(@subset(matchingpeaks,:gvatID .> 0),:pval .< 0.05),:oligo_pval .< 0.05),snpinstances)[!,Symbol("motifdistance",pair[1],"_",pair[2])]), title=string(name," distances ", "(n=",size(matchingpeaks,1),")"), bins=bins, labels="with pbsnps x5 incidence")
    #plot!(fmt=:png)

    # this should produces marks at +/-147 just to see if the nucleosome length is relevant
    output = plot!(size=(1200, 300), xlims=xlims, xticks=25)
    xticker!(output)
    output
end

# sticks those outputs in a dictionary so we have each pair's output on hand
function pairgraph(matchingdict1::Dict{String, DataFrame}, matchingdict2::Dict{String, DataFrame}, pairs::Vector{Vector{String}}; bins=((-1000:10:1000) .+ 0.5), xlims=(-1000, 1000),snpinstances=5,absolute=false)
    plotdict = Dict{String, Plots.Plot{Plots.GRBackend}}()
    for pair in pairs
        name = join(pair,"_")
        plotdict[name] = pairgraph(matchingdict1[name], matchingdict2[name], pair, bins=bins, xlims=xlims,snpinstances=snpinstances,absolute=absolute)
    end
    plotdict
end

function pairsignificancegraph(matchingpeaks1::DataFrame, matchingpeaks2::DataFrame, pair::Vector{String}; bins=(-1500:10:1500), xlims=(-1500, 1500),snpinstances=1,absolute=false,l1="same",l2="opposite",stats1dict=stats1dict,stats2dict=stats2dict)
    name = join(pair,"_")
    if absolute==false
        truecount1 = fit(Histogram, (matchingpeaks1[!,Symbol("motifdistance",pair[1],"_",pair[2])]), bins).weights
        randcount1 = fit(Histogram, (stackdataframe(select(matchingpeaks1[!,:],r"_rand"))), bins).weights
        truecount2 = fit(Histogram, (matchingpeaks2[!,Symbol("motifdistance",pair[1],"_",pair[2])]), bins).weights
        randcount2 = fit(Histogram, (stackdataframe(select(matchingpeaks2[!,:],r"_rand"))), bins).weights
    else
        truecount1 = fit(Histogram, abs.(matchingpeaks1[!,Symbol("motifdistance",pair[1],"_",pair[2])]), bins).weights
        randcount1 = fit(Histogram, abs.(stackdataframe(select(matchingpeaks1[!,:],r"_rand"))), bins).weights
        truecount2 = fit(Histogram, abs.(matchingpeaks2[!,Symbol("motifdistance",pair[1],"_",pair[2])]), bins).weights
        randcount2 = fit(Histogram, abs.(stackdataframe(select(matchingpeaks2[!,:],r"_rand"))), bins).weights
    end
    randcount1 = sum(truecount1) .* randcount1 ./ sum(randcount1)
    randcount2 = sum(truecount2) .* randcount2 ./ sum(randcount2)
    
    

    # keep in mind, these are relative graphs, so they're built out of their own randomization, which is not the same for both
    stats1 = DataFrame(Dist=bins[1:end-1], Count=truecount1, RandCount=randcount1)
    stats1.FC = stats1.Count./stats1.RandCount
    stats1.PoissonPvalue = @with stats1 ccdf.(Poisson.(:RandCount), :Count)
    bar(stats1.Dist, -log10.(stats1.PoissonPvalue),title=string(name," distances ", "(n=",size(matchingpeaks,1),")"),bins=bins,size=(1200, 300), xlims=xlims, ylims=(0,15), alpha=0.6, labels=l1)
    plot!(fmt=:png)
    stats2 = DataFrame(Dist=bins[1:end-1], Count=truecount2, RandCount=randcount2)
    stats2.FC = stats2.Count./stats2.RandCount
    stats2.PoissonPvalue = @with stats2 ccdf.(Poisson.(:RandCount), :Count)
    bar!(stats2.Dist, -log10.(stats2.PoissonPvalue),title=string(name," distances ", "(n=",size(matchingpeaks,1),")"),bins=bins,size=(1200, 300), xlims=xlims, ylims=(0,15), alpha=0.6, labels=l2)
    plot!(fmt=:png)
    stats1dict[name] = stats1
    stats2dict[name] = stats2

    # copying all the stuff from the previous graphing function
    # this down here getting the distance frequency of sites with ascertained binding snps and sites with significantly preferentially binding snps
    #histogram!((instanceincreaser(@subset(@subset(matchingpeaks,:gvatID .> 0),:oligo_pval .< 0.05),snpinstances)[!,Symbol("motifdistance",pair[1],"_",pair[2])]), title=string(name," distances ", "(n=",size(matchingpeaks,1),")"),fillalpha=0.4, bins=bins, labels="with snps x5 incidence")
    #plot!(fmt=:png)
    # remember that the oligo_pval filter means its likely the tf will bind to this oligo and the pval filter means its likely there is a significant difference in the binding affinities of the ref and alt oligos, meaning its a disrupting snp?
    # that is to say, if one of these lands on a specific distance frequency peak, it means cobound peaks where the motifs are this distance apart are more likely to encounter/be affected by a pbsnp(?)
    #histogram!((instanceincreaser(@subset(@subset(@subset(matchingpeaks,:gvatID .> 0),:pval .< 0.05),:oligo_pval .< 0.05),snpinstances)[!,Symbol("motifdistance",pair[1],"_",pair[2])]), title=string(name," distances ", "(n=",size(matchingpeaks,1),")"),fillalpha=0.4, bins=bins, labels="with pbsnps x5 incidence")
    #plot!(fmt=:png)

    # this should produces marks at +/-147 and +/-200, just to see if the nucleosome length is relevant
    output = plot!(size=(1200, 300), xlims=xlims, xticks=25)

    xticker!(output)
    output
end

function pairsignificancegraph(matchingdict1::Dict{String, DataFrame}, matchingdict2::Dict{String, DataFrame}, pairs::Vector{Vector{String}}; bins=(-1500:10:1500), xlims=(-1500, 1500),snpinstances=1,absolute=false,l1="same",l2="opposite")
    plotdict = Dict{String, Plots.Plot{Plots.GRBackend}}()
    stats1dict = Dict{String, DataFrame}()
    stats2dict = Dict{String, DataFrame}()
    for pair in pairs
        name = join(pair,"_")
        plotdict[name] = pairsignificancegraph(matchingdict1[name], matchingdict2[name], pair, bins=bins, xlims=xlims,snpinstances=snpinstances,absolute=absolute,l1=l1,l2=l2,stats1dict=stats1dict,stats2dict=stats2dict)
    end
    plotdict,stats1dict,stats2dict
end

# they end here ------------------------------------------------------------------------------------------------------------------

function directionalitybar(df::DataFrame, pair)
    name = join(pair,"_")
    plusmatch = bettersubset(bettersubset(df,Symbol(pair[1],"strand"), isequal("+")),Symbol(pair[2],"strand"),isequal("+"))
    minusmatch = bettersubset(bettersubset(df,Symbol(pair[1],"strand"), isequal("-")),Symbol(pair[2],"strand"),isequal("-"))
    pluscrossmatch = bettersubset(bettersubset(df,Symbol(pair[1],"strand"), isequal("+")),Symbol(pair[2],"strand"),isequal("-"))
    minuscrossmatch = bettersubset(bettersubset(df,Symbol(pair[1],"strand"), isequal("-")),Symbol(pair[2],"strand"),isequal("+"))
    samedir = append!(copy(plusmatch),minusmatch)
    crossdir = append!(copy(pluscrossmatch),minuscrossmatch)
    bar(["match +", "match -", "cross +", "cross -", "match", "cross"],[size(plusmatch,1),size(minusmatch,1),size(pluscrossmatch,1),size(minuscrossmatch,1),size(samedir,1),size(crossdir,1)],title=string(name)), samedir, crossdir
end

function directionalitybar(dict::Dict{String, DataFrame}, pairs)
    plotdict = Dict{String, Plots.Plot{Plots.GRBackend}}()
    samedirdict = Dict{String, DataFrame}()
    crossdirdict = Dict{String, DataFrame}()
    for pair in pairs
        name = join(pair,"_")
        plotdict[name], samedirdict[name], crossdirdict[name] = directionalitybar(dict[name],pair)
    end
    plotdict, samedirdict, crossdirdict
end


function flatsnpgraph(matchingpeaks::DataFrame, pair::Vector{String}; bins=((-1000:10:1000) .+ 0.5), xlims=(-1000, 1000),absolute=false,snpset="either",nstate=false)

    plot()
    name = join(pair,"_")

    if snpset=="either"
        ydf = matchingpeaks
    elseif snpset=="gvat"
        ydf = append!(@subset(@subset(matchingpeaks, :gvatID .== 0), :gwasID .== 0), @subset(matchingpeaks, :gvatID .!= 0))
    elseif snpset=="gwas"
        ydf = append!(@subset(@subset(matchingpeaks, :gvatID .== 0), :gwasID .== 0), @subset(matchingpeaks, :gwasID .!= 0))
    end

    if absolute==false
        stephist!((ydf[!,Symbol("motifdistance",pair[1],"_",pair[2])]), title=string(name," distances ", "(n=",size(matchingpeaks,1),")"), bins=bins, labels="snps", normalize=nstate)
        plot!(fmt=:png)
        stephist!((@subset(@subset(matchingpeaks, :gvatID .== 0), :gwasID .== 0)[!,Symbol("motifdistance",pair[1],"_",pair[2])]), title=string(name," distances ", "(n=",size(matchingpeaks,1),")"), bins=bins, labels="nosnps", normalize=nstate)
        plot!(fmt=:png)
    else
        stephist!(abs.(ydf[!,Symbol("motifdistance",pair[1],"_",pair[2])]), title=string(name," distances ", "(n=",size(matchingpeaks,1),")"), bins=bins, labels="snps", normalize=nstate)
        plot!(fmt=:png)
        stephist!(abs.(@subset(@subset(matchingpeaks, :gvatID .== 0), :gwasID .== 0)[!,Symbol("motifdistance",pair[1],"_",pair[2])]), title=string(name," distances ", "(n=",size(matchingpeaks,1),")"), bins=bins, labels="nosnps", normalize=nstate)
        plot!(fmt=:png)
    end

    # this should produces marks at +/-147 just to see if the nucleosome length is relevant
    output = plot!(size=(1200, 300), xlims=xlims, xticks=25)
    xticker!(output)
    output
end


function flatsnpgraph(matchingdict::Dict{String, DataFrame}, pairs::Vector{Vector{String}}; bins=((-1000:10:1000) .+ 0.5), xlims=(-1000, 1000),absolute=false,snpset="either",nstate=false)
    plotdict = Dict{String, Plots.Plot{Plots.GRBackend}}()
    for pair in pairs
        name = join(pair,"_")
        plotdict[name] = flatsnpgraph(matchingdict[name], pair, bins=bins, xlims=xlims,absolute=absolute,snpset=snpset,nstate=nstate)
    end
    plotdict
end


function relativegrapher(matchingpeaks::DataFrame, pair::Vector{String}; bins=(-1500:10:1500), xlims=(-1500, 1500),absolute=false,snpset="either")
    name = join(pair,"_")

    if snpset=="either"
        ydf = matchingpeaks
    elseif snpset=="gvat"
        ydf = append!(@subset(@subset(matchingpeaks, :gvatID .== 0), :gwasID .== 0), @subset(matchingpeaks, :gvatID .!= 0))
    elseif snpset=="gwas"
        ydf = append!(@subset(@subset(matchingpeaks, :gvatID .== 0), :gwasID .== 0), @subset(matchingpeaks, :gwasID .!= 0))
    end

    if absolute==false
        ysnps = fit(Histogram, (ydf[!,Symbol("motifdistance",pair[1],"_",pair[2])]), bins)
        nsnps = fit(Histogram, (@subset(@subset(matchingpeaks, :gvatID .== 0), :gwasID .== 0)[!,Symbol("motifdistance",pair[1],"_",pair[2])]), bins)
    else
        ysnps = fit(Histogram, abs.(ydf[!,Symbol("motifdistance",pair[1],"_",pair[2])]), bins).weights
        nsnps = fit(Histogram, abs.(@subset(@subset(matchingpeaks, :gvatID .== 0), :gwasID .== 0)[!,Symbol("motifdistance",pair[1],"_",pair[2])]), bins).weights
    end


    # might be overkill, trying something simpler
    # if absolute==false
    #     ysnps = fit(Histogram, (ydf[!,Symbol("motifdistance",pair[1],"_",pair[2])]), bins).weights
    #     nsnps = fit(Histogram, (@subset(@subset(matchingpeaks, :gvatID .== 0), :gwasID .== 0)[!,Symbol("motifdistance",pair[1],"_",pair[2])]), bins).weights
    # else
    #     ysnps = fit(Histogram, abs.(ydf[!,Symbol("motifdistance",pair[1],"_",pair[2])]), bins).weights
    #     nsnps = fit(Histogram, abs.(@subset(@subset(matchingpeaks, :gvatID .== 0), :gwasID .== 0)[!,Symbol("motifdistance",pair[1],"_",pair[2])]), bins).weights
    # end
    # nsnps = sum(ysnps) .* nsnps ./ sum(nsnps) 
    # stats = DataFrame(Dist=bins[1:end-1], Count=ysnps, RandCount=nsnps)
    # # return stats
    # stats.FC = stats.Count./stats.RandCount
    # stats.PoissonPvalue = @with stats ccdf.(Poisson.(:RandCount), :Count)
    # bar(stats.Dist, -log10.(stats.PoissonPvalue),title=string(name," distances ", "(n=",size(matchingpeaks,1),")"),bins=bins,size=(1200, 300), xlims=xlims, ylims=(0,10))
    # plot!(fmt=:png)

    # this should produces marks at +/-147 and +/-200, just to see if the nucleosome length is relevant
    output = plot!(size=(1200, 300), xlims=xlims, xticks=25)

    xticker!(output)
    output
end

function relativegrapher(matchingdict::Dict{String, DataFrame}, pairs::Vector{Vector{String}}; bins=(-1500:10:1500), xlims=(-1500, 1500),absolute=false,snpset="either")
    plotdict = Dict{String, Plots.Plot{Plots.GRBackend}}()
    for pair in pairs
        name = join(pair,"_")
        plotdict[name] = relativegrapher(matchingdict[name], pair, bins=bins, xlims=xlims,absolute=absolute,snpset=snpset)
    end
    plotdict
end

# just to speed up looking
function distancerangesubset(dict,pair,range,lim; absolute=false)
    key = join(pair,"_")
    if absolute == false
        df = dict[key]
    else
        df = deepcopy(dict[key])
        df[:,Symbol("motifdistance",key)] = abs.(df[:,Symbol("motifdistance",key)])
    end
    @chain df begin
        bettersubset(_, Symbol("motifdistance",key), rangesubset(range[1],range[2]))
        bettersubset(_,Symbol(pair[1],"prmax"),ismorer(lim))
        bettersubset(_,Symbol(pair[2],"prmax"),ismorer(lim))
        sort(_,:FC, rev=true)
    end
    # sort(bettersubset(dict[key], Symbol("motifdistance",key), rangesubset(range[1],range[2])),:score, rev=true)
end

function peaksub!(dict,col,lim;check="more")
    if check == "more"
        for key in keys(dict)
            dict[key] = bettersubset(dict[key],col, ismorer(lim))
        end
    else
        for key in keys(dict)
            dict[key] = bettersubset(dict[key],col, !ismorer(lim))
        end
    end
end

function addpairFCs!(dict,peakdict,pairvec)
    for pair in pairvec
        key = join(pair,"_")
        t1 = peakdict[pair[1]][:,[:ind,:FC]]
        rename!(t1, :ind => :peakind)
        rename!(t1, string.(pair[1], names(t1)))

        t2 = peakdict[pair[2]][:,[:ind,:FC]]
        rename!(t2, :ind => :peakind)
        rename!(t2, string.(pair[2], names(t2)))

        dict[key] = innerjoin(dict[key], t1, on = Symbol(pair[1],"peakind"))
        dict[key] = innerjoin(dict[key], t2, on = Symbol(pair[2],"peakind"))

    end
end

function peaksubpair!(dict,pairvec,col,lim1,lim2;check="more")
    if check == "more"
        for pair in pairvec
            key = join(pair,"_")
            #passcol1=string()
            #passcol2=string()
            dict[key] = bettersubset(dict[key],string(pair[1],col), ismorer(lim1))
            dict[key] = bettersubset(dict[key],string(pair[2],col), ismorer(lim2))
        end
    else
        for pair in pairvec
            key = join(pair,"_")
            dict[key] = bettersubset(dict[key],col, !ismorer(lim1))
            dict[key] = bettersubset(dict[key],col, !ismorer(lim2))
        end
    end
end

function motifminprmax!(dict::Dict,tag::String,lim::Number)
    tag = Regex(tag)
    for key in keys(dict)
        colnamevec = names(dict[key])[occursin.(tag,names(dict[key]))]
        for col in colnamevec
            dict[key] = bettersubset(dict[key],string(col), ismorer(lim))
        end
    end    
end

function matricize(statsdict, pairvec;stats2dict=nothing,removepairvec=["PDX1/NKX6-1","NKX6-1/PDX1"])
    outputmatrix=Matrix
    outputdf=DataFrame
    if stats2dict === nothing
        tempdf = DataFrame(Pair=String[], Dist=Int64[], PoissonPvalue=Float64[])
        for pair in pairvec
            key = join(pair,"_")
            # first we clean up the statsdict tables
            statsdict[key][!,:Pair] .= string(pair[1], "/", pair[2])
            statsdict[key][!,:PoissonPvalue] = -log10.(statsdict[key][!,:PoissonPvalue])
            append!(tempdf,statsdict[key][:,[:Pair,:Dist,:PoissonPvalue]])
        end
        for removepair in removepairvec
            tempdf = @subset(tempdf,.!occursin.(removepair,:Pair))
        end
        outputdf = unstack(tempdf,:Pair,:Dist,:PoissonPvalue)
        outputmatrix = Matrix(unstack(tempdf,:Pair,:Dist,:PoissonPvalue)[!,2:end])
        outputdf,outputmatrix
    else
        tempdf = DataFrame(Pair=String[], Relation=String[], Dist=Int64[], PoissonPvalue=Float64[])
        for pair in pairvec
            key = join(pair,"_")
            # first we clean up the statsdict tables
            statsdict[key][!,:Pair] .= string(pair[1], "/", pair[2], " Cis")
            stats2dict[key][!,:Pair] .= string(pair[1], "/", pair[2], " Trans")
            statsdict[key][!,:Relation] .= string("Cis")
            stats2dict[key][!,:Relation] .= string("Trans")
            statsdict[key][!,:PoissonPvalue] = -log10.(statsdict[key][!,:PoissonPvalue])
            stats2dict[key][!,:PoissonPvalue] = -log10.(stats2dict[key][!,:PoissonPvalue])
            #statsdict[key] = innerjoin(statsdict[key],stats2dict[key], on=[:Pair,:Dist],makeunique=true)
            
            append!(tempdf,statsdict[key][:,[:Pair,:Relation,:Dist,:PoissonPvalue]])
            append!(tempdf,stats2dict[key][:,[:Pair,:Relation,:Dist,:PoissonPvalue]])
        end
        for removepair in removepairvec
            tempdf = @subset(tempdf,.!occursin.(removepair,:Pair))
        end
        outputdf = unstack(tempdf,[:Pair,:Relation],:Dist,:PoissonPvalue)
        outputmatrix = Matrix(unstack(tempdf,[:Pair,:Relation],:Dist,:PoissonPvalue)[!,3:end])
        outputdf,outputmatrix
    end
    
end

function matricizedirectionbars(samedirdict,crossdirdict,pairvec)
    outputmatrix=Matrix
    outputdf=DataFrame
    tempdf = DataFrame(Pair=String[], Relation=String[], Frequency=Float64[])
        for pair in pairvec
            key = join(pair,"_")
            # first we clean up the statsdict tables
            push!(tempdf,[string(pair[1], "/", pair[2]),"Cis",100*(size(samedirdict[key],1)/(size(samedirdict[key],1) + size(crossdirdict[key],1)))])
            push!(tempdf,[string(pair[1], "/", pair[2]),"Trans",100*(size(crossdirdict[key],1)/(size(samedirdict[key],1) + size(crossdirdict[key],1)))])
            #push!(tempdf,[string(pair[1], "/", pair[2]),"Relation","Cis"])
            #push!(tempdf,[string(pair[1], "/", pair[2]),"Frequency",size(samedirdict[key],1)])

            #push!(tempdf,[string(pair[1], "/", pair[2]),"Relation","Trans"])
            #push!(tempdf,[string(pair[1], "/", pair[2]),"Frequency",size(crossdirdict[key],1)])
        end
        
        outputdf = unique(unstack(tempdf,:Pair,:Relation,:Frequency)) 
        outputmatrix = Matrix(unique(unstack(tempdf,:Pair,:Relation,:Frequency)[!,2:end]))
        outputdf,outputmatrix
end





