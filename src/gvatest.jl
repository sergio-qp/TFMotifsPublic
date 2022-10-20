using DataFrames, DataFramesMeta, CSV, Glob, GenomicFeatures, Plots, StatsPlots, Measures
using BioSequences, FASTX, StatsBase
using MotifScanner, GenomicIntersections

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
# this is the end of the convenience function section----------------------------------------------------

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
    gvat = rawgvat
    # gvat = @subset(rawgvat, occursin.(tfreg,:TF))
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

function loadgvatoriginal(rawgvat,tfvector)
    
    tfreg = regexmaker(tfvector)
    gvat = rawgvat
    # gvat = @subset(rawgvat, occursin.(tfreg,:TF))
    loc = split.(String.(gvat.oligo), r"[:-]")
    insertcols!(gvat, 2, :chrom => first.(loc))
    # gvat.chrom = first.(loc)
    insertcols!(gvat, 3, :start => div.(parse.(Int, getindex.(loc,2)) .+ parse.(Int, getindex.(loc,3)), 2))
    insertcols!(gvat, 4, :stop => gvat.start)
    # gvat.start = div.(parse.(Int, getindex.(loc,2)) .+ parse.(Int, getindex.(loc,3)), 2)
    # gvat.stop = gvat.start
    # insertcols!(gvat, 2, :snp => String(join([gvat.chrom, gvat.start, gvat.ref, gvat.alt],"_")))
    # gvat.snp = String(join([:chrom, :start, :ref, :alt],"_"))
    sort!(gvat, [:chrom, :start, :stop])
    uniquefilt = unique(gvat, [:chrom, :start, :stop])
    uniquefilt.ID = 1:size(uniquefilt, 1)
    inefficiency= innerjoin(gvat, uniquefilt, on=[:chrom, :start, :stop], makeunique=true)
    gvat.ID = inefficiency.ID

    gvat.TF = nkxchange.(gvat.TF);
    gvat
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


function scanmachine(snpset,mymotifs,genomelocation)
   
    uniquefilt = unique(snpset, [:chrom, :start, :stop, :ref, :alt])
    seqs = loadrefseqs(uniquefilt, genomelocation, 20); #changed to 20 cuz the gvat people are using 40bp size oligos
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

#this is for comparing gvat pbss and pwm scores
function focusscanmachine(snpset, motifvec, genomelocation)
    uniquefilt = unique(snpset, [:chrom, :start, :stop, :ref, :alt])
    seqs = loadrefseqs(uniquefilt, genomelocation, 20);
    seqs.ID = uniquefilt.ID

    #making an empty fulldf
    fullscan = motifscanall(seqs, [motifvec[1]])
    fullscan.MotifName = uppercase.(fullscan.MotifName);
    fullscan = @subset(innerjoin(fullscan,seqs, on =:ID), :ID .== 0);
    # return seqs
    for motif in motifvec
        # idk how to deal with the combined tf things, im just gonna make an if that removes any empty sets
        # pulls out only the seqs that correspond to the tf by matching snp IDs
        # return snpset
        # @show(motif.name)
        subtbl=@subset(snpset, :TF .== uppercase(motif.name))
        if !isempty(subtbl)
            motifseqs = filter(x -> x.ID ∈ subtbl.ID,seqs)
            scanned = motifscanall(motifseqs, [motif])
            scanned.MotifName = uppercase.(scanned.MotifName);
            scanned = innerjoin(scanned,motifseqs, on =:ID);
            # @show(size(scanned,1))
            append!(fullscan,scanned)
        end
    end
    # return scanned
    fullscan = innerjoin(snpset,fullscan, on =[:TF=>:MotifName,:ID])
    fullscan
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

function motifdistancer!(tf1,tf2,df)
    df[:,Symbol("motifdistance",tf1,"_",tf2)] .= div.(df[:,Symbol(tf1, "motifstart")] .+ df[:,Symbol(tf1, "motifstop")], 2) .- div.(df[:,Symbol(tf2, "motifstart")] .+ df[:,Symbol(tf2, "motifstop")], 2)
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

    mymotifs, motifvec, motifscandict, pairedpeaks
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




# a sanitychecker for makeuniqued tables, hopefully wont need anymore
# function sanitycheckermodv(df, numvec)
#     startbitvec=BitVector[]
#     stopbitvec=BitVector[]
#     # length=""

#     for num in numvec
#         tfdf = df[!, [Symbol.("start_",num),Symbol.("stop_",num)]]
#         startbitvec = df.start .> tfdf[!,1]
#         stopbitvec = df.stop .< tfdf[!,2]
#     end
#     @show sum(startbitvec), sum(stopbitvec)
#     # @show stopbitvec
# end

