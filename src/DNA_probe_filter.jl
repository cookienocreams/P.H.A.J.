using BioAlignments
using BioSequences
using Random
using DataStructures

"""
    nn_params

A structure encapsulating the thermodynamic parameters for DNA nearest-neighbor (NN) interactions. 

The values are based on SantaLucia's 1997 paper titled 
"Thermodynamics and NMR of internal G·T mismatches in DNA." 
[DOI](https://doi.org/10.1021/bi962590c).

# Parameters

The struct contains entropy (ΔS, in cal/mol·K) and enthalpy (ΔH, in kcal/mol) values for 
each of the NN pairs:

| NN Pair | ΔS                | ΔH                |
|:-------:|:-----------------:|:-----------------:|
|   AA    | AA_delta_s        | AA_delta_h        |
|   AC    | AC_delta_s        | AC_delta_h        |
|   AG    | AG_delta_s        | AG_delta_h        |
|   AT    | AT_delta_s        | AT_delta_h        |
|   CA    | CA_delta_s        | CA_delta_h        |
|   CC    | CC_delta_s        | CC_delta_h        |
|   CG    | CG_delta_s        | CG_delta_h        |
|   GA    | GA_delta_s        | GA_delta_h        |
|   GC    | GC_delta_s        | GC_delta_h        |
|   TA    | TA_delta_s        | TA_delta_h        |

# Notes
This struct is essential for functions aiming to compute the thermodynamics of DNA sequences 
or structures.
"""
@kwdef struct nn_params
    AA_delta_s::Float64 = -22.15461879548201
    AA_delta_h::Float64 = -7.972759800039155
    AC_delta_s::Float64 = -20.20944930135789
    AC_delta_h::Float64 = -7.6063649173383086
    AG_delta_s::Float64 = -20.85324697548337
    AG_delta_h::Float64 = -7.7250305886539525
    AT_delta_s::Float64 = -19.931105027077287
    AT_delta_h::Float64 = -7.069032819297647
    CA_delta_s::Float64 = -21.383943453258333
    CA_delta_h::Float64 = -7.9794539347093165
    CC_delta_s::Float64 = -19.677205419525077
    CC_delta_h::Float64 = -7.612190594454272
    CG_delta_s::Float64 = -24.335240882791073
    CG_delta_h::Float64 = -9.348781072782547
    GA_delta_s::Float64 = -21.617116510641825
    GA_delta_h::Float64 = -8.058362530689527
    GC_delta_s::Float64 = -23.896038387305747
    GC_delta_h::Float64 = -9.457805766372232
    TA_delta_s::Float64 = -19.12903239340641
    TA_delta_h::Float64 = -6.6574837353237415
end

"""
    mutable struct Region
        start::Int
        stop::Int
    end

A mutable structure that represents a region in a sequence. `start` and `stop` are the 
1-based indices of the start and end of the region, respectively.
"""
mutable struct Region
    start::Int
    stop::Int
end

"""
    prpy_sequence(sequence::AbstractString)::String

Convert a given DNA sequence into purine (A/G represented by 'R') 
and pyrimidine (C/T represented by 'Y') codes.

# Arguments
- `sequence::AbstractString`: A DNA sequence containing nucleotide bases.

# Returns
- A string where purines are replaced with 'R' and pyrimidines with 'Y'.
"""
function prpy_sequence(sequence::AbstractString)::String
    result = IOBuffer()
    for base in sequence
        if base in ('A', 'G')
            write(result, 'R')
        else
            write(result, 'Y')
        end
    end
    return String(take!(result))
end

"""
    sequence_thermodynamic_sum(delta_h::Vector{Float64}
                                    , delta_s::Vector{Float64}
                                    , nn_pairs_list::Vector{String}
                                    , sequence_nn_list::Vector{String}
                                    )

Calculate the total enthalpy and entropy of each sequence.

# Arguments
- `delta_h::Vector{Float64}`: A vector containing enthalpy values for nearest-neighbor pairs.
- `delta_s::Vector{Float64}`: A vector containing entropy values for nearest-neighbor pairs.
- `nn_pairs_list::Vector{String}`: A vector containing nearest-neighbor pairs.
- `sequence_nn_list`: A list of nearest-neighbor pairs for the given sequence.

# Returns
- Total enthalpy and entropy for the sequence.
"""
function sequence_thermodynamic_sum(delta_h::Vector{Float64}
                                    , delta_s::Vector{Float64}
                                    , nn_pairs_list::Vector{String}
                                    , sequence_nn_list::Vector{String}
                                    )
    sequence_total_dh, sequence_total_ds = 0.0, 0.0

    for (count, nn_pair) in enumerate(nn_pairs_list)
        local_sum_dh, local_sum_ds = 0.0, 0.0
        for nn in sequence_nn_list
            if occursin(nn, nn_pair)
                local_sum_dh += delta_h[count]
                local_sum_ds += delta_s[count]
            end
        end

        sequence_total_dh += local_sum_dh
        sequence_total_ds += local_sum_ds
    end

    return sequence_total_dh, sequence_total_ds
end

"""
    gc_content(sequence::LongSequence{DNAAlphabet{4}})::Float64

Calculate the GC content of a given DNA sequence.

# Arguments
- `sequence::LongSequence{DNAAlphabet{4}}`: A DNA sequence containing nucleotide bases.

# Returns
- The percentage of G and C bases in the given sequence.
"""
function gc_content(sequence::LongSequence{DNAAlphabet{4}})::Float64
    gc_bases = Set([DNA_G, DNA_C])
    gc_count = count(base -> base in gc_bases, sequence)
    return gc_count / length(sequence)
end

"""
    find_homodimers!(seq::LongSequence{DNAAlphabet{4}}
                        , window_size::Int
                        , homodimers::Vector{LongSequence{DNAAlphabet{4}}}
                        )

Find all homodimers of a given size in a DNA sequence.

# Arguments
- `seq::LongSequence{DNAAlphabet{4}}`: The DNA sequence to search.
- `window_size::Int`: The size of the homodimers to search for.
- `homodimers::Vector{LongSequence{DNAAlphabet{4}}}`: A vector to hold the sequences of the homodimers found.

# Returns
- The `Vector` containing the sequences of the homodimers of size = `window_size` found in the sequence.

# Notes
- The homodimer `Vector` is modified in place to avoid unecessary memory allocations.
"""
function find_homodimers!(sequence::LongSequence{DNAAlphabet{4}}
                        , window_size::Int
                        , homodimers::Vector{LongSequence{DNAAlphabet{4}}}
                        )

    # Slide the window along the sequence
    for i in 1:(length(sequence) - window_size + 1)
        # Extract the windowed sequence and its reverse complement
        windowed_seq = sequence[i:(i + window_size - 1)]
        rev_comp = BioSequences.reverse_complement(windowed_seq)

        # Check if the reverse complement appears elsewhere in the sequence
        if  BioSequences.occursin(ExactSearchQuery(rev_comp), sequence[(i + window_size):end])
            # If it does, add it to the vector
            push!(homodimers, windowed_seq)
        end
    end

    return homodimers
end

"""
    find_homodimers(seq::LongSequence{DNAAlphabet{4}})

Find the largest homodimers in a DNA sequence.

# Arguments
- `seq::LongSequence{DNAAlphabet{4}}`: The DNA sequence to search.

# Returns
- A `Vector` containing the sequences of the largest homodimers found in the sequence.
"""
function find_homodimers(sequence::LongSequence{DNAAlphabet{4}})

    # Initialize the output dictionary
    homodimers = Vector{LongSequence{DNAAlphabet{4}}}()

    for num in eachindex(sequence)
        # Decrease window size until a match is found
        window_size = length(sequence) - num
        find_homodimers!(sequence, window_size, homodimers)
        # If match is found, i.e., sequences have been added to `homodimers`, return vector
        if !isempty(homodimers)
            return homodimers
        end
    end

    return homodimers
end

"""
    filter_probes(probes::Vector{String}, temperature_threshold::Integer)

Filter probes based on the adjusted melting temperature with respect to a temperature threshold. 
Calculates the adjusted melting temperature for pairs of sequences and returns those sequences 
with temperatures below the given threshold. Additionally, calculates and returns Gibbs free energy 
for the aligned region of each sequence pair.

# Arguments
- `probes::Vector{String}`: A vector of DNA sequences (probes) to be filtered.
- `temperature_threshold::Integer`: A temperature threshold (in Celsius) to filter the probes.

# Returns
- `promising_probes`: A set of sequences with adjusted melting temperatures below the temperature threshold.
- `match_info`: A dictionary mapping sequences to a tuple of Gibbs free energy for the aligned region and 
   the maximum Gibbs free energy for the entire sequence.

# Notes
The function considers various thermodynamic and environmental parameters, such as mono- and divalent ion concentrations, 
to calculate adjusted melting temperatures and Gibbs free energy for each sequence pair.
"""
function filter_probes(probes::Vector{AbstractString}, temperature_threshold::Integer)
    gas_constant = 1.98720425864083 # cal*K−1*mol−1
    monovalent = 50.0 # mM
    mg = 2.0 # mM
    dntps = 0.0
    oligo_c = .25 * 1e-6 # μM

    # Reuse the same Region object for all comparisons
    region = Region(0, 0)

    # Get all nearest neighbor thermodynamic values by creating instance of nn_params struct
    nn = nn_params()

    promising_probes = Set{AbstractString}()

    number_of_probes = length(probes)
    scoremodel = AffineGapScoreModel(
               match=5,
               mismatch=-4,
               gap_open=-5,
               gap_extend=-3
           )
    
    # Loop through the input sequences and calculate each one's delta G
    # Loops are designed to ensure that each pair of sequences is compared only once
    for i in 1:number_of_probes

        # Get probe percent GC
        gc = gc_content(LongDNA{4}(probes[i]))
        # Skip probes that have a percent GC more than 60% or less than 50%
        if gc > .6 || gc < .5
            continue
        end

        # Homo dimer analysis
        homo_dimers = find_homodimers(LongDNA{4}(probes[i]))
        homo_dimer_delta_g_values = Vector{Float64}()
        if !isempty(homo_dimers)
            for homo_dimer in homo_dimers
                homo_dimer_delta_g, _, _ = calculate_thermodynamic_parameters(string(homo_dimer), nn)
                push!(homo_dimer_delta_g_values, homo_dimer_delta_g)
            end
        end

        # Skip probes that have a stretch of complementary bases with a delta G below -10 kcal/mol
        if any(<(-10), homo_dimer_delta_g_values)
            continue
        end

        for j in (i+1):number_of_probes

            # Hetero-dimer analysis
            alignment_result = pairalign(SemiGlobalAlignment()
                                        , LongDNA{4}(probes[i])
                                        , BioSequences.reverse_complement(LongDNA{4}(probes[j]))
                                        , scoremodel
                                        )
            alignment_anchors = alignment(alignment_result).a.aln.anchors
            longest_aligned_region!(region, alignment_anchors)
            aligned_seq = probes[i][region.start+1:region.stop]

            # Skip probes that have more than 5 bp of alignment to other probes
            if length(aligned_seq) > 5
                continue
            end

            tm_okay = check_probe_tm(probes[i], temperature_threshold, monovalent, nn, oligo_c, dntps, mg)

            # Continue if probe Tm is above the threshold
            if tm_okay
                # Calculate probe enthalpy and entropy values
                _, aligned_total_dh, aligned_total_ds = calculate_thermodynamic_parameters(aligned_seq, nn)

                # Determine the melting temperature of the sequence at 1M monovalent salt, 
                # still needs adjustment to account for actual salt concentration + Mg2+ and dNTPs
                aligned_melting_temperature = (1000 * aligned_total_dh) / (aligned_total_ds + (gas_constant * (log(oligo_c)))) - 273.15
 
                # Get complementary bases percent GC
                aligned_gc = gc_content(LongDNA{4}(aligned_seq))

                # Adjust melting temperature for the specific reactiion conditions
                salt_adj_aligned_melting_temperature = probe_tm_salt_correction(aligned_melting_temperature, monovalent, mg, dntps, aligned_gc, length(aligned_seq))

                # If aligned region has a melting temperature less than 25C, add to list
                if salt_adj_aligned_melting_temperature < 25
                    push!(promising_probes, probes[i])
                end
            end
        end
    end

    return promising_probes
end

"""
    check_probe_tm(probe::String
                 , temperature_threshold::Integer
                 , monovalent::Float64
                 , nn::nn_params
                 , oligo_c::Float64
                 , dntps::Float64
                 , mg::Float64)

Check if all slices of a DNA probe sequence have melting temperatures (`Tm`) above a 
set threshold. The function splits the probe into roughly two equal slices and calculates 
the `Tm` for each slice considering corrections for monovalent ions, magnesium ions, and dNTPs.

# Arguments
- `probe::String`: The DNA probe sequence to check.
- `temperature_threshold::Integer`: The minimum required melting temperature in degrees Celsius.
- `monovalent::Int64`: Concentration of monovalent ions in mM.
- `nn::nn_params`: A structure containing the thermodynamic parameters for all possible nearest-neighbor pairs.
- `oligo_c::Float64`: Oligo concentration in μM.
- `dntps::Float64`: Concentration of dNTPs in mM.
- `mg::Float64`: Concentration of magnesium ions in mM.

# Returns
- `Boolean`: `true` if all slices have a `Tm` above the temperature threshold, otherwise `false`.

# Notes
This function aims to ensure that all subsections of a DNA probe sequence have sufficient 
stability (i.e., high melting temperature) to maintain a double-stranded configuration 
under specific conditions. It's crucial when designing probes for applications where
uniform hybridization across the probe is desired.

# Examples
```julia
probe = "AGTCGATCGA"
threshold = 60
monovalent_concentration = 50
nn_parameters = nn_params()
oligo_concentration = 2 * 1e-6
dNTP_concentration = 0
mg_concentration = 2
result = check_probe_tm(probe
                        , threshold
                        , monovalent_concentration
                        , nn_parameters
                        , oligo_concentration
                        , dNTP_concentration
                        , mg_concentration)
"""
function check_probe_tm(probe::String
                        , temperature_threshold::Integer
                        , monovalent::Float64
                        , nn::nn_params
                        , oligo_c::Float64
                        , dntps::Float64
                        , mg::Float64
                        )

    # Check that the probe melting temperature is above the set temperature threshold
    slice_size = Int(length(probe) / 2)
    slice_tms = Vector{Float64}()

    probe_slices = slice_probe(probe, slice_size)

    # Iterate through each probe slice
    for slice in probe_slices
        sequence_gc = gc_content(LongDNA{4}(slice))
        _, total_dh, total_ds = calculate_thermodynamic_parameters(slice, monovalent, nn)
        melting_temperature = (1000 * total_dh) / (total_ds + (gas_constant * (log(oligo_c)))) - 273.15
        salt_adj_melting_temperature = probe_tm_salt_correction(melting_temperature, monovalent, mg, dntps, sequence_gc, length(slice))
        push!(slice_tms, salt_adj_melting_temperature)

        # Return true if all slices have a melting temperature above the threshold
        if all(>(temperature_threshold), slice_tms)
            return true
        else 
            return false
        end
    end
end

"""
    slice_probe(probe::String, slice_size::Int) -> Vector{String}

Get substring slices from a probe of size `slice_size`.

# Arguments
- `probe`: A string containing the probe sequence to be split.
- `slice_size`: The size of the substring to created. Each probe is divided into
  chunks of this size (from left to right, and without overlap between chunks).

# Returns
Returns the substrings of the input probe sequence.

# Example
```julia
probe = "CGTGCGCCACTAGACTTGGCAAGGCGTGGAACCGATACCTGCTACCGTGTTAGCAACAAACAGCTATCAACACAGCCATG"
slice_size = 40
probe_slices = slice_probe(probe, slice_size)
2-element Vector{String}:
 "CGTGCGCCACTAGACTTGGCAAGGCGTGGAACCGATACCT"
 "GCTACCGTGTTAGCAACAAACAGCTATCAACACAGCCATG"
```
"""
function slice_probe(probe::String, slice_size::Int)
    slices = Vector{String}()
    for i in 1:slice_size:length(probe)
        substring = probe[i:min(i+slice_size-1, end)]
        push!(slices, substring)
    end

    return slices
end

"""
    calculate_thermodynamic_parameters(sequence::AbstractString, monovalent::Float64, nn::nn_params)

Calculate the Gibbs free energy (`ΔG`), total enthalpy (`ΔH`), and total entropy (`ΔS`) of 
a given DNA sequence. The computation is based on the nearest-neighbor (NN) model, taking 
into account the sequence's content and a provided monovalent ion concentration.

# Arguments
- `sequence::AbstractString`: A DNA sequence for which thermodynamic parameters will be calculated.
- `monovalent::Int64`: The concentration of monovalent ions in mM.
- `nn::nn_params`: A structure containing the thermodynamic parameters for all possible nearest-neighbor pairs.

# Returns
- `aligned_delta_g`: The Gibbs free energy (`ΔG`) of the sequence.
- `sequence_dh_total`: The total enthalpy (`ΔH`) of the sequence, considering both internal 
nearest-neighbor interactions and terminal base pair corrections.
- `sequence_total_ds`: The total entropy (`ΔS`) of the sequence based on its nearest-neighbor content.

# Notes
The function uses nearest-neighbor parameters and end compensation (terminal base pair 
corrections) to compute the thermodynamics of the DNA sequence. The computed Gibbs free 
energy is especially important for predicting the stability of DNA duplexes.

# Examples
```julia
seq = "AGTCGA"
monovalent_concentration = 50
nn_parameters = nn_params()
ΔG, ΔH, ΔS = calculate_thermodynamic_parameters(seq, monovalent_concentration, nn_parameters)
"""
function calculate_thermodynamic_parameters(sequence::AbstractString, nn::nn_params)

    # List of all the nearest-neighbor pairs
    nn_pairs_list::Vector{String} = ["AA, TT", "AC, TG", "AG, TC", "AT", "CA, GT"
                                    , "CC, GG", "CG", "GA, CT", "GC", "TA"
    ]

    sequence_length = length(sequence)

    # Split input sequences into doublets for easy NN matching below
    sequence_NN_list = [sequence[i:i+1] for i in 1:sequence_length-1]

    # Lists of each nearest neighbor's entropy and enthalpy values that will be used to calculate the thermodynamic values of each sequence
    delta_s = [nn.AA_delta_s, nn.AC_delta_s, nn.AG_delta_s, nn.AT_delta_s, nn.CA_delta_s, nn.CC_delta_s, nn.CG_delta_s, nn.GA_delta_s, nn.GC_delta_s, nn.TA_delta_s]
    delta_h = [nn.AA_delta_h, nn.AC_delta_h, nn.AG_delta_h, nn.AT_delta_h, nn.CA_delta_h, nn.CC_delta_h, nn.CG_delta_h, nn.GA_delta_h, nn.GC_delta_h, nn.TA_delta_h]

    sequence_total_dh, sequence_total_ds = sequence_thermodynamic_sum(delta_h, delta_s, nn_pairs_list, sequence_NN_list)
    # Calculate Gibb's free energy
    aligned_delta_g = sequence_total_dh - (273.15 * (sequence_total_ds/1000))

    return aligned_delta_g, sequence_total_dh, sequence_total_ds
end

function calculate_thermodynamic_parameters(sequence::AbstractString, monovalent::Float64, nn::nn_params)

    # List of all the nearest-neighbor pairs
    nn_pairs_list::Vector{String} = ["AA, TT", "AC, TG", "AG, TC", "AT", "CA, GT"
                                    , "CC, GG", "CG", "GA, CT", "GC", "TA"
    ]

    sequence_length = length(sequence)
    sequence_prpy_str = prpy_sequence(sequence)

    # Capture the first and last pair of bases in each sequence so that the energies associated with opening the helix can be calculated
    sequence_initial_bases, sequence_terminal_bases = sequence_prpy_str[1:2], sequence_prpy_str[sequence_length-1:end]

    # Split input sequences into doublets for easy NN matching below
    sequence_NN_list = [sequence[i:i+1] for i in 1:sequence_length-1]

    # Lists of each nearest neighbor's entropy and enthalpy values that will be used to calculate the thermodynamic values of each sequence
    delta_s = [nn.AA_delta_s, nn.AC_delta_s, nn.AG_delta_s, nn.AT_delta_s, nn.CA_delta_s, nn.CC_delta_s, nn.CG_delta_s, nn.GA_delta_s, nn.GC_delta_s, nn.TA_delta_s]
    delta_h = [nn.AA_delta_h, nn.AC_delta_h, nn.AG_delta_h, nn.AT_delta_h, nn.CA_delta_h, nn.CC_delta_h, nn.CG_delta_h, nn.GA_delta_h, nn.GC_delta_h, nn.TA_delta_h]

    # End Compensation parameters adjusted for sequence length and changes in monovalent ion concentration
    YY_constant_one = (-.0235 * sequence_length) + .1273
    YY_constant_two = (.1639 * sequence_length) - .895
    YR_constant_one = (-.296 * log(sequence_length)) + .5058
    YR_constant_two = (2.0303 * log(sequence_length)) - 3.4594
    YY_length_adjust_eq = YY_constant_one * log(monovalent) + YY_constant_two
    YR_length_adjust_eq = YR_constant_one * log(monovalent) + YR_constant_two
    YY_h,RY_h = 0.959547884222969 - YY_length_adjust_eq, 0.849386567650066 - YR_length_adjust_eq

    sequence_total_dh, sequence_total_ds = sequence_thermodynamic_sum(delta_h, delta_s, nn_pairs_list, sequence_NN_list)
    sequence_dh_init = (if sequence_initial_bases in ["YY","RR"] sequence_total_dh + YY_h else sequence_total_dh + RY_h end)
    sequence_dh_total = (if sequence_terminal_bases in ["YY","RR"] sequence_dh_init + YY_h else sequence_dh_init + RY_h end)

    # Calculate Gibb's free energy
    aligned_delta_g = sequence_dh_total - (273.15 * (sequence_total_ds/1000))

    return aligned_delta_g, sequence_dh_total, sequence_total_ds
end

"""
    longest_aligned_region!(region::Region, alignment_anchors::Vector{AlignmentAnchor})

Modifies the `region` object in place to represent the longest aligned region in the given 
alignment `aln`.

# Arguments
- `region::Region`: A `Region` object that will be modified in place. After the function call, 
`region.start` and `region.stop` will be the start and end indices of the longest aligned region in `aln`.
- `alignment_anchors::Vector{AlignmentAnchor}`: The alignment anchors to analyze.

# Examples

julia region = Region(0, 0) aln = alignment(alignment_result) longest_aligned_region!(region, aln)
println("The longest aligned region is from index(region.start) to index(region.stop)") 
"""
function longest_aligned_region!(region::Region, alignment_anchors::Vector{AlignmentAnchor})
    longest_region = (start=0, stop=0)
    current_region = (start=0, stop=0)

    # Iterate over the anchors in the alignment
    for index in eachindex(alignment_anchors)
        anchor = alignment_anchors[index]
        # If the operation is a match, extend the current region
        if anchor.op == OP_SEQ_MATCH
            current_region = (start=current_region.start, stop=anchor.seqpos)
            # If the current region is longer than the longest region found so far, update the longest region
            if current_region.stop - current_region.start > longest_region.stop - longest_region.start
                region.start = current_region[1]
                region.stop = current_region[2]
                longest_region = current_region
            end
        else
            current_region = (start=anchor.seqpos, stop=anchor.seqpos)
        end
    end
end

"""
    probe_tm_salt_correction(probe_melting_temperature, monovalent, magnesium, dNTPs, probe_gc, probe_length)

Adjust the provided melting temperature for a DNA probe based on salt and magnesium concentrations in the buffer.

Uses equations from Owczarzy et al. (2008) to modify the melting temperature according to monovalent ion, 
magnesium, and dNTP concentrations. 

Reference: Owczarzy, R., et al. (2008). Biochemistry, https://doi.org/10.1021/bi702363u.

# Arguments
- `probe_melting_temperature::Float64`: Initial melting temperature (in Celsius) calculated at 1 M salt concentration.
- `monovalent::Float64`: Concentration of monovalent ions (e.g., Na+ or K+).
- `magnesium::Float64`: Concentration of magnesium ions.
- `dNTPs::Float64`: Concentration of dNTPs.
- `probe_gc::Float64`: GC content of the probe as a fraction (0.0-1.0).
- `probe_length::Int`: Length of the DNA probe.

# Returns
- `Float64`: Adjusted melting temperature (in Celsius).
"""
function probe_tm_salt_correction(probe_melting_temperature::Float64
                                , monovalent::Float64
                                , magnesium::Float64
                                , dNTPs::Float64
                                , probe_gc::Float64
                                , probe_length::Int
                                ) 
    # Constants and initial adjustments
    ka = 3e4  # Association constant for Mg2+-dNTP complex
    mon = monovalent * 1e-3 / 2  # Convert to mol/L and adjust for counterions
    mg_adj = magnesium * 1e-3
    dntps = dNTPs * 1e-3
    mg = (-(ka * dntps - ka * mg_adj + 1.0) + sqrt((ka * dntps - ka * mg_adj + 1.0)^2 + 4.0 * ka * mg_adj)) / (2.0 * ka)
    R = sqrt(mg) / mon

    # Equations based on the primary factor affecting melting temperature
    if mon == 0 || R > 6.0
        return compute_tm_adjustment(probe_melting_temperature, probe_gc, probe_length, mg)
    elseif R < 0.22
        mon_tm_factor = monovalent > 900 ? 0 : (4.29e-5 * probe_gc - 3.95e-5) * log(mon) + 9.40e-6 * log(mon)^2
        return compute_tm_from_salt(probe_melting_temperature, mon_tm_factor)
    else
        return compute_tm_adjustment(probe_melting_temperature, probe_gc, probe_length, mg, mon)
    end
end

# Helper function to compute melting temperature adjustment
function compute_tm_adjustment(tm, gc, len, mg, mon=nothing)
    a, b, c, d, e, f, g = 3.919e-5, -2.88e-5, 3.603e-5, 2.322e-5, -3.507e-4, 4.711e-4, 6.52e-5 
    const_a,const_b,const_c,const_d,const_e,const_f,const_g,const_h = -0.1156,-2.0725,-0.1445,6.247e-3,6.131e-3,0.0314,0.5308,4.563e-3

    if mon !== nothing
        a *= const_a * log(mg) - (const_b * sqrt(mon) * log(mon))
        d *= const_c * log(mg) - const_d * log(mon) - const_e * log(mon)^2
        g *= const_f * log(mg) - const_g * log(mon) + const_h * log(mon)^3
    end
    
    salty = (1 / (tm + 273.15)) + a + (b * log(mg)) + 
            (gc * (c + d * log(mg))) + (1 / (2.0 * (len - 1))) * 
            (e + f * log(mg) + g * log(mg)^2)

    return round((1 / salty) - 273.15, digits=1)
end

# Helper function for computing temperature from salt factor
function compute_tm_from_salt(tm, mon_tm_factor)
    salty = (1 / (tm + 273.15)) + mon_tm_factor
    return round((1 / salty) - 273.15, digits=1)
end

function write_output(probes_dict::Dict{String, String}, promising_probes::Set{AbstractString})
    open("promising_probes.fa", "w") do output_file
        for probe in promising_probes
            write(output_file, string(">", probes_dict[probe], "\n", probe, "\n"))
        end
    end
end

function write_fasta(promising_probes::Vector{String})
    open("all_probes.fa", "w") do output_file
        for (index, probe) in enumerate(promising_probes)
            write(output_file, string(">", "Probe", index, "\n", probe, "\n"))
        end
    end
end

function read_fasta(input_fasta::String)
    probes = Dict{String, String}()
    open(input_fasta, "r") do input_file
        for line in eachline(input_file)
            if startswith(line, ">")
                probes[readline(input_file)] = line[2:end]
            end
        end
    end

    return probes
end

function main()
    test_probes = [randstring("ACGT", 80) for _ in 1:10000]
    write_fasta(test_probes)
    probes_dict = read_fasta("all_probes.fa")
    probes = collect(keys(probes_dict))
    promising_probes = filter_probes(probes, 65)
    write_output(probes_dict, promising_probes)
end

main()
