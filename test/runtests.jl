using dna_probe_filter
using Test
using BioSequences

@testset "purine_pyrimidine" begin
    @testset "standard prpy" begin
        probe = "TTTAAACCTCAGAGTGATACGCTCAGCTACCTGCAGCCGCCCCAACTAGTACGGCCATTTAAGGCGGAGTTATTTACATT"
        result = dna_probe_filter.prpy_sequence(probe)
        expected = "YYYRRRYYYYRRRRYRRYRYRYYYRRYYRYYYRYRRYYRYYYYRRYYRRYRYRRYYRYYYRRRRYRRRRYYRYYYRYRYY"
        @test result == expected
    end
    
    @testset "all purine" begin
        probe = "AAAAAAAAAAAAAAAAAA"
        result = dna_probe_filter.prpy_sequence(probe)
        expected = "RRRRRRRRRRRRRRRRRR"
        @test result == expected
    end
    
    @testset "all pyrimidine" begin
        probe = "TTTTTTTTTTTTTTTTTT"
        result = dna_probe_filter.prpy_sequence(probe)
        expected = "YYYYYYYYYYYYYYYYYY"
        @test result == expected
    end
end

@testset "thermodynamic_sum" begin
    @testset "standard thermodynamic sum" begin
        delta_h = [-7.972759800039155, -7.6063649173383086, -7.7250305886539525, -7.069032819297647, -7.9794539347093165, -7.612190594454272, -9.348781072782547, -8.058362530689527, -9.457805766372232, -6.6574837353237415]
        delta_s = [-22.15461879548201, -20.20944930135789, -20.85324697548337, -19.931105027077287, -21.383943453258333, -19.677205419525077, -24.335240882791073, -21.617116510641825, -23.896038387305747, -19.12903239340641]
        nn_pairs_list = ["AA, TT", "AC, TG", "AG, TC", "AT", "CA, GT"
                                    , "CC, GG", "CG", "GA, CT", "GC", "TA"
        ]
        sequence_NN_list = ["AA", "TC", "AA", "CA", "AT", "CA"]
        sum = dna_probe_filter.sequence_thermodynamic_sum(delta_h, delta_s, nn_pairs_list, sequence_NN_list)
        expected = (-46.69849087744854, -127.86147650004135)
        @test sum == expected
    end
    
    @testset "empty sequence list" begin
        delta_h = [-7.972759800039155, -7.6063649173383086, -7.7250305886539525, -7.069032819297647, -7.9794539347093165, -7.612190594454272, -9.348781072782547, -8.058362530689527, -9.457805766372232, -6.6574837353237415]
        delta_s = [-22.15461879548201, -20.20944930135789, -20.85324697548337, -19.931105027077287, -21.383943453258333, -19.677205419525077, -24.335240882791073, -21.617116510641825, -23.896038387305747, -19.12903239340641]
        nn_pairs_list = ["AA, TT", "AC, TG", "AG, TC", "AT", "CA, GT"
                                    , "CC, GG", "CG", "GA, CT", "GC", "TA"
        ]
        sequence_NN_list::Vector{String} = []
        sum = dna_probe_filter.sequence_thermodynamic_sum(delta_h, delta_s, nn_pairs_list, sequence_NN_list)
        expected = (0.0, 0.0)
        @test sum == expected
    end
    
    @testset "all N bases" begin
        delta_h = [-7.972759800039155, -7.6063649173383086, -7.7250305886539525, -7.069032819297647, -7.9794539347093165, -7.612190594454272, -9.348781072782547, -8.058362530689527, -9.457805766372232, -6.6574837353237415]
        delta_s = [-22.15461879548201, -20.20944930135789, -20.85324697548337, -19.931105027077287, -21.383943453258333, -19.677205419525077, -24.335240882791073, -21.617116510641825, -23.896038387305747, -19.12903239340641]
        nn_pairs_list = ["AA, TT", "AC, TG", "AG, TC", "AT", "CA, GT"
                                    , "CC, GG", "CG", "GA, CT", "GC", "TA"
        ]
        sequence_NN_list = ["NN", "NN", "NN", "NN", "NN", "NN"]
        sum = dna_probe_filter.sequence_thermodynamic_sum(delta_h, delta_s, nn_pairs_list, sequence_NN_list)
        expected = (0.0, 0.0)
        @test sum == expected
    end

    @testset "delta s vs delta h mismatch" begin
        delta_h = [-7.972759800039155, -7.6063649173383086, -7.7250305886539525, -7.069032819297647, -7.9794539347093165, -7.612190594454272, -9.348781072782547, -8.058362530689527, -9.457805766372232, -6.6574837353237415]
        delta_s = [-22.15461879548201, -20.85324697548337, -19.931105027077287, -21.383943453258333, -19.677205419525077, -24.335240882791073, -21.617116510641825, -23.896038387305747, -19.12903239340641]
        nn_pairs_list = ["AA, TT", "AC, TG", "AG, TC", "AT", "CA, GT"
                                    , "CC, GG", "CG", "GA, CT", "GC", "TA"
        ]
        sequence_NN_list = ["AA", "TC", "AA", "CA", "AT", "CA"]
        sum = dna_probe_filter.sequence_thermodynamic_sum(delta_h, delta_s, nn_pairs_list, sequence_NN_list)
        expected = (-46.69849087744854, -124.97869691034978) # delta s value is low
        @test sum == expected
    end
end

@testset "percent_gc" begin
    @testset "standard gc" begin
        probe = dna"TTTAAACCTCAGAGTGATACGCTCAGCTACCTGCAGCCGCCCCAACTAGTACGGCCATTTAAGGCGGAGTTATTTACATT"
        result = dna_probe_filter.gc_content(probe)
        expected = 0.475
        @test result == expected
    end
    
    @testset "no GC" begin
        probe = dna"AAAAAAAAAAAAAAAAAA"
        result = dna_probe_filter.gc_content(probe)
        expected = 0.0
        @test result == expected
    end
    
    @testset "all GC" begin
        probe = dna"GCGCGCGCGCGCGCGCGC"
        result = dna_probe_filter.gc_content(probe)
        expected = 1.0
        @test result == expected
    end

    @testset "lowercse all GC" begin
        probe = dna"gcgcgcgcgcgcgcgcgc"
        result = dna_probe_filter.gc_content(probe)
        expected = 1.0
        @test result == expected
    end
end

@testset "suffix testing" begin
    @testset "standard suffix array" begin
        probe = "TTTAAACCTC"
        result = dna_probe_filter.build_suffix_array(probe)
        expected = [4, 5, 6, 10, 7, 8, 3, 9, 2, 1]
        @test result == expected
    end
    
    @testset "same base suffix array" begin
        probe = "AAAAAAAAAA"
        result = dna_probe_filter.build_suffix_array(probe)
        expected = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
        @test result == expected
    end
    
    @testset "empty suffix array" begin
        probe = ""
        result = dna_probe_filter.build_suffix_array(probe)
        expected = Int64[]
        @test result == expected
    end
end

@testset "lcp testing" begin
    @testset "standard lcp" begin
        probe = "TTTAAACCTC"
        suffix_array = [4, 5, 6, 10, 7, 8, 3, 9, 2, 1]
        result = dna_probe_filter.build_longest_common_prefix(probe, suffix_array)
        expected = [2, 1, 0, 1, 1, 0, 1, 1, 2]
        @test result == expected
    end
    
    @testset "same base suffix array" begin
        probe = "AAAAAAAAAA"
        suffix_array = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
        result = dna_probe_filter.build_longest_common_prefix(probe, suffix_array)
        expected = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        @test result == expected
    end
end

@testset "homodimer testing" begin
    @testset "standard homodimer" begin
        probe = dna"TTTAAACCTC"
        result = dna_probe_filter.find_homodimers(probe)
        expected = ["AAA", "A", "TAAA", "T", "TTAAA", "TT", "TTTAAA"]
        @test result == expected
    end
    
    @testset "same base suffix array" begin
        probe = dna"AAAAAAAAAA"
        result = dna_probe_filter.find_homodimers(probe)
        expected = String[]
        @test result == expected
    end
end

@testset "tm check testing" begin
    @testset "standard tm probe" begin
        probe = "AATTATGACTGGGAAAGTAAACCGCCTCCACGTAAGCAAGGAAGGCATTCCAATTGTCGAACGGACTGAAGTTTCGGATA"
        threshold = 65
        monovalent_concentration = 50.0
        oligo_concentration = 2 * 1e-6
        dNTP_concentration = 0.0
        mg_concentration = 2.0
        nn_parameters = dna_probe_filter.nn_params()
        result = dna_probe_filter.check_probe_tm(probe
                                                , threshold
                                                , monovalent_concentration
                                                , nn_parameters
                                                , oligo_concentration
                                                , dNTP_concentration
                                                , mg_concentration
        )
        expected = true
        @test result == expected
    end
    
    @testset "below tm cutoff" begin
        probe = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        threshold = 65
        monovalent_concentration = 50.0
        oligo_concentration = 2 * 1e-6
        dNTP_concentration = 0.0
        mg_concentration = 2.0
        nn_parameters = dna_probe_filter.nn_params()
        result = dna_probe_filter.check_probe_tm(probe
                                                , threshold
                                                , monovalent_concentration
                                                , nn_parameters
                                                , oligo_concentration
                                                , dNTP_concentration
                                                , mg_concentration
        )
        expected = false
        @test result == expected
    end
end

@testset "slice_testing" begin
    @testset "standard slicing" begin
        probe = "CGTGCGCCACTAGACTTGGCAAGGCGTGGAACCGATACCTGCTACCGTGTTAGCAACAAACAGCTATCAACACAGCCATG"
        slice_size = 40
        result = dna_probe_filter.slice_probe(probe, slice_size)
        expected = ["CGTGCGCCACTAGACTTGGCAAGGCGTGGAACCGATACCT", "GCTACCGTGTTAGCAACAAACAGCTATCAACACAGCCATG"]
        @test result == expected
    end
    
    @testset "smaller probe" begin
        probe = "CGTGCGCCACTAGACT"
        slice_size = 10
        result = dna_probe_filter.slice_probe(probe, slice_size)
        expected = ["CGTGCGCCAC", "TAGACT"]
        @test result == expected
    end
    
    @testset "slice size larger than probe" begin
        probe = "CGTGCG"
        slice_size = 10
        result = dna_probe_filter.slice_probe(probe, slice_size)
        expected = ["CGTGCG"]
        @test result == expected
    end
end

