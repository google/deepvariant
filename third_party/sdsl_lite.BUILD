# Description:
# (https://github.com/vgteam/sdsl-lite)
# SDSL - Succinct Data Structure Library (vgteam fork)

package(
    default_visibility = ["//visibility:public"],
    features = [
        "-layering_check",
        "-parse_headers",
    ],
)

licenses(["restricted"])

exports_files(["COPYING"])

include_sdsl_lite = "include"

includes = [
    include_sdsl_lite,
    ".",
]

LIBRARY_COPTS = [
    "-w",  # turn off all warnings
]

SRCS = glob(
    ["lib/*cpp"],
    exclude = [
        "lib/lcp_support_tree.cpp",
        "lib/construct_lcp_helper.cpp",
        "lib/construct_lcp.cpp",
        "lib/louds_tree.cpp",
        "lib/bp_support_algorithm.cpp",
        "lib/nn_dict_dynamic.cpp",
    ],
)

PRIVATE_HDRS = [
    "include/sdsl/simple_sds.hpp",
    "include/sdsl/rle_vector.hpp",
    "include/sdsl/bits.hpp",
    "include/sdsl/bit_vector_il.hpp",
    "include/sdsl/bit_vectors.hpp",
    "include/sdsl/coder_comma.hpp",
    "include/sdsl/coder_elias_delta.hpp",
    "include/sdsl/coder_elias_gamma.hpp",
    "include/sdsl/coder_fibonacci.hpp",
    "include/sdsl/coder.hpp",
    "include/sdsl/config.hpp",
    "include/sdsl/construct_bwt.hpp",
    "include/sdsl/construct_config.hpp",
    "include/sdsl/construct.hpp",
    "include/sdsl/construct_isa.hpp",
    "include/sdsl/construct_lcp_helper.hpp",
    "include/sdsl/construct_lcp.hpp",
    "include/sdsl/construct_sa.hpp",
    "include/sdsl/construct_sa_se.hpp",
    "include/sdsl/csa_alphabet_strategy.hpp",
    "include/sdsl/csa_bitcompressed.hpp",
    "include/sdsl/csa_sada.hpp",
    "include/sdsl/csa_sampling_strategy.hpp",
    "include/sdsl/csa_wt.hpp",
    "include/sdsl/dac_vector.hpp",
    "include/sdsl/enc_vector.hpp",
    "include/sdsl/fast_cache.hpp",
    "include/sdsl/hyb_vector.hpp",
    "include/sdsl/int_vector_buffer.hpp",
    "include/sdsl/int_vector.hpp",
    "include/sdsl/inv_perm_support.hpp",
    "include/sdsl/io.hpp",
    "include/sdsl/iterators.hpp",
    "include/sdsl/memory_management.hpp",
    "include/sdsl/qsufsort.hpp",
    "include/sdsl/ram_filebuf.hpp",
    "include/sdsl/ram_fs.hpp",
    "include/sdsl/rank_support.hpp",
    "include/sdsl/rank_support_scan.hpp",
    "include/sdsl/rank_support_v5.hpp",
    "include/sdsl/rank_support_v.hpp",
    "include/sdsl/rrr_helper.hpp",
    "include/sdsl/rrr_vector_15.hpp",
    "include/sdsl/rrr_vector.hpp",
    "include/sdsl/sdsl_concepts.hpp",
    "include/sdsl/sd_vector.hpp",
    "include/sdsl/select_support.hpp",
    "include/sdsl/select_support_mcl.hpp",
    "include/sdsl/select_support_scan.hpp",
    "include/sdsl/sfstream.hpp",
    "include/sdsl/structure_tree.hpp",
    "include/sdsl/suffix_array_algorithm.hpp",
    "include/sdsl/suffix_array_helper.hpp",
    "include/sdsl/suffix_arrays.hpp",
    "include/sdsl/uint128_t.hpp",
    "include/sdsl/uint256_t.hpp",
    "include/sdsl/uintx_t.hpp",
    "include/sdsl/util.hpp",
    "include/sdsl/vectors.hpp",
    "include/sdsl/vlc_vector.hpp",
    "include/sdsl/wavelet_trees.hpp",
    "include/sdsl/wm_int.hpp",
    "include/sdsl/wt_algorithm.hpp",
    "include/sdsl/wt_ap.hpp",
    "include/sdsl/wt_blcd.hpp",
    "include/sdsl/wt_gmr.hpp",
    "include/sdsl/wt_helper.hpp",
    "include/sdsl/wt_huff.hpp",
    "include/sdsl/wt_hutu.hpp",
    "include/sdsl/wt_int.hpp",
    "include/sdsl/wt_pc.hpp",
    "include/sdsl/wt_rlmn.hpp",
]

cc_library(
    name = "sdsl_lite",
    srcs = PRIVATE_HDRS + SRCS,
    hdrs = [
        "include/sdsl/int_vector.hpp",
        "include/sdsl/sd_vector.hpp",
        "include/sdsl/simple_sds.hpp",
        "include/sdsl/suffix_arrays.hpp",
    ],
    copts = LIBRARY_COPTS,
    includes = includes,
    deps = [
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/log:absl_log",
        "@libdivsufsort",
        "@libdivsufsort//:libdivsufsort64",
        # "@zstdlib//:dict-builder-lib",
    ],
)

cc_binary(
    name = "test_fmindex",
    srcs = ["test_fmindex.cc"],
    copts = ["-fexceptions"],
    deps = [":sdsl_lite"],
)
