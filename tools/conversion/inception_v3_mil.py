"""DeepVariant Inception-v3 expressed in coremltools MIL.

Builds a coremltools MIL program from weights read directly by
TensorBundle (no TF, no PyTorch, no ONNX).

Architecture: standard Keras InceptionV3, first Conv2D adapted from
3-channel to 7-channel pileup input (kernel shape [3,3,7,32]).

Weight key convention:
    layer_with_weights-N/kernel/.ATTRIBUTES/VARIABLE_VALUE  → Conv (HWIO)
    layer_with_weights-N/beta/.ATTRIBUTES/VARIABLE_VALUE    → BN β
    layer_with_weights-N/moving_mean/.ATTRIBUTES/VARIABLE_VALUE   → BN μ
    layer_with_weights-N/moving_variance/.ATTRIBUTES/VARIABLE_VALUE → BN σ²
    (gamma is NOT stored — frozen at 1.0, we use np.ones)

BN indices are NOT always conv+1 inside Inception blocks (branches are
interleaved in the Keras TrackableObjectGraph). All (conv_n, bn_n) pairs
below are verified against the real WGS SavedModel bundle dump.

Input spec:  (N, 100, 221, 7)  NHWC — matches the original SavedModel.
             A mb.transpose at the top converts to NCHW before the convs.
Output spec: (N, 3)            softmax probabilities.

Usage:
    from tensor_bundle_reader import TensorBundle
    import inception_v3_mil as iv3
    bundle = TensorBundle("models/wgs/variables/variables")
    prog = iv3.build_program(bundle)
    # Then: ct.convert(prog, ...) in convert_coreml.py
"""

from __future__ import annotations

import numpy as np

try:
    from coremltools.converters.mil import Builder as mb
except ImportError as exc:
    raise ImportError(
        "coremltools not installed — activate venv-coreml first: "
        "source tools/conversion/venv-coreml/bin/activate"
    ) from exc

_ATTR = ".ATTRIBUTES/VARIABLE_VALUE"


# ---------------------------------------------------------------------------
# Weight helpers
# ---------------------------------------------------------------------------

def _k(bundle, n: int) -> np.ndarray:
    """Conv kernel for layer_with_weights-N in Core ML OIHW layout."""
    raw = bundle.read_tensor(f"layer_with_weights-{n}/kernel/{_ATTR}")
    return np.transpose(raw, (3, 2, 0, 1)).astype(np.float32)  # HWIO→OIHW


# Expected stem-conv kernel geometry in OIHW: 7-channel pileup input, 3×3
# spatial. A regression that flips the HWIO→OIHW perm in _k() would build a
# valid-but-wrong model, so we assert this once at build time.
_STEM_IN_CHANNELS = 7
_STEM_SPATIAL = (3, 3)

# Number of output classes for the WGS model (hom-ref / het / hom-alt).
_NUM_CLASSES = 3


def _bn_params(bundle, n: int) -> tuple:
    """Return (gamma, beta, mean, var) for layer_with_weights-N.

    gamma is not stored in the DV checkpoint (frozen=1.0), so we supply
    a vector of ones with the same shape as beta.
    """
    beta = bundle.read_tensor(
        f"layer_with_weights-{n}/beta/{_ATTR}"
    ).astype(np.float32)
    mean = bundle.read_tensor(
        f"layer_with_weights-{n}/moving_mean/{_ATTR}"
    ).astype(np.float32)
    var = bundle.read_tensor(
        f"layer_with_weights-{n}/moving_variance/{_ATTR}"
    ).astype(np.float32)
    return np.ones_like(beta), beta, mean, var


def _cbr(
    x,
    bundle,
    conv_n: int,
    bn_n: int,
    strides: list[int] | None = None,
    padding: str = "same",
    name: str = "",
) -> object:
    """Conv2D (OIHW kernel) + BatchNorm + ReLU."""
    strides = strides or [1, 1]
    kernel = _k(bundle, conv_n)
    gamma, beta, mean, var = _bn_params(bundle, bn_n)
    x = mb.conv(
        x=x, weight=kernel, strides=strides, pad_type=padding,
        name=f"{name}_c",
    )
    x = mb.batch_norm(
        x=x, mean=mean, variance=var, gamma=gamma, beta=beta,
        # Keras BatchNormalization default epsilon is 1e-3 (NOT 1e-4).
        # Inception-v3 SavedModels are trained with epsilon=1e-3.
        # See metal_inference.mm::kBNEpsilon and CLAUDE.md "Pitfalls".
        # Wrong epsilon → subtle scale mismatch on channels with small
        # variance → INDEL F1 collapse.
        epsilon=1e-3, name=f"{name}_bn",
    )
    return mb.relu(x=x, name=f"{name}_r")


def _avg_cbr(x, bundle, conv_n: int, bn_n: int, name: str) -> object:
    """AvgPool 3×3 (same, no count boundary) then CBR."""
    x = mb.avg_pool(
        x=x, kernel_sizes=[3, 3], strides=[1, 1], pad_type="same",
        exclude_padding_from_average=True, name=f"{name}_ap",
    )
    return _cbr(x, bundle, conv_n, bn_n, name=name)


# ---------------------------------------------------------------------------
# InceptionA blocks  (layers 10-51)
# Input 192 → 256 → 288 → 288
# Branches: 1×1 | 1×1→5×5 | 1×1→3×3→3×3 | AvgPool→1×1
# ---------------------------------------------------------------------------

def _mixed_5b(x, bundle) -> object:
    """Mixed_5b: layers 10-23, 192→256.

    BUG FIX (2026-05-24): b1 ↔ b3_3a pairs were swapped — wrong
    layer-with-weights indices for the branch1x1 conv vs the
    branch3x3dbl reduce conv. Authoritative pairs from
    metal_inference.mm:Mixed_5b (Phase 5.5a 2026-04-28 fix).
    """
    b1 = _cbr(x, bundle, 16, 20, name="5b_1")    # was (10,11)
    b5 = _cbr(x, bundle, 12, 14, name="5b_5a")
    b5 = _cbr(b5, bundle, 17, 21, [1, 1], "same", "5b_5b")
    b3 = _cbr(x, bundle, 10, 11, name="5b_3a")   # was (16,20)
    b3 = _cbr(b3, bundle, 13, 15, name="5b_3b")
    b3 = _cbr(b3, bundle, 18, 22, name="5b_3c")
    bp = _avg_cbr(x, bundle, 19, 23, "5b_p")
    return mb.concat(values=(b1, b5, b3, bp), axis=1, name="5b")


def _mixed_5c(x, bundle) -> object:
    """Mixed_5c: layers 24-37, 256→288. Same swap as Mixed_5b."""
    b1 = _cbr(x, bundle, 30, 34, name="5c_1")    # was (24,25)
    b5 = _cbr(x, bundle, 26, 28, name="5c_5a")
    b5 = _cbr(b5, bundle, 31, 35, [1, 1], "same", "5c_5b")
    b3 = _cbr(x, bundle, 24, 25, name="5c_3a")   # was (30,34)
    b3 = _cbr(b3, bundle, 27, 29, name="5c_3b")
    b3 = _cbr(b3, bundle, 32, 36, name="5c_3c")
    bp = _avg_cbr(x, bundle, 33, 37, "5c_p")
    return mb.concat(values=(b1, b5, b3, bp), axis=1, name="5c")


def _mixed_5d(x, bundle) -> object:
    """Mixed_5d: layers 38-51, 288→288. Same swap as Mixed_5b."""
    b1 = _cbr(x, bundle, 44, 48, name="5d_1")    # was (38,39)
    b5 = _cbr(x, bundle, 40, 42, name="5d_5a")
    b5 = _cbr(b5, bundle, 45, 49, [1, 1], "same", "5d_5b")
    b3 = _cbr(x, bundle, 38, 39, name="5d_3a")   # was (44,48)
    b3 = _cbr(b3, bundle, 41, 43, name="5d_3b")
    b3 = _cbr(b3, bundle, 46, 50, name="5d_3c")
    bp = _avg_cbr(x, bundle, 47, 51, "5d_p")
    return mb.concat(values=(b1, b5, b3, bp), axis=1, name="5d")


# ---------------------------------------------------------------------------
# Reduction-A  (layers 52-59)
# 288 → 768
# ---------------------------------------------------------------------------

def _mixed_6a(x, bundle) -> object:
    """Mixed_6a (Reduction-A): layers 52-59, 288→768."""
    b3 = _cbr(x, bundle, 56, 58, [2, 2], "valid", "6a_3")
    bd = _cbr(x, bundle, 52, 53, name="6a_da")
    bd = _cbr(bd, bundle, 54, 55, name="6a_db")
    bd = _cbr(bd, bundle, 57, 59, [2, 2], "valid", "6a_dc")
    bp = mb.max_pool(
        x=x, kernel_sizes=[3, 3], strides=[2, 2],
        pad_type="valid", name="6a_mp",
    )
    return mb.concat(values=(b3, bd, bp), axis=1, name="6a")


# ---------------------------------------------------------------------------
# InceptionB blocks  (layers 60-139)
# 768 → 768 (factorized 7×1 + 1×7)
# Branches: 1×1 | 1×1→1×7→7×1 | 1×1→7×1→1×7→7×1→1×7 | AvgPool→1×1
# ---------------------------------------------------------------------------

def _mixed_6b(x, bundle) -> object:
    """Mixed_6b: layers 60-79, 768→768 (128-ch factorized).

    BUG FIX (2026-05-24): b7a_b ↔ b7b_c pairs were swapped — wrong
    layer-with-weights indices for the b7a 1×7 conv vs b7b 1×7 conv.
    Authoritative pairs from metal_inference.mm:Mixed_6b.
    """
    b1 = _cbr(x, bundle, 72, 76, name="6b_1")
    b7a = _cbr(x, bundle, 64, 66, name="6b_7aa")
    b7a = _cbr(b7a, bundle, 68, 70, name="6b_7ab")   # was (65,67)
    b7a = _cbr(b7a, bundle, 73, 77, name="6b_7ac")
    b7b = _cbr(x, bundle, 60, 61, name="6b_7ba")
    b7b = _cbr(b7b, bundle, 62, 63, name="6b_7bb")
    b7b = _cbr(b7b, bundle, 65, 67, name="6b_7bc")   # was (68,70)
    b7b = _cbr(b7b, bundle, 69, 71, name="6b_7bd")
    b7b = _cbr(b7b, bundle, 74, 78, name="6b_7be")
    bp = _avg_cbr(x, bundle, 75, 79, "6b_p")
    return mb.concat(values=(b1, b7a, b7b, bp), axis=1, name="6b")


def _mixed_6c(x, bundle) -> object:
    """Mixed_6c: layers 80-99, 768→768 (160-ch factorized). Same swap as 6b."""
    b1 = _cbr(x, bundle, 92, 96, name="6c_1")
    b7a = _cbr(x, bundle, 84, 86, name="6c_7aa")
    b7a = _cbr(b7a, bundle, 88, 90, name="6c_7ab")   # was (85,87)
    b7a = _cbr(b7a, bundle, 93, 97, name="6c_7ac")
    b7b = _cbr(x, bundle, 80, 81, name="6c_7ba")
    b7b = _cbr(b7b, bundle, 82, 83, name="6c_7bb")
    b7b = _cbr(b7b, bundle, 85, 87, name="6c_7bc")   # was (88,90)
    b7b = _cbr(b7b, bundle, 89, 91, name="6c_7bd")
    b7b = _cbr(b7b, bundle, 94, 98, name="6c_7be")
    bp = _avg_cbr(x, bundle, 95, 99, "6c_p")
    return mb.concat(values=(b1, b7a, b7b, bp), axis=1, name="6c")


def _mixed_6d(x, bundle) -> object:
    """Mixed_6d: layers 100-119, 768→768 (160-ch factorized). Same swap as 6b."""
    b1 = _cbr(x, bundle, 112, 116, name="6d_1")
    b7a = _cbr(x, bundle, 104, 106, name="6d_7aa")
    b7a = _cbr(b7a, bundle, 108, 110, name="6d_7ab")  # was (105,107)
    b7a = _cbr(b7a, bundle, 113, 117, name="6d_7ac")
    b7b = _cbr(x, bundle, 100, 101, name="6d_7ba")
    b7b = _cbr(b7b, bundle, 102, 103, name="6d_7bb")
    b7b = _cbr(b7b, bundle, 105, 107, name="6d_7bc")  # was (108,110)
    b7b = _cbr(b7b, bundle, 109, 111, name="6d_7bd")
    b7b = _cbr(b7b, bundle, 114, 118, name="6d_7be")
    bp = _avg_cbr(x, bundle, 115, 119, "6d_p")
    return mb.concat(values=(b1, b7a, b7b, bp), axis=1, name="6d")


def _mixed_6e(x, bundle) -> object:
    """Mixed_6e: layers 120-139, 768→768 (192-ch factorized). Same swap as 6b."""
    b1 = _cbr(x, bundle, 132, 136, name="6e_1")
    b7a = _cbr(x, bundle, 124, 126, name="6e_7aa")
    b7a = _cbr(b7a, bundle, 128, 130, name="6e_7ab")  # was (125,127)
    b7a = _cbr(b7a, bundle, 133, 137, name="6e_7ac")
    b7b = _cbr(x, bundle, 120, 121, name="6e_7ba")
    b7b = _cbr(b7b, bundle, 122, 123, name="6e_7bb")
    b7b = _cbr(b7b, bundle, 125, 127, name="6e_7bc")  # was (128,130)
    b7b = _cbr(b7b, bundle, 129, 131, name="6e_7bd")
    b7b = _cbr(b7b, bundle, 134, 138, name="6e_7be")
    bp = _avg_cbr(x, bundle, 135, 139, "6e_p")
    return mb.concat(values=(b1, b7a, b7b, bp), axis=1, name="6e")


# ---------------------------------------------------------------------------
# Reduction-B  (layers 140-151)
# 768 → 1280
# ---------------------------------------------------------------------------

def _mixed_7a(x, bundle) -> object:
    """Mixed_7a (Reduction-B): layers 140-151, 768→1280.

    BUG FIX (2026-05-24): b3_a ↔ b7_a swapped — wrong indices for the
    branch3x3 reduce conv vs branch7x7 reduce conv. Authoritative
    pairs from metal_inference.mm:Mixed_7a.
    """
    b3 = _cbr(x, bundle, 144, 146, name="7a_3a")   # was (140,141)
    b3 = _cbr(b3, bundle, 148, 150, [2, 2], "valid", "7a_3b")
    b7 = _cbr(x, bundle, 140, 141, name="7a_7a")   # was (144,146)
    b7 = _cbr(b7, bundle, 142, 143, name="7a_7b")
    b7 = _cbr(b7, bundle, 145, 147, name="7a_7c")
    b7 = _cbr(b7, bundle, 149, 151, [2, 2], "valid", "7a_7d")
    bp = mb.max_pool(
        x=x, kernel_sizes=[3, 3], strides=[2, 2],
        pad_type="valid", name="7a_mp",
    )
    return mb.concat(values=(b3, b7, bp), axis=1, name="7a")


# ---------------------------------------------------------------------------
# InceptionC blocks  (layers 152-187)
# 1280→2048→2048
# Branches: 1×1 | 1×1→{1×3,3×1} | 1×1→3×3→{1×3,3×1} | AvgPool→1×1
# ---------------------------------------------------------------------------

def _inception_c(x, bundle, idx, name: str) -> object:
    """Generic InceptionC block.

    idx = (b1_c, b1_bn,
           b3a_c1, b3a_bn1, b3a_1x3_c, b3a_1x3_bn, b3a_3x1_c, b3a_3x1_bn,
           b3b_c1, b3b_bn1, b3b_c2, b3b_bn2,
               b3b_1x3_c, b3b_1x3_bn, b3b_3x1_c, b3b_3x1_bn,
           bp_c, bp_bn)
    """
    (b1c, b1n, b3ac1, b3an1, b3a1c, b3a1n, b3a3c, b3a3n,
     b3bc1, b3bn1, b3bc2, b3bn2, b3b1c, b3b1n, b3b3c, b3b3n,
     bpc, bpn) = idx

    b1 = _cbr(x, bundle, b1c, b1n, name=f"{name}_1")
    b3a = _cbr(x, bundle, b3ac1, b3an1, name=f"{name}_3aa")
    b3a_1x3 = _cbr(b3a, bundle, b3a1c, b3a1n, name=f"{name}_3a1x3")
    b3a_3x1 = _cbr(b3a, bundle, b3a3c, b3a3n, name=f"{name}_3a3x1")
    b3b = _cbr(x, bundle, b3bc1, b3bn1, name=f"{name}_3ba")
    b3b = _cbr(b3b, bundle, b3bc2, b3bn2, name=f"{name}_3bb")
    b3b_1x3 = _cbr(b3b, bundle, b3b1c, b3b1n, name=f"{name}_3b1x3")
    b3b_3x1 = _cbr(b3b, bundle, b3b3c, b3b3n, name=f"{name}_3b3x1")
    bp = _avg_cbr(x, bundle, bpc, bpn, f"{name}_p")
    return mb.concat(
        values=(b1, b3a_1x3, b3a_3x1, b3b_1x3, b3b_3x1, bp),
        axis=1, name=name,
    )


def _mixed_7b(x, bundle) -> object:
    """Mixed_7b: layers 152-169, 1280→2048."""
    return _inception_c(x, bundle, (
        162, 168,
        154, 156, 158, 163, 159, 164,
        152, 153, 155, 157, 160, 165, 161, 166,
        167, 169,
    ), "7b")


def _mixed_7c(x, bundle) -> object:
    """Mixed_7c: layers 170-187, 2048→2048."""
    return _inception_c(x, bundle, (
        180, 186,
        172, 174, 176, 181, 177, 182,
        170, 171, 173, 175, 178, 183, 179, 184,
        185, 187,
    ), "7c")


# ---------------------------------------------------------------------------
# Full MIL program
# ---------------------------------------------------------------------------

def build_program(
    bundle, batch_min: int = 1, batch_max: int = 4096,
) -> object:
    """Return the coremltools MIL program for the DeepVariant WGS model.

    Input:  x  shape (N, 100, 221, 7)  NHWC float32 — matches original.
    Output: classification shape (N, 3)  float32 softmax.
    """
    from coremltools.converters.mil.mil.program import get_new_symbol

    # Validate the stem conv geometry up front: OIHW must be (O, 7, 3, 3) for
    # the 7-channel pileup input. Guards against a perm regression in _k().
    stem_kernel = _k(bundle, 0)
    if stem_kernel.ndim != 4 or stem_kernel.shape[1] != _STEM_IN_CHANNELS or (
        stem_kernel.shape[2], stem_kernel.shape[3]
    ) != _STEM_SPATIAL:
        raise ValueError(
            "stem conv kernel 'layer_with_weights-0/kernel' has OIHW shape "
            f"{stem_kernel.shape}, expected "
            f"(O, {_STEM_IN_CHANNELS}, {_STEM_SPATIAL[0]}, {_STEM_SPATIAL[1]}) — "
            "the HWIO→OIHW transpose in _k() may be wrong"
        )

    # Use a symbol for the batch dim so ct.convert() can override it with
    # a RangeDim (see convert_coreml.py inputs= parameter).
    batch_sym = get_new_symbol()
    @mb.program(
        input_specs=[mb.TensorSpec(shape=(batch_sym, 100, 221, 7))]
    )
    def prog(x):
        # NHWC (N,100,221,7) → NCHW (N,7,100,221) for Core ML convs.
        x = mb.transpose(x=x, perm=[0, 3, 1, 2], name="nhwc2nchw")

        # Stem
        x = _cbr(x, bundle, 0, 1, [2, 2], "valid", "s1a")
        x = _cbr(x, bundle, 2, 3, [1, 1], "valid", "s2a")
        x = _cbr(x, bundle, 4, 5, [1, 1], "same",  "s2b")
        x = mb.max_pool(
            x=x, kernel_sizes=[3, 3], strides=[2, 2],
            pad_type="valid", name="mp3a",
        )
        x = _cbr(x, bundle, 6, 7, [1, 1], "valid", "s3b")
        x = _cbr(x, bundle, 8, 9, [1, 1], "valid", "s4a")
        x = mb.max_pool(
            x=x, kernel_sizes=[3, 3], strides=[2, 2],
            pad_type="valid", name="mp5a",
        )

        # InceptionA
        x = _mixed_5b(x, bundle)
        x = _mixed_5c(x, bundle)
        x = _mixed_5d(x, bundle)

        # Reduction-A
        x = _mixed_6a(x, bundle)

        # InceptionB (factorized)
        x = _mixed_6b(x, bundle)
        x = _mixed_6c(x, bundle)
        x = _mixed_6d(x, bundle)
        x = _mixed_6e(x, bundle)

        # Reduction-B
        x = _mixed_7a(x, bundle)

        # InceptionC
        x = _mixed_7b(x, bundle)
        x = _mixed_7c(x, bundle)

        # Global avg pool → (N, 2048, 1, 1) → (N, 2048)
        x = mb.reduce_mean(x=x, axes=[2, 3], keep_dims=True, name="gap")
        x = mb.squeeze(x=x, axes=[2, 3], name="squeeze")

        # Dense 2048 → 3 (kernel + bias — the bias is the only one in the
        # whole model since every conv is fused with BN).
        w = bundle.read_tensor(
            f"layer_with_weights-188/kernel/{_ATTR}"
        ).T.astype(np.float32)
        b = bundle.read_tensor(
            f"layer_with_weights-188/bias/{_ATTR}"
        ).astype(np.float32)
        # mb.linear weight is [Dout, Din]; the Dense kernel is stored [Din, Dout]
        # so the .T above is load-bearing. Assert the final layer is a 2-D
        # weight whose output dim matches the bias and the expected class count
        # — a dropped .T would surface here rather than as silent miscalibration.
        if w.ndim != 2:
            raise ValueError(
                "final Dense weight 'layer_with_weights-188/kernel' must be 2-D "
                f"after transpose, got shape {w.shape}"
            )
        if b.ndim != 1 or b.shape[0] != w.shape[0]:
            raise ValueError(
                "final Dense bias 'layer_with_weights-188/bias' shape "
                f"{b.shape} is inconsistent with weight output dim {w.shape[0]}"
            )
        if w.shape[0] != _NUM_CLASSES:
            raise ValueError(
                "final Dense weight 'layer_with_weights-188/kernel' has output "
                f"dim {w.shape[0]}, expected {_NUM_CLASSES} classes — the "
                ".T transpose may have been dropped"
            )
        x = mb.linear(x=x, weight=w, bias=b, name="logits")

        return mb.softmax(x=x, axis=1, name="classification")

    return prog
