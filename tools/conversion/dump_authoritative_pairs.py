"""Dump authoritative (conv_n, bn_n) pairs for every Keras conv2d_M /
batch_normalization_M in a DeepVariant Inception-v3 SavedModel.

For each Keras `conv2d_M` op in the frozen graph, byte-matches its
kernel const against the bundle's `layer_with_weights-K` entries to
recover the canonical (M → K) mapping. Same for `batch_normalization_M`
beta. Writes one line per matched conv to stdout:

    M  conv_n  bn_n  shape

Run inside the conversion Docker:

    docker run --rm --platform linux/amd64 \\
      -v $(realpath models/wgs):/in:ro \\
      -v $(realpath tools/conversion):/work:ro \\
      google/deepvariant:1.10.0 \\
      python3 /work/dump_authoritative_pairs.py

Output is consumed by `generate_inception_blocks.py` (TBD) which
emits the corresponding C++ Mixed_* functions in metal_inference.mm.
"""
import numpy as np
import tensorflow as tf
from tensorflow.python.framework.convert_to_constants import (
    convert_variables_to_constants_v2,
)

sm = tf.saved_model.load('/in')
fn = sm.signatures['serving_default']
frozen = convert_variables_to_constants_v2(fn)
gd = frozen.graph.as_graph_def()

# Bundle entries
ckpt = tf.train.load_checkpoint('/in/variables/variables')
bundle_kernel = {}  # n -> (shape, first4)
bundle_beta = {}
for n in range(0, 200):
    kname = f'layer_with_weights-{n}/kernel/.ATTRIBUTES/VARIABLE_VALUE'
    bname = f'layer_with_weights-{n}/beta/.ATTRIBUTES/VARIABLE_VALUE'
    if kname in ckpt.get_variable_to_shape_map():
        a = ckpt.get_tensor(kname)
        bundle_kernel[n] = (a.shape, tuple(a.flatten()[:4].tolist()))
    if bname in ckpt.get_variable_to_shape_map():
        a = ckpt.get_tensor(bname)
        bundle_beta[n] = (a.shape, tuple(a.flatten()[:4].tolist()))

def get_value(node_name):
    for nn in gd.node:
        if nn.name == node_name:
            if nn.op == 'Const':
                t = nn.attr['value'].tensor
                if not t.tensor_content:
                    return None
                shape = tuple(d.size for d in t.tensor_shape.dim)
                return np.frombuffer(t.tensor_content, dtype=np.float32).reshape(shape)
            elif nn.op == 'Identity':
                return get_value(nn.input[0].split(':')[0])
    return None

# Find max M (number of conv2d / bn ops)
all_convs = sorted({n.name for n in gd.node if 'inceptionv3/conv2d' in n.name and n.op == 'Conv2D'})
all_bns = sorted({n.name for n in gd.node if 'inceptionv3/batch_normalization' in n.name and n.op == 'FusedBatchNormV3'})
print(f'-- found {len(all_convs)} conv ops, {len(all_bns)} bn ops --')

def match_kernel(arr):
    sh = arr.shape
    f4 = tuple(arr.flatten()[:4].tolist())
    for n, (s, t) in bundle_kernel.items():
        if s == sh and all(abs(a - b) < 1e-9 for a, b in zip(t, f4)):
            return n
    return None

def match_beta(arr):
    sh = arr.shape
    f4 = tuple(arr.flatten()[:4].tolist())
    for n, (s, t) in bundle_beta.items():
        if s == sh and all(abs(a - b) < 1e-9 for a, b in zip(t, f4)):
            return n
    return None

# Build mapping conv2d_M -> (conv_n, bn_n)
pairs = []  # (conv2d_idx, conv_n_in_bundle, bn_n_in_bundle, kernel_shape)
for M in range(0, 100):
    target_c = f'StatefulPartitionedCall/inceptionv3/conv2d{"" if M==0 else f"_{M}"}/Conv2D'
    target_b = f'StatefulPartitionedCall/inceptionv3/batch_normalization{"" if M==0 else f"_{M}"}/FusedBatchNormV3'
    cn, bn, kshape = None, None, None
    for nn in gd.node:
        if nn.name == target_c:
            karr = get_value(nn.input[1].split(':')[0])
            if karr is not None:
                cn = match_kernel(karr)
                kshape = karr.shape
        elif nn.name == target_b:
            barr = get_value(nn.input[2].split(':')[0])
            if barr is not None:
                bn = match_beta(barr)
    if cn is not None and bn is not None:
        pairs.append((M, cn, bn, kshape))
    elif cn is not None or bn is not None:
        print(f'  PARTIAL: conv2d_{M} cn={cn} bn={bn}')

print(f'-- complete pairs: {len(pairs)} --')
print('M\tconv_n\tbn_n\tshape')
for M, cn, bn, sh in pairs:
    print(f'{M}\t{cn}\t{bn}\t{sh}')
