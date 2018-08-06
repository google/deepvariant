"""Bazel rules for nucleus_py_* targets that can depend on C++ code."""

# A provider with one field, transitive_deps.
CppFilesInfo = provider(fields = ["transitive_deps"])

def get_transitive_deps(deps):
    """Return all the transitive dependencies of deps."""
    return depset(
        deps,
        transitive = [dep[CppFilesInfo].transitive_deps for dep in deps],
    )

def _py_cdeps_impl(ctx):
    return [CppFilesInfo(transitive_deps = get_transitive_deps(ctx.attr.deps))]

# A rule for storing the C libraries a python library or extension depends on.
py_cdeps = rule(
    implementation = _py_cdeps_impl,
    attrs = {
        "deps": attr.label_list(),
    },
)

def _classify_dependencies(deps = []):
    """Return a 3-tuple with the C, Python, & Python extension subsets of deps."""
    c_deps = []
    py_deps = []
    py_ext_deps = []
    for dep in deps:
        kind = ""
        er = native.existing_rule(dep)
        if er:
            kind = er["kind"]
        if kind == "cc_library":
            c_deps.append(dep)
        elif kind == "nucleus_py_extension":
            py_ext_deps.append(dep)
        else:
            py_deps.append(dep)
    return (c_deps, py_deps, py_ext_deps)

def nucleus_py_library(name, deps = [], **kwargs):
    """A py_library that can depend on cc_library's.

    Args:
      name: The name of the py_library target.
      deps: The python and C++ dependencies of the target.
      **kwargs:  Any additional arguments to py_library.
    """
    c_deps, py_deps, py_ext_deps = _classify_dependencies(deps)
    native.py_library(name = name, deps = py_deps, **kwargs)
    py_cdeps(
        name = name + "_cdeps",
        deps = c_deps + py_ext_deps,
    )

def nucleus_py_extension(name, srcs = [], deps = [], **kwargs):
    """Create a C++ library that extends Python.

    Args:
      name: The name of the extension.  It must be the actual module name.
        For example, it should be the foo in "import foo" used to load the
        module.
      srcs: The C++ files implementing the module.
      deps: The C++ libraries the extension depends on.
      **kwargs:  Any additional arguments to cc_binary.
    """
    so_name = name + ".so"
    native.cc_binary(
        name = so_name,
        linkstatic = 0,
        linkshared = 1,
        srcs = srcs,
        deps = ["//external:python_headers"],
        **kwargs
    )

    # Don't put the dependencies into the binary, but instead propagate
    # them for the nucleus_py_binary rules to use later.
    py_cdeps(
        name = name + "_cdeps",
        deps = deps,
    )

def _add_header_impl(ctx):
    header_loc = ctx.outputs.out + "_header"
    out = ctx.outputs.out
    ctx.actions.write(
        output = header_loc,
        content = ctx.attr.header,
        is_executable = True,
    )

    ctx.actions.run_shell(
        inputs = [header_loc, ctx.attrs.src],
        outputs = [out],
        command = "cat $1 $2 > $3",
        arguments = [header_loc, ctx.attrs.src, out],
    )

add_header = rule(
    implementation = _add_header_impl,
    attrs = {
        "src": attr.label(allow_single_file = True),
        "header": attr.string(),
    },
)

def nucleus_py_binary(name, srcs = [], deps = [], data = [], **kwargs):
    """A py_binary whose C++ dependencies are put into a single .so file.

    Args:
      name: The name of the py_binary.
      srcs: The list of Python source files for the binary.
      deps: The python and C++ dependencies of the py_binary.
      data: The data files used by the py_binary.
      **kwargs:  Any additional arguments to py_binary.
    """
    c_deps, py_deps, py_ext_deps = _classify_dependencies(deps)
    trans_deps = get_transitive_deps(c_deps + py_ext_deps)
    extended_data = data[:]
    new_srcs = srcs[:]
    if len(trans_deps) > 0:
        # Create a .so containing all of the C++ dependencies.
        so_name = name + ".so"
        extended_data.append(so_name)
        native.cc_binary(
            name = so_name,
            linkstatic = 0,
            linkshared = 1,
            deps = trans_deps,
        )
        prelude = "import ctypes\nctypes.CDLL(\"" + so_name
        prelude += "\", ctypes.RTLD_GLOBAL)"
        if len(srcs) > 1:
            fail("nucleus_py_binary currently only supports one src")
        new_srcs[0] = "load_then_" + srcs[0]
        add_header(
            name = new_srcs[0],
            src = srcs[0],
            header = prelude,
        )

    native.py_binary(
        name = name,
        srcs = new_srcs,
        data = extended_data,
        deps = py_deps,
        **kwargs
    )
