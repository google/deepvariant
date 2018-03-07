"""Bazel rule for making a zip file."""

def zip_dir(name, srcs, zipname, **kwargs):
  """Zips up an entire directory or Fileset.

  Args:
    name: The name of the target
    srcs: A single-item list with a directory or fileset
    zipname: The name of the output zip file
    **kwargs: Further generic arguments to pass to genrule, e.g. visibility.
  """
  if len(srcs) > 1:
    fail("More than one directory is not supported by zip_dir yet", attr=srcs)
  native.genrule(
      name=name,
      srcs=srcs,
      outs=[zipname],
      cmd="zip $(OUTS) $(SRCS)",
      **kwargs)
