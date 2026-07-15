// dv_signpost.h — Apple os_signpost wrappers for Instruments profiling.
//
// Wraps `os_signpost_interval_begin/end` and `os_signpost_event_emit` so
// hot paths in deepvariant native code can be instrumented for
// Time Profiler / Points of Interest in Instruments.app without polluting
// Linux/non-Apple builds (the macros become no-ops on non-__APPLE__).
//
// Usage:
//   #include "deepvariant/native/dv_signpost.h"
//   ...
//   DV_SIGNPOST_INTERVAL_BEGIN(MakeExamples, "chr20:10M-10.1M");
//   ... heavy work ...
//   DV_SIGNPOST_INTERVAL_END(MakeExamples);
//
//   DV_SIGNPOST_EVENT(CallVariants, "batch=512");
//
// View in Instruments:
//   xctrace record --template 'Points of Interest' \
//     --launch -- ./build-macos/bin/deepvariant run [args...]
//   open *.trace
// Each DV_SIGNPOST_* call appears in the "Points of Interest" track.
//
// Subsystem identifier "com.demaille.deepvariant" used uniformly so
// Instruments groups all our signposts together.

#pragma once

#if defined(__APPLE__)

#include <os/log.h>
#include <os/signpost.h>

namespace deepvariant {
namespace signpost {

// Singleton log handle. Created on first use; thread-safe via static
// local init (C++11 magic-static guarantees one-time init).
inline os_log_t Logger() {
  static os_log_t log = os_log_create("com.demaille.deepvariant", "perf");
  return log;
}

}  // namespace signpost
}  // namespace deepvariant

// Begin an interval. Use a unique name (becomes the C++ identifier of a
// stack-local os_signpost_id_t variable). Must be paired with END.
#define DV_SIGNPOST_INTERVAL_BEGIN(name, fmt_or_str)                           \
  os_signpost_id_t _dv_sp_##name =                                              \
      os_signpost_id_generate(::deepvariant::signpost::Logger());               \
  os_signpost_interval_begin(::deepvariant::signpost::Logger(),                 \
                              _dv_sp_##name, #name, "%s", fmt_or_str)

// End an interval started with DV_SIGNPOST_INTERVAL_BEGIN(name, ...).
#define DV_SIGNPOST_INTERVAL_END(name)                                          \
  os_signpost_interval_end(::deepvariant::signpost::Logger(),                   \
                            _dv_sp_##name, #name)

// One-shot event marker (no duration). For batch boundaries, queue
// fills, etc.
#define DV_SIGNPOST_EVENT(name, fmt_or_str)                                     \
  os_signpost_event_emit(::deepvariant::signpost::Logger(),                     \
                          OS_SIGNPOST_ID_EXCLUSIVE, #name, "%s", fmt_or_str)

#else  // !__APPLE__

// No-op on non-Apple platforms. Compile to nothing.
#define DV_SIGNPOST_INTERVAL_BEGIN(name, fmt_or_str) ((void)0)
#define DV_SIGNPOST_INTERVAL_END(name)               ((void)0)
#define DV_SIGNPOST_EVENT(name, fmt_or_str)          ((void)0)

#endif
