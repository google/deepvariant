#pragma once
namespace deepvariant {

// Subcommand entry points — each parses its own flags from argv.
int RunMakeExamples(int argc, char** argv);
int RunCallVariants(int argc, char** argv);
int RunPostprocessVariants(int argc, char** argv);

// "run" subcommand: chains all three stages in process.
int RunAll(int argc, char** argv);

}  // namespace deepvariant
