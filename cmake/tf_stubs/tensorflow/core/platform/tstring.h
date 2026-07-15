// tensorflow::tstring is just std::string in our TF-free build.
#pragma once
#include <string>
namespace tensorflow { using tstring = std::string; }
