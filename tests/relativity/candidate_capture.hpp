#ifndef COMPACTSTAR_TESTS_RELATIVITY_CANDIDATE_CAPTURE_HPP
#define COMPACTSTAR_TESTS_RELATIVITY_CANDIDATE_CAPTURE_HPP
// Optional test-only evidence sink. It records the actual output of a baseline
// comparison without changing which reference is compared or its acceptance.
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
namespace unit_candidate_evidence {
template<class Writer> void Capture(const char* filename, Writer writer) {
 const char* root=std::getenv("COMPACTSTAR_UNIT1_CAPTURE_DIR");if(!root||!*root)return;
 namespace fs=std::filesystem;const auto dir=fs::canonical(root);
 if(!fs::is_directory(dir))throw std::runtime_error("candidate capture requires an existing directory");
 auto prev=dir.end();
 for(auto it=dir.begin();it!=dir.end();++it){if(prev!=dir.end()&&*prev=="tests"&&*it=="baselines")throw std::runtime_error("candidate capture cannot write governed baselines");prev=it;}
 const auto path=dir/filename;
 if(fs::exists(path))throw std::runtime_error("candidate capture requires a fresh file");
 writer(path);
 if(!fs::is_regular_file(path)||fs::file_size(path)==0)throw std::runtime_error("candidate capture writer failed");
}
}
#endif
