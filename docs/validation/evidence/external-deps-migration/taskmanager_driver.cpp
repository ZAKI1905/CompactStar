// Test-only caller of the unchanged CompactStar TaskManager API.
#include <chrono>
#include <filesystem>
#include <iostream>
#include "CompactStar/Core/TaskManager.hpp"

int main(int argc, char** argv)
{
    if (argc != 3) return 2;
    const std::string root = argv[1], mode = argv[2];
    if (!std::filesystem::path(root).is_absolute() || root.size() > 100) return 3;
    if (mode != "T1" && mode != "T2") return 4;
    CompactStar::Core::TaskManager task(1);
    task.SetWrkDir(root);
    task.SetVisEOSDir("EOS/DS(CMF)-1_with_crust.eos");
    task.SetDarEOSDir("EOS/Fermi_Gas_0.8mn.eos");
    task.SetChiMass(0.8);
    task.SetGrid({{8e14, 5e15}, 19, "Log"}, {{5e13, 5e16}, 19, "Log"});
    const auto start = std::chrono::steady_clock::now();
    if (mode == "T2") { task.FindDarkEOS(0.8); task.Work(); }
    const auto work = std::chrono::steady_clock::now();
    task.ImportSequence("NStar/Dark_Core/0.8/0.8_19x19_Sequence.tsv");
    task.FindCriticalCurve();
    task.FindMtotContour(2.01);
    const auto contours = std::chrono::steady_clock::now();
    if (mode == "T2") { task.Precision_Task(2.01); task.FindLimits(2.01); }
    const auto finish = std::chrono::steady_clock::now();
    std::cout << "Z201_TIMING " << mode << " work="
              << std::chrono::duration<double>(work-start).count()
              << " contours=" << std::chrono::duration<double>(contours-work).count()
              << " downstream=" << std::chrono::duration<double>(finish-contours).count()
              << " total=" << std::chrono::duration<double>(finish-start).count() << '\n';
}
