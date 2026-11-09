#include "../rotochemical/fixture.hpp"
#include "stellar_archive.hpp"
int main(int argc,char** argv) {
 try {
    require(argc==3,"profile fresh-output");
    const std::filesystem::path profile=argv[1],work=argv[2];
    require(!std::filesystem::exists(work),"refuse repeat sequence: fresh output required");
    std::filesystem::create_directories(work);
    auto src=source(profile/"freegas.tsv");
    require(RC::FrozenSource(src->table_path).sha256=="f6b9c4068a430cc356b54e39c21a0ca17e17b8c2b35d17b836638298d1112d7f","EOS authentication failed");
    NumberSequenceRecipe recipe{src,species,whole,1.10e15,1.095e15,1.105e15,
      "Structure-1 monotone midpoint smooth central branch",{.001,.0005,.00025},80000,
      (work/"sequence").string(),tail_bound};
    recipe.tail_policy_identity="positive-source pe comparison inequalities";recipe.tail_policy_revision="PB13-comparison-v1";
    auto sequence=EquilibriumSequenceNumberDerivative::Compute(recipe);sequence.RequireCurrent();
    // Persist the entire authorized producer result before any downstream assembly.
    const auto path=(work/"sequence.archive").string();Core::StellarArchive::Save(path,sequence);
    const auto hash=RC::FrozenSource(path).sha256;
    {std::ofstream out(work/"SHA256SUMS");out<<hash<<"  sequence.archive\n";}
    auto loaded=Core::StellarArchive::Load(path,hash);
    Core::StellarArchive::Save((work/"roundtrip.archive").string(),*loaded);
    require(RC::FrozenSource(work/"roundtrip.archive").sha256==hash,"archive roundtrip changed numerical/provenance bytes");
    std::cout<<"PASS fifteen-star capture and lossless reconstruction "<<hash<<'\n';
    return 0;
 } catch(const std::exception& e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}
}
