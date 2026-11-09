#pragma once
#include "../bnv/breakpoint_fixture.hpp"
#include "fixture.hpp"
BreakpointTest::Pair NativePair(const Inputs& inputs, const std::shared_ptr<EquilibriumSequenceNumberDerivative>& sequence)
{
    const auto& card=Campaign::RunCards().back();
    require(card.identity=="CPL-P2-LINEAR-QSS-v1"&&card.partition=="P2"&&card.fractional_drive_per_year==-1.0e-12
      &&card.duration_year==5.0e5&&card.checkpoints==8193&&!card.reaction_free,"bounded P2 card changed");
    auto fixture=NativeFixture(inputs.profile,inputs.certificate.string(),inputs.work/"owning-star",80000,sequence);
    auto channels=Channels(fixture);
    for(int i=0;i<2;++i)for(int j=0;j<2;++j)
      require(fixture.z->Values()[i][j]==Campaign::Column(inputs.coefficients,"Z"+std::to_string(i)+std::to_string(j)),"native Z coefficient mismatch");
    require(channels->Entry(RC::UrcaProcess::Me).luminosity_erg_s_Kq==Campaign::Column(inputs.coefficients,"Ltilde_Me"),"native electron thermal coefficient mismatch");
    require(channels->Entry(RC::UrcaProcess::Mmu).luminosity_erg_s_Kq==Campaign::Column(inputs.coefficients,"Ltilde_Mmu"),"native muon thermal coefficient mismatch");
    EOS::CompOSE_Thermo::Options options;options.Tmin_for_derivative_MeV=0;options.clamp_to_domain=false;
    auto thermal=std::make_shared<const RC::FrozenThermalSource>(inputs.thermal,options,
      "controlled mathematical fixed-background free-gas entropy; qualified radial80000");
    auto src=source(inputs.profile/"freegas.tsv");
    auto numbers=ParticleNumbers::Compute(input(*fixture.central,src));numbers.RequireCurrent();
    const double B0=numbers.Values()[0]+numbers.Values()[1];
    require(B0==Campaign::Column(inputs.coefficients,"B_solved"),"native owning B0 mismatch");
    auto tangent=std::make_shared<const BNV::EquilibriumBaryonTangent>(BNV::EquilibriumBaryonTangent::Compute(
      sequence,B0,"Structure-1 rho_c=1.10e15 radial80000 EOS8192",BNV::EquilibriumBaryonTangent::SerializeDomainIdentity(whole)));
    const auto tangent_table=inputs.coefficients.parent_path()/"tangent.tsv";
    const std::array<std::string,3> axes{{"n","e","mu"}};
    for(int i=0;i<3;++i) {
      require(tangent->ClosedValues()[i]==Campaign::Column(tangent_table,"t_"+axes[i]),"native tangent mismatch");
      require(tangent->NumericalErrors()[i]==Campaign::Column(tangent_table,"u_t_"+axes[i]),"native tangent error mismatch");
    }
    require(tangent->RawClosureResidual()==Campaign::Column(tangent_table,"raw_closure")&&tangent->ClosureBudget()==Campaign::Column(tangent_table,"closure_budget"),"native tangent closure mismatch");
    auto bnv_token=std::make_shared<BNV::BnvDependencyToken>();
    auto monitor=Phase6A1Test::LoadFrozenMonitor(inputs.frozen,tangent,bnv_token);
    auto run_token=std::make_shared<RC::RunDependencyToken>();
    auto spin=std::make_shared<const BNV::StaticZeroSpinHistory>(run_token);
    auto q=Campaign::Qualification(inputs.profile,inputs.certificate,inputs.entry);
    q.background_inputs=RC::QualifiedBackgroundInputs::NativeStructure1;
    auto ordinary=Context(fixture,channels,thermal,spin,q,run_token);
    auto historical=q;historical.background_inputs=RC::QualifiedBackgroundInputs::HistoricalStructure1;
    MustRefuse([&]{Context(fixture,channels,thermal,spin,historical,run_token);},"historical Mac binding cannot substitute native owner");
    auto missing=q;missing.certificate_path=(inputs.work/"missing-certificate").string();
    MustRefuse([&]{Context(fixture,channels,thermal,spin,missing,run_token);},"missing chemical certificate");
    auto foreign=q;foreign.certificate_path=inputs.entry.string();
    MustRefuse([&]{Context(fixture,channels,thermal,spin,foreign,run_token);},"foreign chemical certificate");
    MustRefuse([&]{Phase6A1Test::LoadFrozenMonitor(inputs.work/"missing-frozen",tangent,bnv_token);},"missing frozen certificate");
    auto bad_rows=Phase6A1Test::ReadFrozenCertificateRows(inputs.frozen);bad_rows[0].B_solved_count*=1.001;
    MustRefuse([&]{BNV::FrozenSensitivityCertificate bad(B0,tangent->StarIdentity(),tangent->DomainIdentity(),"mismatch test",bnv_token,bad_rows);},"mismatched frozen owner");
    const double mu_B=Campaign::Column(inputs.coefficients,"mu_B_inf");
    const double Bdot=card.fractional_drive_per_year*tangent->B0Count()/Year;
    require(Bdot==-2.413652026364133e37,"native P2 Bdot changed");
    const std::string fate_id="phase6a1-P2-terminal-fate-v1";
    BNV::ProductFateLedger fate(fate_id,{{"controlled-terminal-product",BNV::TerminalProductFate::BoundInert,1,"abstract-neutron-disappearance"}});
    auto partition=Partition(card,fixture,bnv_token,fate_id);
    auto history=std::make_shared<const Phase6A1Test::UniformProperNeutronSinkHistory>(fixture.central,tangent->B0Count(),Bdot,
      tangent->DomainIdentity(),tangent->StarIdentity(),tangent->SequenceStateIdentity(),partition->Identity(),fate_id,bnv_token);
    const std::string potential="frozen Structure-1 B0 equilibrium mu_B plus governed moving-reference actual-potential correction";
    auto source=std::make_shared<const BNV::FrozenControlledBnvRunContext>(ordinary,tangent,history,partition,fate,monitor,channels,mu_B,potential,"RF-BNV-P2-NATIVE-BOUNDED-v1-SOURCE");
    auto zero_history=std::make_shared<const Campaign::ExactZeroHistory>(tangent->B0Count(),tangent->DomainIdentity(),
      tangent->StarIdentity(),tangent->SequenceStateIdentity(),bnv_token);
    auto zero_partition=std::make_shared<const Campaign::ExactZeroPartition>(bnv_token);
    BNV::ProductFateLedger zero_fate("zero-control-fate",{{"none",BNV::TerminalProductFate::BoundInert,1,"zero-control"}});
    auto control=std::make_shared<const BNV::FrozenControlledBnvRunContext>(ordinary,tangent,zero_history,zero_partition,zero_fate,
      monitor,channels,mu_B,potential,"RF-BNV-P2-NATIVE-BOUNDED-v1-CONTROL");
    return {source,control};
}
