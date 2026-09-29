# CompactStar plotting removal ledger

Authority: canonical `812463ac9ed374f64ac9cadd500066ab723d3a6c`, owner implementation and resume authorization.
75 active Zaki/CONFIND plot calls removed from nine compiled translation units;
zero remain. Twelve new table export sites retain previously unexported derived
data. `PlotFermiE` uses its existing `ExportFermiE` companion. All arithmetic,
MakeSmooth calls, contour ordering, numerical exports and TaskManager threading
remain in place. Visualization belongs in external tools.

Classification B: data already retained/exported or returned to caller; C: useful
derived data now exported. Every listed plotting operation was read-only with
respect to numerical dataset/curve contents. Setup only changed removed plotting
state. No numerical-side-effect plot deletions were found. The arm64 CONFIND
Plot implementation was a no-op. Source positions below refer to the canonical
pre-migration source. The JSON companion retains exact actions and replacements.

Two unbuilt legacy demos (`Table_5-8_Glenn.cpp`, `rotating_ns.cpp`) also lose their
embedded matplotlib branches (12 calls); their existing sequence/profile TSV
outputs remain the visualization inputs. Historical vendored headers are preserved.

| File:line | Function | API | Class/action | Data retention / replacement |
|---|---|---|---|---|
| `CompactStar/Core/src/TaskManager.cpp:347` | `TaskManager::FindCriticalCurve` | `con.Plot` | B / DELETE | Existing critical/intersection/contour exports; historical arm64 CONFIND Plot is a no-op.  |
| `CompactStar/Core/src/TaskManager.cpp:415` | `TaskManager::FindCriticalCurve` | `critical_curve.Plot` | B / DELETE | Existing critical/intersection/contour exports; historical arm64 CONFIND Plot is a no-op.  |
| `CompactStar/Core/src/TaskManager.cpp:519` | `TaskManager::FindMtotContour` | `con.Plot` | B / DELETE | Existing critical/intersection/contour exports; historical arm64 CONFIND Plot is a no-op.  |
| `CompactStar/Core/src/TaskManager.cpp:623` | `TaskManager::FindBtotContour` | `Zaki::Math::Curve2D::Plot` | B / DELETE | Existing critical/intersection/contour exports; historical arm64 CONFIND Plot is a no-op.  |
| `CompactStar/Core/src/TaskManager.cpp:628` | `TaskManager::FindBtotContour` | `con.Plot` | B / DELETE | Existing critical/intersection/contour exports; historical arm64 CONFIND Plot is a no-op.  |
| `CompactStar/Extensions/MixedStar/src/DarkCore_Analysis.cpp:373` | `DarkCore_Analysis::ExportBNV` | `neutron_out.Plot` | B / DELETE | BNV_tau species export retains inputs; gamma curves are fixed rescalings of its columns.  |
| `CompactStar/Extensions/MixedStar/src/DarkCore_Analysis.cpp:374` | `DarkCore_Analysis::ExportBNV` | `neutron_out.Plot` | B / DELETE | BNV_tau species export retains inputs; gamma curves are fixed rescalings of its columns.  |
| `CompactStar/Extensions/MixedStar/src/DarkCore_Analysis.cpp:382` | `DarkCore_Analysis::ExportBNV` | `n_gamma.Plot` | B / DELETE | BNV_tau species export retains inputs; gamma curves are fixed rescalings of its columns.  |
| `CompactStar/Extensions/MixedStar/src/DarkCore_Analysis.cpp:383` | `DarkCore_Analysis::ExportBNV` | `n_gamma.Plot` | B / DELETE | BNV_tau species export retains inputs; gamma curves are fixed rescalings of its columns.  |
| `CompactStar/Extensions/MixedStar/src/DarkCore_Analysis.cpp:384` | `DarkCore_Analysis::ExportBNV` | `n_gamma.Plot` | B / DELETE | BNV_tau species export retains inputs; gamma curves are fixed rescalings of its columns.  |
| `CompactStar/Extensions/MixedStar/src/DarkCore_Analysis.cpp:386` | `DarkCore_Analysis::ExportBNV` | `n_gamma.SemiLogYPlot` | B / DELETE | BNV_tau species export retains inputs; gamma curves are fixed rescalings of its columns.  |
| `CompactStar/Extensions/MixedStar/src/DarkCore_Analysis.cpp:395` | `DarkCore_Analysis::ExportBNV` | `lambda_out.Plot` | B / DELETE | BNV_tau species export retains inputs; gamma curves are fixed rescalings of its columns.  |
| `CompactStar/Extensions/MixedStar/src/DarkCore_Analysis.cpp:396` | `DarkCore_Analysis::ExportBNV` | `lambda_out.Plot` | B / DELETE | BNV_tau species export retains inputs; gamma curves are fixed rescalings of its columns.  |
| `CompactStar/Extensions/MixedStar/src/DarkCore_Analysis.cpp:402` | `DarkCore_Analysis::ExportBNV` | `lambda_gamma.Plot` | B / DELETE | BNV_tau species export retains inputs; gamma curves are fixed rescalings of its columns.  |
| `CompactStar/Extensions/MixedStar/src/DarkCore_Analysis.cpp:403` | `DarkCore_Analysis::ExportBNV` | `lambda_gamma.Plot` | B / DELETE | BNV_tau species export retains inputs; gamma curves are fixed rescalings of its columns.  |
| `CompactStar/Extensions/MixedStar/src/DarkCore_Analysis.cpp:404` | `DarkCore_Analysis::ExportBNV` | `lambda_gamma.Plot` | B / DELETE | BNV_tau species export retains inputs; gamma curves are fixed rescalings of its columns.  |
| `CompactStar/Extensions/MixedStar/src/DarkCore_Analysis.cpp:406` | `DarkCore_Analysis::ExportBNV` | `lambda_gamma.SemiLogYPlot` | B / DELETE | BNV_tau species export retains inputs; gamma curves are fixed rescalings of its columns.  |
| `CompactStar/Extensions/MixedStar/src/DarkCore_Analysis.cpp:415` | `DarkCore_Analysis::ExportBNV` | `sigmam_out.Plot` | B / DELETE | BNV_tau species export retains inputs; gamma curves are fixed rescalings of its columns.  |
| `CompactStar/Extensions/MixedStar/src/DarkCore_Analysis.cpp:416` | `DarkCore_Analysis::ExportBNV` | `sigmam_out.Plot` | B / DELETE | BNV_tau species export retains inputs; gamma curves are fixed rescalings of its columns.  |
| `CompactStar/Extensions/MixedStar/src/DarkCore_Analysis.cpp:422` | `DarkCore_Analysis::ExportBNV` | `sigmam_gamma.Plot` | B / DELETE | BNV_tau species export retains inputs; gamma curves are fixed rescalings of its columns.  |
| `CompactStar/Extensions/MixedStar/src/DarkCore_Analysis.cpp:423` | `DarkCore_Analysis::ExportBNV` | `sigmam_gamma.Plot` | B / DELETE | BNV_tau species export retains inputs; gamma curves are fixed rescalings of its columns.  |
| `CompactStar/Extensions/MixedStar/src/DarkCore_Analysis.cpp:424` | `DarkCore_Analysis::ExportBNV` | `sigmam_gamma.Plot` | B / DELETE | BNV_tau species export retains inputs; gamma curves are fixed rescalings of its columns.  |
| `CompactStar/Extensions/MixedStar/src/DarkCore_Analysis.cpp:426` | `DarkCore_Analysis::ExportBNV` | `sigmam_gamma.SemiLogYPlot` | B / DELETE | BNV_tau species export retains inputs; gamma curves are fixed rescalings of its columns.  |
| `CompactStar/Microphysics/BNV/Analysis/src/BNV_Sequence.cpp:441` | `MicroBNVAna::BNV_Sequence::Find_b_factors` | `b_o_ds.Plot` | B / DELETE | Existing B_Factors/evolution export or shared derived dataset retained by one table export.  |
| `CompactStar/Microphysics/BNV/Analysis/src/BNV_Sequence.cpp:634` | `MicroBNVAna::BNV_Sequence::Plot_Dimless_O` | `Scaled_O.SemiLogXPlot` | B / DELETE | Existing B_Factors/evolution export or shared derived dataset retained by one table export.  |
| `CompactStar/Microphysics/BNV/Analysis/src/BNV_Sequence.cpp:641` | `MicroBNVAna::BNV_Sequence::Plot_Dimless_O` | `Scaled_O.Plot` | C / EXPORT_DATA | New deterministic table; calculations retained verbatim. `Scaled_O.Export("Dimensionless_Observables.tsv");` |
| `CompactStar/Microphysics/BNV/Analysis/src/BNV_Sequence.cpp:690` | `MicroBNVAna::BNV_Sequence::EvalBeta` | `tmp_ds.SemiLogYPlot` | C / EXPORT_DATA | New deterministic table; calculations retained verbatim. `tmp_ds.Export("Derived_B_Factors.tsv");` |
| `CompactStar/Microphysics/BNV/Analysis/src/BNV_Sequence.cpp:1501` | `MicroBNVAna::BNV_Sequence::Solve` | `omega_ds.SemiLogXPlot` | C / EXPORT_DATA | New deterministic table; calculations retained verbatim. `omega_ds.Export(omega_t_evol_f_name + file_stamps_str + ".tsv");` |
| `CompactStar/Microphysics/BNV/Analysis/src/BNV_Sequence.cpp:1516` | `MicroBNVAna::BNV_Sequence::Solve` | `br_idx_ds.SemiLogXPlot` | C / EXPORT_DATA | New deterministic table; calculations retained verbatim. `br_idx_ds.Export(br_idx_evol_f_name + file_stamps_str + ".tsv");` |
| `CompactStar/Microphysics/BNV/Analysis/src/BNV_Sequence.cpp:1536` | `MicroBNVAna::BNV_Sequence::Solve` | `P_Pdot_ds.LogLogPlot` | B / DELETE | Existing B_Factors/evolution export or shared derived dataset retained by one table export.  |
| `CompactStar/Microphysics/BNV/Analysis/src/BNV_Sequence.cpp:1569` | `MicroBNVAna::BNV_Sequence::Solve` | `M_Omega_ds.SemiLogYPlot` | B / DELETE | Existing B_Factors/evolution export or shared derived dataset retained by one table export.  |
| `CompactStar/Microphysics/BNV/Analysis/src/Decay_Analysis.cpp:160` | `MicroBNVAna::Decay_Analysis::ImportEffMass` | `m_eff_ds.Plot` | B / DELETE | Imported microphysics input or existing Lambda/Neutron Plot Points.tsv export.  |
| `CompactStar/Microphysics/BNV/Analysis/src/Decay_Analysis.cpp:186` | `MicroBNVAna::Decay_Analysis::ImportVSelfEnergy` | `V_self_E_ds.Plot` | B / DELETE | Imported microphysics input or existing Lambda/Neutron Plot Points.tsv export.  |
| `CompactStar/Microphysics/BNV/Analysis/src/Decay_Analysis.cpp:550` | `MicroBNVAna::Decay_Analysis::AttachPulsar` | `lam_plt_pt.Plot` | B / DELETE | Imported microphysics input or existing Lambda/Neutron Plot Points.tsv export.  |
| `CompactStar/Microphysics/BNV/Analysis/src/Decay_Analysis.cpp:552` | `MicroBNVAna::Decay_Analysis::AttachPulsar` | `lam_plt_pt.SemiLogYPlot` | B / DELETE | Imported microphysics input or existing Lambda/Neutron Plot Points.tsv export.  |
| `CompactStar/Microphysics/BNV/Analysis/src/Decay_Analysis.cpp:558` | `MicroBNVAna::Decay_Analysis::AttachPulsar` | `neu_plt_pt.Plot` | B / DELETE | Imported microphysics input or existing Lambda/Neutron Plot Points.tsv export.  |
| `CompactStar/Microphysics/BNV/Analysis/src/Decay_Analysis.cpp:560` | `MicroBNVAna::Decay_Analysis::AttachPulsar` | `neu_plt_pt.SemiLogYPlot` | B / DELETE | Imported microphysics input or existing Lambda/Neutron Plot Points.tsv export.  |
| `CompactStar/Microphysics/BNV/Analysis/src/BNV_Analysis.cpp:273` | `MicroBNVAna::BNV_Analysis::Evolve` | `n_out.SemiLogXPlot` | B / DELETE | Existing n/Lambda/sigma evolution.tsv export.  |
| `CompactStar/Microphysics/BNV/Analysis/src/BNV_Analysis.cpp:318` | `MicroBNVAna::BNV_Analysis::Evolve` | `lam_out.SemiLogXPlot` | B / DELETE | Existing n/Lambda/sigma evolution.tsv export.  |
| `CompactStar/Microphysics/BNV/Analysis/src/BNV_Analysis.cpp:363` | `MicroBNVAna::BNV_Analysis::Evolve` | `sig_out.SemiLogXPlot` | B / DELETE | Existing n/Lambda/sigma evolution.tsv export.  |
| `CompactStar/Microphysics/BNV/Internal/src/BNV_Chi.cpp:542` | `MicroBNVInt::BNV_Chi::Plot_Meff_Radius` | `ds_meff_r.Plot` | C / EXPORT_DATA | New deterministic table; calculations retained verbatim. `ds_meff_r.Export(pulsar.GetName() + "_m_vs_R.tsv");` |
| `CompactStar/Microphysics/BNV/Internal/src/BNV_Chi.cpp:658` | `MicroBNVInt::BNV_Chi::Plot_RestEnergy_Radius` | `ds_E_r.Plot` | C / EXPORT_DATA | New deterministic table; calculations retained verbatim. `ds_E_r.Export(pulsar.GetName() + "_E_vs_R.tsv");` |
| `CompactStar/Microphysics/BNV/Internal/src/BNV_Chi.cpp:804` | `MicroBNVInt::BNV_Chi::Plot_EF_Radius` | `ds_EF_r.Plot` | C / EXPORT_DATA | New deterministic table; calculations retained verbatim. `ds_EF_r.Export(pulsar.GetName() + "_EF_vs_R.tsv");` |
| `CompactStar/Microphysics/BNV/Internal/src/BNV_Chi.cpp:879` | `MicroBNVInt::BNV_Chi::Plot_Estar_Radius` | `ds_E0_EF_r.Plot` | C / EXPORT_DATA | New deterministic table; calculations retained verbatim. `ds_E0_EF_r.Export(pulsar.GetName() + "_Estar_vs_R_" + B.short_name + ".tsv");` |
| `CompactStar/Microphysics/BNV/Internal/src/BNV_Chi.cpp:975` | `MicroBNVInt::BNV_Chi::Plot_CM_E_Radius` | `ds_E0_EF_r.Plot` | B / DELETE | Existing rate/epsilon/radial table export following the removed plot.  |
| `CompactStar/Microphysics/BNV/Internal/src/BNV_Chi.cpp:1102` | `MicroBNVInt::BNV_Chi::Plot_RestE_EF_Radius` | `ds_E0_EF_r.Plot` | C / EXPORT_DATA | New deterministic table; calculations retained verbatim. `ds_E0_EF_r.Export(pulsar.GetName() + "_EBand_vs_R_" + B.short_name + ".tsv");` |
| `CompactStar/Microphysics/BNV/Internal/src/BNV_Chi.cpp:1229` | `MicroBNVInt::BNV_Chi::PlotVacuumBrLim` | `dec_lim_ds.SemiLogYPlot` | C / EXPORT_DATA | New deterministic table; calculations retained verbatim. `dec_lim_ds.Export("Br(" + B.short_name + ")_vs_m_chi.tsv");` |
| `CompactStar/Microphysics/BNV/Internal/src/BNV_Chi.cpp:1311` | `MicroBNVInt::BNV_Chi::PlotRate_Eps` | `rate.SemiLogYPlot` | B / DELETE | Existing rate/epsilon/radial table export following the removed plot.  |
| `CompactStar/Microphysics/BNV/Internal/src/BNV_Chi.cpp:1318` | `MicroBNVInt::BNV_Chi::PlotRate_Eps` | `rate.SemiLogYPlot` | B / DELETE | Existing rate/epsilon/radial table export following the removed plot.  |
| `CompactStar/Microphysics/BNV/Internal/src/BNV_Chi.cpp:1392` | `MicroBNVInt::BNV_Chi::PlotRate_Eps` | `rate.SemiLogYPlot` | B / DELETE | Existing rate/epsilon/radial table export following the removed plot.  |
| `CompactStar/Microphysics/BNV/Internal/src/BNV_Chi.cpp:1402` | `MicroBNVInt::BNV_Chi::PlotRate_Eps` | `rate.SemiLogYPlot` | B / DELETE | Existing rate/epsilon/radial table export following the removed plot.  |
| `CompactStar/Microphysics/BNV/Internal/src/BNV_Chi.cpp:1577` | `MicroBNVInt::BNV_Chi::Rate_vs_R` | `rate_vs_r.SemiLogYPlot` | B / DELETE | Existing rate/epsilon/radial table export following the removed plot.  |
| `CompactStar/Microphysics/BNV/Internal/src/BNV_Chi.cpp:1753` | `MicroBNVInt::BNV_Chi::hidden_Plot_Rate_vs_R` | `ds->SemiLogYPlot` | B / DELETE | Dataset remains returned to the numerical caller; deletion changes visualization only.  |
| `CompactStar/Microphysics/BNV/Internal/src/BNV_Chi.cpp:1783` | `MicroBNVInt::BNV_Chi::hidden_Plot_Rate_vs_Density` | `ds->SemiLogYPlot` | B / DELETE | Dataset remains returned to the numerical caller; deletion changes visualization only.  |
| `CompactStar/Microphysics/BNV/Channels/src/BNV_B_Chi_Photon.cpp:705` | `MicroBNVCh::BNV_B_Chi_Photon::hidden_Plot_Thermal_Hole_E_Rate_vs_Density` | `ds->SemiLogYPlot` | B / DELETE | Dataset remains returned to the numerical caller; deletion changes visualization only.  |
| `CompactStar/Microphysics/BNV/Channels/src/BNV_B_Chi_Photon.cpp:837` | `MicroBNVCh::BNV_B_Chi_Photon::Thermal_Hole_E_Rate_vs_R` | `rate_vs_r.SemiLogYPlot` | B / DELETE | Existing rate/epsilon/radial table export following the removed plot.  |
| `CompactStar/Microphysics/BNV/Channels/src/BNV_B_Chi_Photon.cpp:909` | `MicroBNVCh::BNV_B_Chi_Photon::hidden_Plot_Thermal_Hole_E_Rate_vs_R` | `ds->SemiLogYPlot` | B / DELETE | Dataset remains returned to the numerical caller; deletion changes visualization only.  |
| `CompactStar/Microphysics/BNV/Channels/src/BNV_B_Chi_Photon.cpp:959` | `MicroBNVCh::BNV_B_Chi_Photon::hidden_Plot_Thermal_Photon_E_Rate_vs_R` | `ds->SemiLogYPlot` | B / DELETE | Dataset remains returned to the numerical caller; deletion changes visualization only.  |
| `CompactStar/Microphysics/BNV/Channels/src/BNV_B_Chi_Photon.cpp:994` | `MicroBNVCh::BNV_B_Chi_Photon::hidden_Plot_Thermal_Photon_E_Rate_vs_Density` | `ds->SemiLogYPlot` | B / DELETE | Dataset remains returned to the numerical caller; deletion changes visualization only.  |
| `CompactStar/Microphysics/BNV/Channels/src/BNV_B_Chi_Photon.cpp:1169` | `MicroBNVCh::BNV_B_Chi_Photon::Thermal_Photon_E_Rate_vs_R` | `rate_vs_r.SemiLogYPlot` | B / DELETE | Existing rate/epsilon/radial table export following the removed plot.  |
| `CompactStar/Microphysics/BNV/Channels/src/BNV_B_Chi_Photon.cpp:1277` | `MicroBNVCh::BNV_B_Chi_Photon::Plot_Thermal_Photon_E_Rate` | `rate.SemiLogYPlot` | B / DELETE | Existing rate/epsilon/radial table export following the removed plot.  |
| `CompactStar/Microphysics/BNV/Channels/src/BNV_B_Chi_Photon.cpp:1335` | `MicroBNVCh::BNV_B_Chi_Photon::Plot_Limited_Thermal_E_Rate` | `default_rate.SemiLogYPlot` | B / DELETE | Existing rate/epsilon/radial table export following the removed plot.  |
| `CompactStar/Microphysics/BNV/Channels/src/BNV_B_Chi_Photon.cpp:1819` | `MicroBNVCh::BNV_B_Chi_Photon::Thermal_Total_E_Rate_vs_R` | `rate_vs_r.SemiLogYPlot` | B / DELETE | Existing rate/epsilon/radial table export following the removed plot.  |
| `CompactStar/Microphysics/BNV/Channels/src/BNV_B_Chi_Photon.cpp:1926` | `MicroBNVCh::BNV_B_Chi_Photon::Plot_Thermal_Total_E_Rate` | `rate.SemiLogYPlot` | B / DELETE | Existing rate/epsilon/radial table export following the removed plot.  |
| `CompactStar/Microphysics/BNV/Channels/src/BNV_B_Chi_Transition.cpp:429` | `MicroBNVCh::BNV_B_Chi_Transition::PlotTransCond` | `ds_Cond_a.Plot` | C / EXPORT_DATA | New deterministic table; calculations retained verbatim. `ds_Cond_a.Export(model + "/" + process.name + "/Transit_Lim_a_" + model + "_" + B.short_name + "_" + m_chi_str + ".tsv");` |
| `CompactStar/Microphysics/BNV/Channels/src/BNV_B_Chi_Transition.cpp:506` | `MicroBNVCh::BNV_B_Chi_Transition::PlotTransBand` | `ds_trans_range.Plot` | C / EXPORT_DATA | New deterministic table; calculations retained verbatim. `ds_trans_range.Export(model + "/" + process.name + "/Transit_Den_Ranges_" + model + "_" + B.short_name + ".tsv");` |
| `CompactStar/EOS/src/CompOSE_EOS.cpp:332` | `CompactStar::CompOSE_EOS::ImportThermo` | `eos.LogLogPlot` | B / DELETE | Existing eos/micro exports and ExportFermiE; PlotFermiE delegates to the existing numerical exporter.  |
| `CompactStar/EOS/src/CompOSE_EOS.cpp:458` | `CompactStar::CompOSE_EOS::ImportCompo` | `eos.LogLogPlot` | B / DELETE | Existing eos/micro exports and ExportFermiE; PlotFermiE delegates to the existing numerical exporter.  |
| `CompactStar/EOS/src/CompOSE_EOS.cpp:757` | `CompactStar::CompOSE_EOS::ImportMicro` | `m_eff.Plot` | B / DELETE | Existing eos/micro exports and ExportFermiE; PlotFermiE delegates to the existing numerical exporter.  |
| `CompactStar/EOS/src/CompOSE_EOS.cpp:786` | `CompactStar::CompOSE_EOS::ImportMicro` | `V_eff.LogLogPlot` | B / DELETE | Existing eos/micro exports and ExportFermiE; PlotFermiE delegates to the existing numerical exporter.  |
| `CompactStar/EOS/src/CompOSE_EOS.cpp:815` | `CompactStar::CompOSE_EOS::ImportMicro` | `U.Plot` | B / DELETE | Existing eos/micro exports and ExportFermiE; PlotFermiE delegates to the existing numerical exporter.  |
| `CompactStar/EOS/src/CompOSE_EOS.cpp:921` | `CompactStar::CompOSE_EOS::ExtendMeffToCrust` | `m_eff.Plot` | B / DELETE | Existing eos/micro exports and ExportFermiE; PlotFermiE delegates to the existing numerical exporter.  |
| `CompactStar/EOS/src/CompOSE_EOS.cpp:1020` | `CompactStar::CompOSE_EOS::ExtendVeffToCrust` | `V_eff.Plot` | B / DELETE | Existing eos/micro exports and ExportFermiE; PlotFermiE delegates to the existing numerical exporter.  |
| `CompactStar/EOS/src/CompOSE_EOS.cpp:1122` | `CompactStar::CompOSE_EOS::ExtendUeffToCrust` | `U.Plot` | B / DELETE | Existing eos/micro exports and ExportFermiE; PlotFermiE delegates to the existing numerical exporter.  |
| `CompactStar/EOS/src/CompOSE_EOS.cpp:1279` | `CompactStar::CompOSE_EOS::PlotFermiE` | `fermi_ds.SemiLogXPlot` | B / DELETE | Existing eos/micro exports and ExportFermiE; PlotFermiE delegates to the existing numerical exporter.  |

Numerical side-effect assessment for every row: no numerical mutation in the removed call; subsequent arithmetic preserved. Exact TaskManager and governed-suite results are recorded in the migration qualification report. Unexercised legacy presentation paths have source-level rather than runtime evidence.

The final cleanup also removes five CompOSE legend-selection branches whose
only result was the deleted plot-index vector. Their Max/Min tests filtered
legend entries only and never wrote EOS data. The historical ImportThermo
working-directory assignment is retained. Existing public boolean arguments
and method names remain source-compatible; they do not recreate plotting.
