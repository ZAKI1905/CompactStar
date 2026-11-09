# Proposed bounded thermal qualification — separate authorization required

Native assembly now supplies the missing complete owners. The previous
native-assembly/NONLINEAR_PLAN.md remains the fixed accuracy specification,
including every norm, uncertainty allocation and stopping rule. Its former
missing-input/launch status is superseded by this package, not its numerical
criteria. BA12 and BA12R remain FAIL. No trajectory has been executed here.

The compiled `native_p2_stage` collects actual-path evidence using the existing
UninterruptedBnvTrajectory, BreakpointStepper, passive schedules, isolated
Rk8pdCheckpointReconstructor, and ComputeCheckpointR20. This collection is not
an automatic accuracy certificate. In particular, independent accumulated
roundoff, U_main, reference interpolation uncertainty, event uncertainty, and
matched source-minus-control signal resolution must be assessed on the actual
outputs. Missing terms remain UNQUALIFIED, never zero. No static RHS evaluation
can establish them. See tests/native/stage.cpp:1 and the earlier plan:82.

Stage A: both roles start at original (t,x,eta_e,eta_mu)=(0,0,0,0), end
31,557,600 s (1 Julian year). Stage B: both roles independently restart from the
same original initial conditions, end 64,629,964,800 s (2048 Julian years).
Stage B requires Stage A's capped-policy/reference assessment to pass and owner
authorization; it is never submitted automatically. No extrapolation to 14,648
or 500,000 years. No automatic extension to find a breakpoint.

Per stage: 24 RKF45 candidate trajectories (two roles, three tolerance tiers,
four caps), ten independent references, 24 endpoint-only schedule repeats,
and six restart tails. Three rtol/atol tiers are 1e-7/(1e-12,1e-18,1e-18),
1e-9/(1e-14,1e-20,1e-20), 1e-11/(1e-16,1e-22,1e-22), initial step 1 s.
Caps: infinity, 38,723,807,757, 19,361,903,878.5, 9,680,951,939.25 s.
References: rk8pd 1e-14/(1e-19,1e-25,1e-25) at H=duration/256;
rk8pd 1e-15/(1e-20,1e-26,1e-26) at H/2,H/4; independent RKF45 at
that tighter tolerance and H/4,H/8. Every main/reference begins at original t0.
Passive scientific grids contain 240,480,960,1920 intervals. Preserve all
accepted states, event records, reconstruction error contributions, diagnostics,
and original R20 grid trapezoids. Quarter-cap restart tails start at the first
accepted endpoint at or beyond half-time, retaining the original prefix and
role/time/state. They do not promise bitwise continuation of adaptive step history.

Assessment retains the earlier plan exactly: normalized state error including
reference/output uncertainty <=1; reference allocation <=0.05 of ULTRA budget;
reconstruction <=0.20 Fi; powers <=1e-7 gross power (BASELINE), <=1e-9
(REFINED/ULTRA); R18/gross power <=1e-10. R20 including U_main,
U_reconstruction, U_event and U_quad <=2e-4 N, reconstruction/event <=5e-6 N.
Nested quadrature and references must contract at <=0.5 above independently
bounded roundoff. Use at least64 ulps without treating that floor as a bound on
accumulated roundoff. Add source/control uncertainties; signal resolution needs
U_delta<=0.10|delta_ref|. Any unknown uncertainty prevents qualification.

The collector stops at the first refusal or unresolved checkpoint, preserving
already-written artifacts; it does not retry, retune or automatically continue
the remaining matrix. Require rejection <=10%, finite progress, final endpoint,
Tinf within [1e6,1e9] K, historical and conditional frozen utilization <=1,
|DeltaB/B0|<=1e-6 and unchanged scientific guards. A timeout is INCOMPLETE.
No actual crossing means event coverage remains UNQUALIFIED. An explicit
subsequent decision is needed for any missing/failing matrix members.

Resource ceiling per stage: phyall/phyn004, one task/CPU, 4 GiB, aggregate30min;
single-thread libraries. Assembly measured about109s/522MiB in initial static
verification; trajectory runtime is unmeasured. Build separately4CPU/8GiB/15min.
All data remain under /mnt/sdd. Exact commands and snapshot hashes are recorded
in REPRODUCE.md and snapshots/identities.json after final native verification.
