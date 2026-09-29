# Latent undefined spin diagnostics

Removing Python3::Module's `-undefined dynamic_lookup` made the existing
`psr_debug_main` link fail on `Physics::Spin::CharacteristicAge` and
`Physics::Spin::DipoleFieldEstimate`. Canonical 812463a declares these helpers
but contains no definitions. They were unresolved runtime traps, not qualified
numerical implementations. The dipole declaration explicitly calls its
normalization a placeholder. No governed fixture calls either helper.

The migration supplies explicit `logic_error` failures for these unavailable
entry points and tests that neither silently returns a fabricated value. No
spin equation, constant, normalization, or evolutionary driver is changed.
A future scientific implementation must qualify its contract separately.
This is a nonblocking pre-existing API limitation surfaced by strict linking.
