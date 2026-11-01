# Panel in rounds (Amendment D94)

For a study whose panel block uses the keys of experiments Amendment D94. The keys, their meaning and their checks have one home: `workspaces/module-build/stages/01-spec/output/experiments-spec.md`, "Amendment D94", R22 (panel input) and R25 (run block). This note only says what changes at this stage.

- New panel keys: `screen_source`, `screen_points`, `r_max`, `gate_rounds`, `quota`, `d_fraction`, `db_size`, `newton_grid_side`, `controls`. Values come from the preregistration only.
- All `r_max` rounds are drawn at build time into one panel file (P70); which rounds run is decided later by the round log at stage 03.
- The build also writes `<panel>.pools.toml` (class bits per screening point id, no screening value) for the oracle; keep it next to the panel file.
- Checkpoint: a round cell may be `short` or `empty` when its pool is spent; it is never filled from another class or round. Read the per-cell status in the panel manifest (`round_cells`) together with the portraits.
- The `controls = "none"` panel holds no control orbit, so the audit check "Classes" covers the classes of the preregistration only.
