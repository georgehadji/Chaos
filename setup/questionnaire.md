# Onboarding Questionnaire

<!-- Agent instructions: read this file when the user types "setup". Ask ALL questions in one
     conversational pass; the user answers in one message. Replace each placeholder in the files
     listed. Then scan the repo (except _core/, workspace-builder/ and this file) for remaining {{ patterns
     and ask for anything missing. Setup is done only when none remain. -->

Answer what you know; skip any question to accept its default.

### Q1: Your name as it should appear on papers
- Placeholder: `{{AUTHOR_NAME}}`
- Files: `_config/publication.md`
- Type: free text
- Default: the git `user.name` of this repo

### Q2: Your affiliation (institution, department, city, country)
- Placeholder: `{{AFFILIATION}}`
- Files: `_config/publication.md`
- Type: free text
- Default: "Independent researcher"

### Q3: Your ORCID iD
- Placeholder: `{{ORCID}}`
- Files: `_config/publication.md`
- Type: free text
- Default: "none yet"

### Q4: Journals you would like to publish in, in order of preference
- Placeholder: `{{TARGET_JOURNALS}}`
- Files: `_config/publication.md`
- Type: free text (a list)
- Default: decided per paper at paper/01-outline
- Example: Physica D; Chaos (AIP); Communications in Nonlinear Science and Numerical Simulation; Celestial Mechanics and Dynamical Astronomy; Journal of Open Source Software (for the library)

### Q5: Citation style
- Placeholder: `{{CITATION_STYLE}}`
- Files: `_config/publication.md`
- Type: selection
- Options: the target journal's own style; author-year (natbib); numeric
- Default: the target journal's own style

### Q6: Where will long runs happen?
- Placeholder: `{{COMPUTE_TARGET}}`
- Files: `_config/environment.md`
- Type: selection
- Options: this computer; a workstation; a cluster (name it); cloud
- Default: this computer
- Derived: `{{MAX_CORES}}` in `_config/environment.md`. For this computer, the agent reads the logical core count (`python -c "import os; print(os.cpu_count())"`) and writes that number minus one; for other targets it asks for the core count once.

### Q7: Where should run data and paper full texts be stored?
- Placeholder: `{{DATA_ROOT}}`
- Files: `_config/environment.md`
- Type: free text (a folder path)
- Default: `data/` inside this repo (git-ignored). Pick a larger disk if you expect runs over about 50 GB.

### Q8: Extended-precision type for precision studies
- Placeholder: `{{HIGH_PRECISION_BACKEND}}`
- Files: `_config/environment.md`
- Type: selection
- Options: Double64 (DoubleFloats.jl, about 32 digits, fast); MultiFloats (Float64x2 to Float64x8, fast, adjustable); BigFloat only (slow, arbitrary digits)
- Default: Double64

### Q9: Will you write papers with this repo?
- Type: yes/no
- Default: yes
- If NO: remove `workspaces/paper/` and its rows in `CLAUDE.md` and `CONTEXT.md`; remove `_config/setup-guides/latex.md`
- If YES: keep it

---

## After Onboarding

Configured: author and publication settings (`_config/publication.md`), compute and storage (`_config/environment.md`). Next: milestone M0 in `_config/build-order.md`. Its first step is installing Julia (`_config/setup-guides/julia.md`).
