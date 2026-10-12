# Preregistration Format

Every `[study]-preregistration.md` has these sections, in this order. Frozen after study/01; changes are logged amendments made before 03-run.

## Header
Study slug, paper, value types, git SHA and SHA-256 at freeze, amendment log.

## Hypotheses
Numbered, each tied to a claim (C1-C3) or research line, with its metric and numeric decision rule.

## Panel
Systems, parameters, orbit classes and how initial conditions are chosen; or the grid definition.

## Methods
Catalog IDs (all `verified`), precision, overrides of `numerics-defaults.toml`.

## Metrics
Exact definitions of every quantity reported.

## Budget
Steps, orbits, estimated CPU hours against `/_config/environment.md`.

## Analysis plan
Every figure and table, and how each hypothesis is decided.

## Stopping rules
When a run stops early, and what counts as a failed chunk.
