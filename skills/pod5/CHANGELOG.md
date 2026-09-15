# Changelog — pod5

## [Unreleased]

## [2026-03-30]
- Skill version: 1.1.0
- Fixed: corrected binary path from `bin/bin/pod5` to `bin/pod5` (double-bin was a typo)
- Updated: re-installed conda environment (`rm -rf venvs/pod5 && bash install.sh`); tool version remains 0.3.23
- Changed: Environment section now uses allowlist-safe `readlink -f` probe pattern instead of `$VAR` shell variable export
- Added: `## Allowlist entries` section with per-machine resolution instructions
- Changed: `## Full flag reference` and `## Patterns` grep paths updated to use skill-relative `$(dirname ...)` form

## [2026-03-18] — Initial skill
- Tool version: 0.3.23 (run `readlink -f "$(dirname <SKILL.md path>)/bin/pod5"` then `<path> --version` to check)
- Skill version: 1.0.0
- Added: initial SKILL.md, reference.md, patterns.md, environment.yaml
