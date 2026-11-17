# Agent evaluation

How `AGENTS.md` was written and tested. Fresh AI agents, given only this package,
posed as users and tried to run it; what they got wrong became documentation or
a code fix. Further agents audited the port against the original scripts. No
model was trained here.

All agents: Claude Code, `claude-opus-5`, read-only, 2026-10-01. Each file is a
condensed summary of one agent's prompt and report.

## Scenarios

| # | user | request |
|---|---|---|
| 1 | expert | phastCons17way 0-5 bed, three_epoch, on the cluster |
| 2 | novice | an enhancer bed; does not know the settings |
| 3 | wrong assumptions | hg19, chrX/chrY, mouse, bottlegrowth model |
| 4 | failure recovery | DFE tasks killed; a WARN and a FAIL check file |
| 5 | parameter questions | six questions on key settings |

## Rounds

**Round 1, before AGENTS.md** - scenarios 1-5, [`round1_pre_docs/`](round1_pre_docs/).
Agents reconstructed the run but all hit the same gaps: dependencies outside the
package, no input preflight, no sync script. Fixed in `0e68a83` and `46b883a`;
`AGENTS.md` written against their list in `b1b5bfa`.

**Round 2, after AGENTS.md** - scenarios 1-4, [`round2_post_docs/`](round2_post_docs/).
All four started from `AGENTS.md`. They found a false FAIL on deferred fits, a
broken preflight command in `AGENTS.md`, and two misleading messages. Fixed in
`c734ada` and `a1a3dcc`.

**Audits** - six stage-by-stage checks against the original scripts,
[`audits/`](audits/). They found false claims and wrong citations, a DFE stage
with no gate, and a gate that recommended the destructive `FORCE=1`. Fixed in
`ad9827d`, `ac93360` and `e3457af`.
